import sys
sys.path.append(r'C:\Users\44894\OneDrive\Coding\Python\oopt-gnpy')

from gnpy.core.utils import lin2db,db2lin, automatic_nch, power_dbm_to_psd_mw_ghz
from gnpy.core.info import create_input_spectral_information, ReferenceCarrier
from gnpy.core.network import build_network
from gnpy.core.elements import Roadm
from gnpy.core.parameters import Parameters
from gnpy.tools.json_io import load_equipment,network_from_params
from support.gobal_control import GlobalControl
import numpy as np
import matplotlib.pyplot as plt

gc = GlobalControl()
data_path = gc.gnpy_path/ 'data'
EQPT_FILENAME = data_path / 'eqpt_config.json'

def solve_network(path,si):
    net_res = Parameters()  # 仿真结果是一个parameters类型的，不知道这种是不是合理？
    sig_psd_list = []
    sig_power_list = []
    net_res.frequency = si.frequency
    for i, el in enumerate(path):
        sig_power_list.append(lin2db(1e3*np.sum(si.signal)))
        if isinstance(el, Roadm):
            sig_psd_list.append(lin2db(1e3*si.signal))
            si = el(si, degree=path[i + 1].uid)
            sig_psd_list.append(lin2db(1e3*si.signal))
        else:
            si = el(si)
        print(el.uid)
    sig_power_list.append(lin2db(1e3*np.sum(si.signal)))
    net_res.sig_psd = np.array(sig_psd_list).T
    net_res.sig_power_dB = np.array(sig_power_list)

    return net_res

def c120_l120_sim(input_power,target, equalization):
    """check that power target on roadm is correct for these cases; check on booster
    - SI : target_pch_out_db / target_psd_out_mWperGHz
    - node : target_pch_out_db / target_psd_out_mWperGHz
    - per degree : target_pch_out_db / target_psd_out_mWperGHz
    for these cases with and without power from user
    """
    equipment = load_equipment(EQPT_FILENAME)
    setattr(equipment['Roadm']['default'], 'target_pch_out_db', target)
    target_psd = power_dbm_to_psd_mw_ghz(target, 32e9)

    delattr(equipment['Roadm']['default'], 'target_pch_out_db')
    setattr(equipment['Roadm']['default'], equalization, target_psd)

    link_params = Parameters()
    link_params.spans_per_oms = 4
    link_params.span_num = 10
    el_params = Parameters()
    el_params.oa_params = Parameters()
    el_params.fiber_params = Parameters()
    el_params.wss_params = Parameters()

    el_params.oa_params.gain_target = 10
    el_params.fiber_params.length = 80

    network = network_from_params(link_params,equipment,el_params)
    spectrum = equipment['SI']['default']
    p_db = spectrum.power_dbm
    p_total_db = p_db + lin2db(automatic_nch(spectrum.f_min, spectrum.f_max, spectrum.spacing))
    build_network(network, equipment, p_db, p_total_db, no_insert_edfas=True)
    path = []

    path = [n for n in network.nodes]

    p = input_power
    p = db2lin(p) * 1e-3
    spacing = 50e9  # THz
    si = create_input_spectral_information(f_min=191.3e12, f_max=191.3e12 + 79 * spacing, roll_off=0.15,
                                           baud_rate=32e9, power=p, spacing=spacing, tx_osnr=None,
                                           ref_carrier=ReferenceCarrier(baud_rate=32e9, slot_width=50e9))
    net_res = solve_network(path,si)

    plt.figure()
    plt.plot(net_res.frequency,net_res.sig_psd)
    plt.legend()
    
    plt.figure()
    plt.plot(net_res.sig_power_dB,'-o')
    plt.show()

if __name__ == '__main__':
    c120_l120_sim(input_power=-5,target=-2,equalization='target_psd_out_mWperGHz')