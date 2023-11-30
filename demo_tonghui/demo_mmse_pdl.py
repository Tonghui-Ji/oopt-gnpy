from pathlib import Path
import pytest
from numpy.testing import assert_allclose, assert_array_equal, assert_raises
from numpy import array
from copy import deepcopy

from gnpy.core.utils import lin2db, automatic_nch, dbm2watt, power_dbm_to_psd_mw_ghz, watt2dbm, psd2powerdbm,generate_random_numbers
from gnpy.core.network import build_network
from gnpy.core.elements import Roadm
from gnpy.core.info import create_input_spectral_information, Pref, create_arbitrary_spectral_information, \
    ReferenceCarrier
from gnpy.core.equipment import trx_mode_params
from gnpy.core.exceptions import ConfigurationError
from gnpy.tools.json_io import network_from_json, load_equipment, load_network, _spectrum_from_json, load_json, \
    Transceiver, requests_from_json
from gnpy.topology.request import PathRequest, compute_constrained_path, propagate, propagate_and_optimize_mode
from support.gobal_control import GlobalControl
from gnpy.core.parameters import Parameters
import numpy as np
import matplotlib.pyplot as plt
from gnpy.core.elements import Transceiver, Fiber, Edfa, Roadm,Fused
from gnpy.core.utils import db2lin,lin2db
from gnpy.core.mmse_pdl import pdl_pen
from networkx import dijkstra_path

gc = GlobalControl()
data_path = gc.gnpy_path/ 'data'
EQPT_FILENAME = data_path / 'eqpt_config.json'
NETWORK_FILENAME = data_path / 'testTopology_expected.json'
con_in = 1
con_out = 0
input_power = 0
dest = 'trx F'
baud_rate = 32e9
grid_hz = 50e9

def net_setup(equipment):
    """common setup for tests: builds network, equipment and oms only once"""
    network = load_network(NETWORK_FILENAME, equipment)
    spectrum = equipment['SI']['default']
    p_db = spectrum.power_dbm
    p_total_db = p_db + lin2db(automatic_nch(spectrum.f_min, spectrum.f_max, spectrum.spacing))
    build_network(network, equipment, p_db, p_total_db)
    return network

def create_voyager_req(equipment, source, dest, bidir, nodes_list, loose_list, mode, spacing, power_dbm):
    """create the usual request list according to parameters"""
    params = {'request_id': 'test_request',
              'source': source,
              'bidir': bidir,
              'destination': dest,
              'trx_type': 'Voyager',
              'trx_mode': mode,
              'format': mode,
              'spacing': spacing,
              'nodes_list': nodes_list,
              'loose_list': loose_list,
              'path_bandwidth': 100.0e9,
              'effective_freq_slot': None}
    trx_params = trx_mode_params(equipment, params['trx_type'], params['trx_mode'], True)
    params.update(trx_params)
    params['power'] = dbm2watt(power_dbm) if power_dbm else dbm2watt(equipment['SI']['default'].power_dbm)
    f_min = params['f_min']
    f_max_from_si = params['f_max']
    params['nb_channel'] = automatic_nch(f_min, f_max_from_si, params['spacing'])
    return PathRequest(**params)

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

if __name__ == '__main__':
    target = -8
    mode = 'mode 1'
    slot_width = 50e9
    deltap = 0
    equalization = 'target_psd_out_mWperGHz'
    equipment = load_equipment(EQPT_FILENAME)
    setattr(equipment['Roadm']['default'], 'target_pch_out_db', target)
    req = create_voyager_req(equipment, 'trx Brest_KLA', 'trx Rennes_STA', False,
                             ['east edfa in Brest_KLA to Quimper', 'roadm Lannion_CAS', 'trx Rennes_STA'],
                             ['STRICT', 'STRICT', 'STRICT'],
                             mode, slot_width, deltap)
    roadms = ['roadm Brest_KLA', 'roadm Lorient_KMA', 'roadm Lannion_CAS', 'roadm Rennes_STA']
    # degree = {'roadm Brest_KLA': 'east edfa in Brest_KLA to Quimper',
    #           'roadm Lorient_KMA': 'east edfa in Lorient_KMA to Loudeac'}
    # boosters = ['east edfa in Brest_KLA to Quimper', 'east edfa in Lorient_KMA to Loudeac',
    #             'east edfa in Lannion_CAS to Stbrieuc']
    target_psd = power_dbm_to_psd_mw_ghz(target, 32e9)
    ref = ReferenceCarrier(baud_rate=32e9, slot_width=50e9)
    
    delattr(equipment['Roadm']['default'], 'target_pch_out_db')
    setattr(equipment['Roadm']['default'], equalization, target_psd)
    network = net_setup(equipment)

    # for e in network.nodes():
    #     if isinstance(e, Fiber):
    #         loss = e.params.loss_coef * e.params.length
    #         e.params.con_in = con_in
    #         e.params.con_out = con_out
    #     if isinstance(e, Edfa):
    #         e.operational.gain_target = loss + con_in + con_out

    path = compute_constrained_path(network, req)
    si = create_input_spectral_information(
        f_min=req.f_min, f_max=req.f_max, roll_off=req.roll_off, baud_rate=req.baud_rate, power=req.power,
        spacing=req.spacing, tx_osnr=req.tx_osnr, ref_carrier=ref)
    # net_res = solve_network(path,si)

    # for i, el in enumerate(path):
    #     if isinstance(el, Roadm):
    #         si = el(si, degree=path[i + 1].uid)
    #     else:
    #         si = el(si)

    params = Parameters()
    params.snr_trx_dB = 17
    params.sig_power_dBm = 0
    params.Rs = baud_rate
    # FIXME: 根据si.frequency选择对应的channel
    params.Fc = 193e12

    link_params = Parameters()
    link_config = []
    wss_params_list = []
    oa_params_list = []
    fiu_params_list = []
    fiber_params_list = []
    voa_params_list = []
    fused_params_list = []
    for i, el in enumerate(path):
        if isinstance(el,Roadm):
            si = el(si, degree=path[i + 1].uid)
            link_config.append("wss")
            wss_params = Parameters()
            # wss_params.loss_dB = el.loss
            # wss_params.pdl_dB = el.params.pdl
            wss_params.loss_dB = 6
            wss_params.pdl = 0.5
            wss_params_list.append(wss_params)
        elif isinstance(el,Edfa):
            si = el(si)
            link_config.append("oa")
            oa_params = Parameters()
            oa_params.nf_dB = el.nf[0]
            oa_params.gain_dB = el.gprofile[0]
            # oa_params.pdl_dB = el.params.pdl
            oa_params.pdl_dB  = 0.3
            oa_params_list.append(oa_params)
        elif isinstance(el,Fiber):
            si = el(si)
            link_config.append("Fiber")
            fiber_params = Parameters()
            fiber_params.loss_dB = el.loss
            fiber_params_list.append(fiber_params)
        # elif isinstance(el,Fused):
        #     si = el(si)
        #     link_config.append("Fused")
        #     fused_params = Parameters()
        #     fused_params.loss_dB = el.loss
        #     fused_params_list.append(fused_params)

    link_params.link_config = link_config
    link_params.num_oa = link_config.count("oa")
    link_params.num_wss = link_config.count("wss")
    link_params.num_fiber = link_config.count("fiber")
    link_params.num_fiu = link_config.count("fiu")
    link_params.num_voa = link_config.count("voa")
    params.link_params = link_params
    params.oa_params_list = oa_params_list[::-1]
    params.wss_params_list = wss_params_list[::-1]
    params.fiu_params_list = fiu_params_list[::-1]
    params.fiber_params_list = fiber_params_list[::-1]
    params.voa_params_list = voa_params_list[::-1]
    params.alpha_list = generate_random_numbers('pdf',(lambda x:np.sin(2*x),0,np.pi/2),len(link_config))
    params.beta_list = generate_random_numbers('uniform',(0,np.pi*2),len(link_config))
    snr_esti_dB_x,snr_esti_dB_y = pdl_pen(params)

    pass