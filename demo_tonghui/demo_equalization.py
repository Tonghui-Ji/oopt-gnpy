import sys
sys.path.append(r'C:\Users\44894\OneDrive\Coding\Python\oopt-gnpy')

from gnpy.core.utils import lin2db,db2lin, automatic_nch, dbm2watt, power_dbm_to_psd_mw_ghz, watt2dbm, psd2powerdbm
from gnpy.core.network import build_network
from gnpy.core.elements import Roadm
from gnpy.core.info import create_input_spectral_information, Pref, create_arbitrary_spectral_information, \
    ReferenceCarrier
from gnpy.core.equipment import trx_mode_params
from gnpy.core.exceptions import ConfigurationError
from gnpy.tools.json_io import network_from_json, load_equipment, load_network, _spectrum_from_json, load_json, \
    Transceiver, requests_from_json,network_from_params
from gnpy.topology.request import PathRequest, compute_constrained_path, propagate, propagate_and_optimize_mode
from support.gobal_control import GlobalControl
from gnpy.core.parameters import Parameters
import numpy as np
import matplotlib.pyplot as plt

from networkx import (dijkstra_path, NetworkXNoPath,
                      all_simple_paths, shortest_simple_paths)

gc = GlobalControl()
data_path = gc.gnpy_path/ 'data'
EQPT_FILENAME = data_path / 'eqpt_config.json'
NETWORK_FILENAME = data_path / 'testTopology_expected.json'

def net_setup(equipment):
    """common setup for tests: builds network, equipment and oms only once"""
    link_params = Parameters
    link_params.spans_per_oms = 4
    link_params.span_num = 10
    network = load_network(NETWORK_FILENAME, equipment)
    spectrum = equipment['SI']['default']
    p_db = spectrum.power_dbm
    p_db = -8
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

def test_equalization(case, deltap, target, mode, slot_width, equalization):
    """check that power target on roadm is correct for these cases; check on booster
    - SI : target_pch_out_db / target_psd_out_mWperGHz
    - node : target_pch_out_db / target_psd_out_mWperGHz
    - per degree : target_pch_out_db / target_psd_out_mWperGHz
    for these cases with and without power from user
    """
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
    if case == 'SI':
        delattr(equipment['Roadm']['default'], 'target_pch_out_db')
        setattr(equipment['Roadm']['default'], equalization, target_psd)
        network = net_setup(equipment)
    elif case == 'nodes':
        json_data = load_json(NETWORK_FILENAME)
        for el in json_data['elements']:
            if el['uid'] in roadms:
                el['params'] = {equalization: target_psd}
        network = network_from_json(json_data, equipment)
        spectrum = equipment['SI']['default']
        p_db = spectrum.power_dbm
        p_total_db = p_db + lin2db(automatic_nch(spectrum.f_min, spectrum.f_max, spectrum.spacing))
        build_network(network, equipment, p_db, p_total_db)
        # check that nodes not in roadms have target_pch_out_db not None
        pw_roadms = [r for r in network.nodes() if r.uid not in roadms and isinstance(r, Roadm)]
        for roadm in pw_roadms:
            assert roadm.target_psd_out_mWperGHz is None
            assert roadm.target_pch_out_dbm == target
        for roadm in [r for r in network.nodes() if r.uid in roadms and isinstance(r, Roadm)]:
            assert roadm.target_pch_out_dbm is None
            assert getattr(roadm, equalization) == target_psd
    path = compute_constrained_path(network, req)
    si = create_input_spectral_information(
        f_min=req.f_min, f_max=req.f_max, roll_off=req.roll_off, baud_rate=req.baud_rate, power=req.power,
        spacing=req.spacing, tx_osnr=req.tx_osnr, ref_carrier=ref)
    net_res = solve_network(path,si)

    plt.figure()
    plt.plot(net_res.frequency,net_res.sig_psd)
    plt.legend()
    
    plt.figure()
    plt.plot(net_res.sig_power_dB,'-o')
    plt.show()

def test_equalization_from_params(input_power,target, equalization):
    """check that power target on roadm is correct for these cases; check on booster
    - SI : target_pch_out_db / target_psd_out_mWperGHz
    - node : target_pch_out_db / target_psd_out_mWperGHz
    - per degree : target_pch_out_db / target_psd_out_mWperGHz
    for these cases with and without power from user
    """
    equipment = load_equipment(EQPT_FILENAME)
    setattr(equipment['Roadm']['default'], 'target_pch_out_db', target)
    target_psd = power_dbm_to_psd_mw_ghz(target, 32e9)
    ref = ReferenceCarrier(baud_rate=32e9, slot_width=50e9)

    delattr(equipment['Roadm']['default'], 'target_pch_out_db')
    setattr(equipment['Roadm']['default'], equalization, target_psd)

    link_params = Parameters
    link_params.spans_per_oms = 4
    link_params.span_num = 10
    network = network_from_params(link_params,equipment)
    spectrum = equipment['SI']['default']
    p_db = spectrum.power_dbm
    p_db = -8
    p_total_db = p_db + lin2db(automatic_nch(spectrum.f_min, spectrum.f_max, spectrum.spacing))
    build_network(network, equipment, p_db, p_total_db)
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
    # test_equalization(case='SI', deltap=-8, target=0, mode="mode 1", slot_width=50e9, equalization='target_psd_out_mWperGHz')
    test_equalization_from_params(input_power=0,target=0,equalization='target_psd_out_mWperGHz')