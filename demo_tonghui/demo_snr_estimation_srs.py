from pathlib import Path
from gnpy.core.parameters import SimParams
from gnpy.core.utils import db2lin,lin2db
from gnpy.core.elements import Transceiver, Fiber, Edfa, Roadm, RamanFiber
from gnpy.core.exceptions import NetworkTopologyError
from gnpy.core.network import span_loss,build_network
from gnpy.tools.json_io import load_equipment, load_network
from support.gobal_control import GlobalControl
from networkx import dijkstra_path
from numpy import mean, sqrt, ones
from gnpy.core.info import create_input_spectral_information, ReferenceCarrier
import networkx as nx
import matplotlib.pyplot as plt
from gnpy.tools.json_io import load_json

gc = GlobalControl()
data_path = gc.gnpy_path/ 'data'
EQPT_FILENAME = data_path / 'eqpt_config.json'
NETWORK_FILENAME = data_path / 'LinkforSnrEstimationTest.json'
SimParams.set_params(load_json(data_path / 'sim_params.json'))

from pandas import read_csv
from gnpy.core.parameters import SimParams
from numpy.testing import assert_allclose
from gnpy.tools.json_io import load_json
TEST_DIR = gc.gnpy_path

def test_raman_fiber():
    """ Test the accuracy of propagating the RamanFiber."""
    # spectral information generation
    spectral_info_input = create_input_spectral_information(f_min=191.3e12, f_max=196.1e12, roll_off=0.15,
                                                            baud_rate=32e9, power=2e-3, spacing=50e9, tx_osnr=40.0,
                                                            ref_carrier=ReferenceCarrier(baud_rate=32e9, slot_width=50e9))
    SimParams.set_params(load_json(TEST_DIR / 'data' / 'sim_params.json'))
    fiber = RamanFiber(**load_json(TEST_DIR / 'data' / 'test_science_utils_fiber_config.json'))

    # propagation
    spectral_info_out = fiber(spectral_info_input)

    p_signal = spectral_info_out.signal
    p_ase = spectral_info_out.ase
    p_nli = spectral_info_out.nli

    expected_results = read_csv(TEST_DIR / 'data' / 'test_raman_fiber_expected_results.csv')
    assert_allclose(p_signal, expected_results['signal'], rtol=1e-3)
    assert_allclose(p_ase, expected_results['ase'], rtol=1e-3)
    assert_allclose(p_nli, expected_results['nli'], rtol=1e-3)

def propagation(input_power, con_in, con_out, dest):
    equipment = load_equipment(EQPT_FILENAME)
    network = load_network(NETWORK_FILENAME, equipment)
    nx.draw(network, with_labels=True)

    build_network(network, equipment, 0, 20)

    # parametrize the network elements with the con losses and adapt gain
    # (assumes all spans are identical)
    for e in network.nodes():
        if isinstance(e, Fiber):
            loss = e.params.loss_coef * e.params.length
            e.params.con_in = con_in
            e.params.con_out = con_out
        if isinstance(e, Edfa):
            e.operational.gain_target = loss + con_in + con_out

    transceivers = {n.uid: n for n in network.nodes() if isinstance(n, Transceiver)}

    p = input_power
    p = db2lin(p) * 1e-3
    spacing = 50e9  # THz
    si = create_input_spectral_information(f_min=191.3e12, f_max=191.3e12 + 79 * spacing, roll_off=0.15,
                                           baud_rate=32e9, power=p, spacing=spacing, tx_osnr=None,
                                           ref_carrier=ReferenceCarrier(baud_rate=32e9, slot_width=50e9))
    source = next(transceivers[uid] for uid in transceivers if uid == 'trx A')
    sink = next(transceivers[uid] for uid in transceivers if uid == dest)
    path = dijkstra_path(network, source, sink)
    plt.figure()
    for el in path:
        si = el(si)
        plt.plot(si.frequency,lin2db(1e3*si.signal))
        print(el)  # remove this line when sweeping across several powers
    plt.show()
    edfa_sample = next(el for el in path if isinstance(el, Edfa))
    nf = mean(edfa_sample.nf)

    print(f'pw: {input_power} conn in: {con_in} con out: {con_out}',
          f'OSNR@0.1nm: {round(mean(sink.osnr_ase_01nm),2)}',
          f'SNR@bandwitdth: {round(mean(sink.snr),2)}')
    return sink, nf, path


if __name__ == '__main__':
    sink, nf, _ = propagation(input_power=0, con_in=1, con_out=0, dest='trx F')
    osnr = round(mean(sink.osnr_ase), 3)
    snr_nli = lin2db(1/(1.0 / db2lin(round(mean(sink.snr), 3)) - 1.0 / db2lin(osnr)))

    print('snr nli is:'+str(snr_nli)+' dB')
    # test_raman_fiber()

    pass