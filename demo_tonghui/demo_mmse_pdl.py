from gnpy.core.parameters import Parameters
from gnpy.core.mmse_pdl import pdl_pen, h_pdl, prod_mat
from gnpy.core.utils import lin2db,db2lin,generate_random_numbers
import numpy as np
import matplotlib.pyplot as plt
from support.gobal_control import GlobalControl
from gnpy.tools.json_io import load_equipment, load_network
from gnpy.core.network import span_loss,build_network
from gnpy.core.elements import Transceiver, Fiber, Edfa, Roadm
from gnpy.core.info import create_input_spectral_information, ReferenceCarrier
import networkx as nx
from networkx import dijkstra_path

gc = GlobalControl()
data_path = gc.gnpy_path/ 'data'
EQPT_FILENAME = data_path / 'eqpt_config.json'
NETWORK_FILENAME = data_path / 'LinkforSnrEstimationTest.json'
con_in = 1
con_out = 0
input_power = 0
dest = 'trx F'
baud_rate = 32e9
grid_hz = 50e9

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
si = create_input_spectral_information(f_min=191.3e12, f_max=191.3e12 + 79 * grid_hz, roll_off=0.15,
                                        baud_rate=baud_rate, power=p, spacing=grid_hz, tx_osnr=None,
                                        ref_carrier=ReferenceCarrier(baud_rate=baud_rate, slot_width=grid_hz))
source = next(transceivers[uid] for uid in transceivers if uid == 'trx A')
sink = next(transceivers[uid] for uid in transceivers if uid == dest)
path = dijkstra_path(network, source, sink)

params = Parameters()
params.snr_trx_dB = 17
params.sig_power_dB = 3
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
for el in path:
    if isinstance(el,Roadm):
        link_config.append("wss")
    elif isinstance(el,Edfa):
        link_config.append("oa")
        oa_params = Parameters()
        oa_params.nf_dB = 4.5
        oa_params.gain_dB = 15
        oa_params.pdl_dB = 0.1
        oa_params_list.append(oa_params)
    elif isinstance(el,Fiber):
        link_config.append("Fiber")

        



link_config = ["wss","wss","oa","voa","fiber","fiu","oa"]
link_params.link_config = link_config
link_params.num_oa = link_config.count("oa")
link_params.num_wss = link_config.count("wss")
link_params.num_fiber = link_config.count("fiber")
link_params.num_fiu = link_config.count("fiu")
link_params.num_voa = link_config.count("voa")
params.link_params = link_params

wss_params_list = []
oa_params_list = []
fiu_params_list = []
fiber_params_list = []
voa_params_list = []
for i in range(link_params.num_wss):
    wss_params = Parameters()
    wss_params.loss_dB = 8
    wss_params.type = 'lcos_model'
    wss_params.bw_otf = 10e9
    wss_params.bw = 100e9
    wss_params.pdl_dB = 0.55
    wss_params_list.append(wss_params)

for i in range(link_params.num_oa):
    oa_params = Parameters()
    oa_params.nf_dB = 4.5
    oa_params.gain_dB = 15
    oa_params.pdl_dB = 0.1
    oa_params_list.append(oa_params)

for i in range(link_params.num_fiber):
    fiber_params = Parameters()
    fiber_params.loss_dB = 16
    fiber_params_list.append(fiber_params)

for i in range(link_params.num_fiu):
    fiu_params = Parameters()
    fiu_params.loss_dB = 0.7
    fiu_params.pdl_dB = 0.1
    fiu_params_list.append(fiu_params)

for i in range(link_params.num_voa):
    voa_params = Parameters()
    voa_params.loss_dB = 1
    voa_params.pdl_dB = 0.1
    voa_params_list.append(voa_params)

params.oa_params_list = oa_params_list[::-1]
params.wss_params_list = wss_params_list[::-1]
params.fiu_params_list = fiu_params_list[::-1]
params.fiber_params_list = fiber_params_list[::-1]
params.voa_params_list = voa_params_list[::-1]

snr_esti_dB_x,snr_esti_dB_y = pdl_pen(params)

pass