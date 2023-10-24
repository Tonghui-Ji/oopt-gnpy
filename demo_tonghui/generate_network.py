from pathlib import Path
from gnpy.core.exceptions import NetworkTopologyError
from gnpy.core.network import build_network,span_loss
from gnpy.tools.json_io import load_equipment, load_network

import os, sys

script_path = os.path.dirname(os.path.abspath(__file__ if '__file__' in globals() else sys.executable))
gnpy_path = os.path.dirname(script_path)

eqpt_filename = Path(gnpy_path + '/data/eqpt_config.json')
network_filename = Path(gnpy_path + '/data/bugfixiteratortopo.json')

equipment = load_equipment(eqpt_filename)
network = load_network(network_filename, equipment)

build_network(network, equipment, 0, 20)

pass