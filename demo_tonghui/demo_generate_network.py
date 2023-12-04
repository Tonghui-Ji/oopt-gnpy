from pathlib import Path
import pytest
from gnpy.core.exceptions import NetworkTopologyError
from gnpy.core.network import span_loss
from gnpy.tools.json_io import load_equipment, load_network
from support.gobal_control import GlobalControl
import networkx as nx
import matplotlib.pyplot as plt
gc = GlobalControl()
data_path = gc.gnpy_path/ 'data'

EQPT_FILENAME = data_path / 'eqpt_config.json'
NETWORK_FILENAME = data_path / 'bugfixiteratortopo.json'

equipment = load_equipment(EQPT_FILENAME)
network = load_network(NETWORK_FILENAME, equipment)

subax1 = plt.subplot(121)
nx.draw(network, with_labels=True)
plt.show()
pass