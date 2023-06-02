# WAN Example Data

The example data is generated from the original WAN data (provided by OAG, [Official Airline Guide](https://www.oag.com/airline-schedules-data)) by
* 1. select 30% random nodes
* 2. shuffle the fluxes + add noise to them
* 3. select largest component + symmetrize flow

The code on how the sample data was generated is shown below, thereby is the Graph `G_wan` the original graph provided by OAG
```python
import numpy as np
import ImportRisk as ir

# note that G_wan is a networkx.Graph 
nodes = np.array(list(G_wan.nodes))
# create a subgraph using a randomly selected thrid fraction
N0 = len(nodes)
sub_nodes = np.unique(np.random.randint(0, N0, int(3*N0/10))).astype(int)
sub_nodes = nodes[sub_nodes]
G_example = G_wan.subgraph(sub_nodes).copy()

edges_data = sorted(list(G_example.edges(data=True)))

e_tuples = [(e[0], e[1]) for e in edges_data]
e_values = np.array([e[2]['flux'] * 1 for e in edges_data])
# shuffle and randomize the values
np.random.shuffle(e_values)
Ne = len(e_tuples)
e_values += e_values * 0.2 * np.random.standard_normal(Ne)
e_values[e_values<0] = 50
e_values = np.floor(e_values)
assert (e_values < 0).sum() == 0, 'no negative fluxes allowed'

# now set new edge attributes
new_edge_attris = {e: {'flux': f} for e, f in zip(e_tuples, e_values)}
G_exa =  G_example.copy()
nx.set_edge_attributes(G_exa, new_edge_attris)
G_exa = ir.get_largest_component(ir.get_symmetric_graph(G_exa))
```
