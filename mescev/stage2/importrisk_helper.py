from sklearn.metrics.pairwise import haversine_distances  # 4 geodist
from math import radians  # 4 geodist
import pandas as pd
import numpy as np
from pathlib import Path
import networkx as nx
import json


def load_wan_example():
    '''
    for publication here an example WAN is loaded
    it is modified from the original WAN by:
        * 1. select 30% random nodes
        * 2. shuffle the fluxes + add noise to them
        * 3. select largest component + symmetrize flow
    '''
    d_data = Path(__file__).parent.parent.parent / 'data/stage2/wan_example_dat'
    # load edgelist
    f_name = d_data / 'nx_edgesList.csv'
    G = nx.read_edgelist(f_name)
    # load and add node attributes
    f_name = d_data / 'nx_nodesAttr.csv'
    nodes_attributes = json.load(f_name.open('r'))
    nx.set_node_attributes(G, nodes_attributes)
    return G


def get_airports_of_country_and_outflux(outbreak_country, G,
                                        out_of_outbreak_set=True,
                                        level='country_id'):
    """
    returns dictionary with key=outbreak airport
                            and
                            value=weight of airport
                                 (=influx of airport)
    INPUT:
        outbreak_country
            iso_alpha2 specification of outbreak country
        G WAN-networkx
            created by e.g.: G = load_wan(year=2019, top_traffic=0.95)
        out_of_outbreak_set bool
            flux of international flights to countries not part of outbreak_cntr
        level str
            defines the outbreak aggregation level if:
                -"country_id" -> outbreak_country="DE" (for example)
                -"region_id" -> outbreak_country=0 (for example)
                -"region_name" -> outbreak_country='Eastern Africa' (for example)
                -"continent_id" -> outbreak_country=0 (for example)
                -"continent_name" -> outbreak_country='Africa' (for example)
    OUTPUT:
        source_weight dict
            keys = iata of airports in the country
            values = outflux
            {'TXL' : 123123123,
             'BER' : 123123123,
             ...}
    """
    outbreak_airports = [
        key for key, info in G.nodes.items() if info[level] == outbreak_country
    ]
    # only use those outbreak airports that are in the network
    nodes = np.array(G.nodes())
    nodes_to_int = {node: i for i, node in enumerate(nodes)}
    outbreak_nds = [nodes_to_int[ap] for ap in outbreak_airports]
    from importrisk import get_flux_matrix
    flux_matrix = get_flux_matrix(G)
    if out_of_outbreak_set:
        others = list(set(range(len(flux_matrix))) - set(outbreak_nds))
        outbreak_wghts = flux_matrix[:, outbreak_nds][others].sum(axis=0)
    else:
        outbreak_wghts = flux_matrix.sum(axis=0)[outbreak_nds]
    source_weights = {apt: outbreak_wghts[i] for i, apt in enumerate(outbreak_airports)}
    return source_weights


def get_outflux_of_cntr(outbreak_country, G, out_of_outbreak_set=False):
    '''
    OUTPUT:
        out float
            summed outflux of all airports in the country
    '''
    out = get_airports_of_country_and_outflux(outbreak_country,
                                              G,
                                              out_of_outbreak_set=out_of_outbreak_set)
    return np.sum(list(out.values()))


def get_outflux_of_all_cntrs(G, out_of_outbreak_set=False):
    all_cntrs = np.unique([dicc['country_id']
                           for arpt, dicc in
                           dict(G.nodes(data=True)).items()])
    dat_outflux_cntrs = [{'source_iso2': c,
                          'outflux': get_outflux_of_cntr(c, G,
                                                         out_of_outbreak_set=out_of_outbreak_set)
                          }
                         for c in all_cntrs]
    df = pd.DataFrame(dat_outflux_cntrs)
    return df


def wan_2_geo_distance(network):
    '''
    wrapper for wan_2_geo_dist
    '''
    lat_lon = np.empty((len(network.nodes), 2))
    apts = []
    for i, [apt, info] in enumerate(network.nodes.items()):
        apts.append(apt)
        lat_lon[i] = [info["lat"], info["lon"]]
    df = wan_2_geo_dist(apts, lat_lon)
    return df


def wan_2_geo_dist(apts, lat_lon):
    """
    INPUT:
        apts list
            list of airport names (or country names...)
        lat_lon shape=[N, 2]
            array of latitute and longitude
    """
    assert len(apts) == len(lat_lon)
    lat_lon_rad = np.array([[radians(_) for _ in __] for __ in lat_lon])
    result = (
        haversine_distances(lat_lon_rad) * 6371
    )  # multiply by Earth radius to get kilometers
    df = pd.DataFrame(result, columns=apts, index=apts)
    return df
