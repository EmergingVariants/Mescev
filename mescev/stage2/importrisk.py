# -*- coding: utf-8 -*-
"""
Efficient computation of import probabilities given a flux network.
"""

import networkx as nx
import numpy as np
from functools import partial
from pathlib import Path
import pandas as pd
import multiprocessing as mp
from importrisk_helper import wan_2_geo_distance, get_airports_of_country_and_outflux, get_outflux_of_all_cntrs

def get_weight_name(network):
    example_edge = list(network.edges(data=True))[0]  # = [start, end, {'weight': 11111}]
    example_edge_attributes = example_edge[-1]  # = {'weight': 11111}
    attr_name = list(example_edge_attributes.keys())[0]  # = 'weight'
    return attr_name


def get_flux_matrix(network):
    """
    Obtain the flux matrix from a networkx.DiGraph or networkx.Graph object.

    Parameters
    ----------
    network : nx.Graph or nx.DiGraph
        The flux network

    Returns
    -------
    flux_matrix : numpy.ndarray
        The weighted adjacency matrix where the entry ``flux_matrix[j,i]`` gives
        the flux from node `i` to node `j`. Will have shape ``(N, N)`` with `N`
        being the number of nodes.
    """
    attr_name = get_weight_name(network)
    flux_matrix = nx.to_numpy_array(network, weight=attr_name)
    return flux_matrix.T


def get_transition_matrix(flux_matrix):
    """
    Computes the random walk transition matrix for a given flux matrix.

    Parameters
    ----------
    flux_matrix : numpy.ndarray
        The weighted adjacency matrix where the entry ``flux_matrix[j,i]`` gives
        the flux from node `i` to node `j`. Must have shape ``(N, N)`` with `N`
        being the number of nodes.

    Returns
    -------
    transition_matrix : numpy.ndarray
        The transition matrix where the entry ``transition_matrix[j,i]`` gives
        the probability to jump from node `i` to node `j`. Will
        have shape ``(N, N)`` with `N` being the number of nodes.
    """
    k = get_node_outflux(flux_matrix)
    k[k == 0.0] = 1.0  # this is done in case there are columns with no nonzero entry
    transition_matrix = flux_matrix / k[None, :]
    return transition_matrix


def get_node_influx(flux_matrix):
    """
    Computes the total node influx for a given flux matrix.

    Parameters
    ----------
    flux_matrix : numpy.ndarray
        The weighted adjacency matrix where the entry ``flux_matrix[j,i]`` gives
        the flux from node `i` to node `j`. Must have shape ``(N, N)`` with `N`
        being the number of nodes.

    Returns
    -------
    node_influx : numpy.ndarray
        The transition matrix where the entry ``node_influx[i]`` gives
        the total influx to node `i`. Will have shape ``(N,)`` with `N`
        being the number of nodes.
    """
    return flux_matrix.sum(axis=1)


def get_node_outflux(flux_matrix):
    r"""
    Computes the total node outflux for a given flux matrix.

    Parameters
    ----------
    flux_matrix : numpy.ndarray
        The weighted adjacency matrix where the entry ``flux_matrix[j,i]`` gives
        the flux from node `i` to node `j`. Must have shape ``(N, N)`` with `N`
        being the number of nodes.

    Returns
    -------
    node_outflux : numpy.ndarray
        The transition matrix where the entry ``node_outflux[i]`` gives
        the total outflux of node `i`. Will have shape ``(N,)`` with `N`
        being the number of nodes.
    """
    return flux_matrix.sum(axis=0)


def get_effective_distance_graph(transition_matrix, effective_distance_offset=1):
    r"""
    Evaluates the effective distance graph tree given a random walk
    transition matrix.

    Parameters
    ----------
    transition_matrix : numpy.ndarray
        The transition matrix where the entry ``transition_matrix[j,i]`` gives
        the probability to jump from node `i` to node `j`. Must
        have shape ``(N, N)`` with `N` being the number of nodes.
    outbreak_location : :obj:`int`
        The node at which the pandemic broke out.
    effective_distance_offset : :obj:`float`, default : 1.
        The offset parameter :math:`d_0` in the equation :math:`d_{nm} = d_0 - \mathrm{ln}P_{nm}`.

    Returns
    -------
    effective_distance_graph : networkx.DiGraph
        Each edge ``effective_distance_graph.edge(source,target)`` contains the effective
        distance from node ``source`` to node ``target``.
    """

    N = transition_matrix.shape[0]
    D = nx.empty_graph(N, create_using=nx.DiGraph)
    for target, source in zip(*transition_matrix.nonzero()):
        d = effective_distance_offset - np.log(transition_matrix[target, source])
        D.add_edge(source, target, weight=d)

    return D


def get_shortest_path_tree(
    transition_matrix, outbreak_location, effective_distance_offset=1.0, geo_distances=None
):
    r"""
    Evaluates the shortest path tree on the effective distance matrix
    of the given transition matrix.

    Parameters
    ----------
    transition_matrix : numpy.ndarray
        The transition matrix where the entry ``transition_matrix[j,i]`` gives
        the probability to jump from node `i` to node `j`. Must
        have shape ``(N, N)`` with `N` being the number of nodes.
    outbreak_location : :obj:`int`
        The node at which the pandemic broke out.
    effective_distance_offset : :obj:`float`, default : 1.
        The offset parameter :math:`d_0` in the equation :math:`d_{nm} = d_0 - \mathrm{ln}P_{nm}`.
    geo_distances : numpy.ndarray
        geo_distances[0, 1] = goegraphic 'crow flight' distance between airport 0 and 1
        note that matrix is symmetric

    Returns
    -------
    shortest_path_tree : networkx.DiGraph
        The shortest path tree from the outbreak location to all other nodes.
    """

    D = get_effective_distance_graph(transition_matrix, effective_distance_offset)
    paths = nx.single_source_dijkstra_path(D, outbreak_location, weight="weight")

    N = D.number_of_nodes()
    T = nx.empty_graph(N, create_using=nx.DiGraph)

    T.nodes[outbreak_location]["distance"] = 0.0
    T.nodes[outbreak_location]["step_distance"] = 0
    compute_geo_dist = geo_distances is not None
    if compute_geo_dist:
        T.nodes[outbreak_location]["geo_path_distance"] = 0.0
    for p in paths.values():
        nx.add_path(T, p)
        for j in range(1, len(p)):
            T.nodes[p[j]]["distance"] = (
                T.nodes[p[j - 1]]["distance"] + D.edges[p[j - 1], p[j]]["weight"]
            )
            T.nodes[p[j]]["step_distance"] = T.nodes[p[j - 1]]["step_distance"] + 1
            if compute_geo_dist:
                T.nodes[p[j]]["geo_path_distance"] =  (T.nodes[p[j - 1]]["geo_path_distance"] +
                                                       geo_distances[p[j - 1], p[j]])

    return T


def estimate_total_descendant_population(
    shortest_path_tree, current_location, flux_vector, flux_scaling_exponent=1.0
):
    r"""
    Recursively estimate the cumulated population size of all descendants at the demanded location and write
    the result to the tree.

    Parameters
    ----------
    shortest_path_tree : networkx.DiGraph
        The shortest path tree from the outbreak location to all other nodes.
    current_location : :obj:`int`
        The location for which to compute the cumulated descendant population size.
    flux_vector : numpy.ndarray
        The flux vector from which the population sizes will be estimated.
    flux_scaling_exponent : :obj:`float`, default : 1.
        The scaling exponent :math:`\xi` to relate node flux to a population size
        (where :math:`N_n \propto F_n^\xi`).

    Notes
    -----

    Call this function once with the current location set to the outbreak location,
    i.e.

        >>> estimate_total_descendant_population(T, outbreak_location, flux_vector)

    The rest will be evaluated recursively.
    """

    T = shortest_path_tree
    T.nodes[current_location]["descendant_population"] = (
        flux_vector[current_location] ** flux_scaling_exponent
    )

    for child in T.successors(current_location):
        estimate_total_descendant_population(
            T, child, flux_vector, flux_scaling_exponent
        )
        T.nodes[current_location]["descendant_population"] += T.nodes[child][
            "descendant_population"
        ]


def get_exit_probabilities(
    shortest_path_tree, outbreak_location, flux_vector,
    flux_scaling_exponent=1.0
):
    r"""
    Get the exit probabilities at each airport, given the effective-distance shortest path tree and
    the outbreak location. Will be estimated from the node flux.

    Parameters
    ----------
    shortest_path_tree : networkx.DiGraph
        The shortest path tree from the outbreak location to all other nodes.
    outbreak_location : :obj:`int`
        The node at which the pandemic broke out.
    flux_vector : numpy.ndarray
        The flux vector from which the population sizes will be estimated.
    flux_scaling_exponent : :obj:`float`, default : 1.
        The scaling exponent :math:`\xi` to relate node flux to a population size
        (where :math:`N_n \propto F_n^\xi`).

    Returns
    -------
    exit_probabilities : numpy.ndarray
        The entry ``exit_probabilities[i]`` will contain the probability that
        a passenger exits at node `i` given it entered a plane at the given
        outbreak location.
    """

    T = shortest_path_tree
    estimate_total_descendant_population(
        T, outbreak_location, flux_vector, flux_scaling_exponent
    )

    N = T.number_of_nodes()
    q = np.zeros((N,))
    for node in range(N):
        # the root node will have zero exit probability
        if node == outbreak_location:
            continue
        elif (
            "descendant_population" not in T.nodes[node].keys()
        ):  # happens if node has no influx
            continue
        current_pop = flux_vector[node] ** flux_scaling_exponent
        q[node] = (
            current_pop
            / T.nodes[node]["descendant_population"]
        )

    return q


def estimate_fractions_of_transit_exit_spt(shortest_path_tree, transition_matrix_spt,
                                      exit_probs, current_location):
    r"""
    Recursively computes the fraction of transiting and exiting passengers from
    the outbreak region
    - assigns to each node of the shortest_path_tree:
        -"exit_fraction": fraction of passengers that started in
                            outbreak-location and leave at the node
        -"transit_fraction": fraction of passengers that started in
                            outbreak-location and use the node as transit
                            to other downstream  nodes

    Parameters
    ----------
    shortest_path_tree : networkx.DiGraph
        The shortest path tree from the outbreak location to all other nodes.
    current_location : :obj:`int`
        The location for which to compute the cumulated descendant population size.
    transition_matrix_spt : numpy.ndarray
        The transition matrix of the shortest_path_tree
            needed because the edges of the SPT have no attributes
    exit_probs : numpy.ndarray
        array of exit probabilities that were generated via "get_exit_probabilities"

    Notes
    -----

    Call this function once with the current location set to the outbreak location,
    i.e.

        >>> estimate_fractions_of_transit_exit_spt(shortest_path_tree, transition_matrix_spt,
        >>>                                        exit_probs, outbreak_location)

    The rest will be evaluated recursively.
    """
    T = shortest_path_tree
    parents = list(T.predecessors(current_location))
    if len(parents) == 0:
        T.nodes[current_location]['transit_fraction'] = 1
        T.nodes[current_location]['exit_fraction'] = 0
    else:
        arriver = (T.nodes[parents[0]]['transit_fraction'] *
                   transition_matrix_spt[parents[0], current_location])
        T.nodes[current_location]['exit_fraction'] =  arriver * exit_probs[current_location]
        T.nodes[current_location]['transit_fraction'] = arriver * (1 - exit_probs[current_location])
    for child in T.successors(current_location):
        estimate_fractions_of_transit_exit_spt(T, transition_matrix_spt, exit_probs, child)


def get_exit_and_transit_fraction_spt(
    shortest_path_tree, flux_matrix_of_original_graph, exit_probabilities, outbreak_location
):
    r"""
    Get the exit fraction of the shortest path tree (SPT) at each airport.
    The exit fraction is the fraction of people that started in the outbreak-location leaving at node n
    if they only move on the SPT.
    It is therefore analog to the import risk with the only difference that the exit fracion is only
    considering shortest paths and not all possible paths as the import risk.

    Parameters
    ----------
    shortest_path_tree : networkx.DiGraph
        The shortest path tree from the outbreak location to all other nodes.
    flux_matrix_of_original_graph : np.ndarray
        The flux matrix of the graph from which the distance graph and the shortest path is created
        i.e. G ---effectiveDistance---> D ---shortestPathFromOutbreak---> T
            flux_matrix_of_original_graph = flux matrix of G
    exit_probabilities : numpy.ndarray
        The exit probability computed via "get_exit_probabilities"
    outbreak_location : int
        node where the outbreak started (origin of SPT)

    Returns
    -------
    df : pandas.DataFrame
        index = node-id
        columns = ['exit_fraction', 'transit_fraction']
    """

    T = shortest_path_tree
    flux_matrix_spt = flux_matrix_of_original_graph.copy()
    flux_matrix_spt[nx.adjacency_matrix(T).toarray() != 1] = 0

    # create transition matrix
    outflux_spt = flux_matrix_spt.sum(axis=1)
    outflux_spt[outflux_spt == 0] = 1  # to prevent division by zero
    transition_matrix_spt = flux_matrix_spt / outflux_spt[:, None]

    # recursively compute
    estimate_fractions_of_transit_exit_spt(T, transition_matrix_spt,
                                           exit_probabilities,
                                           outbreak_location)
    df = pd.DataFrame(dict(T.nodes(data=True))
                     ).T[['exit_fraction', 'transit_fraction']]
    return df


def get_arrival_matrix(transition_matrix, exit_probabilities):
    """
    Given a transition matrix and the conditional
    exit probabilities, compute the arrival
    matrix.

    Parameters
    ----------
    transition_matrix : numpy.ndarray
        The transition matrix where the entry ``transition_matrix[j,i]`` gives
        the probability to jump from node `i` to node `j`. Must
        have shape ``(N, N)`` with `N` being the number of nodes.
    exit_probabilities : numpy.ndarray
        The entry ``exit_probabilities[i]`` will contain the probability that
        a passenger exits at node `i` given they started their walk at the defined
        outbreak location.

    Returns
    -------
    arrival_matrix : numpy.ndarray
        The arrival matrix where the entry ``arrival_matrix[j,i]`` gives
        the conditional probability that a walker jumped from node `i` to node `j`
        given that it did not exit at node `j` before.
        Will have shape ``(N, N)`` with `N` being the number of nodes.
    """
    return transition_matrix * (1 - exit_probabilities[None, :])


def get_limit_arrival_matrix(transition_matrix, exit_probabilities):
    """
    Given a transition matrix and the conditional
    exit probabilities, compute the arrival
    matrix after an infinite amount of steps.

    Parameters
    ----------
    transition_matrix : numpy.ndarray
        The transition matrix where the entry ``transition_matrix[j,i]`` gives
        the probability to jump from node `i` to node `j`. Must
        have shape ``(N, N)`` with `N` being the number of nodes.
    exit_probabilities : np.ndarray
        The entry ``exit_probabilities[i]`` will contain the probability that
        a passenger exits at node `i` given they started their walk at the defined
        outbreak location.

    Returns
    -------
    limit_arrival_matrix : numpy.ndarray
        The limit arrival matrix where the entry ``arrival_matrix[j,i]`` gives
        the conditional probability that a walker jumped from node `i` to node `j`
        given that it did not exit at node `j` before, for any time step `t`.
    """

    S = get_arrival_matrix(transition_matrix, exit_probabilities)
    return get_geometric_series(S)


def get_geometric_series(S):
    """
    Computes the geometric series from a matrix numerically
    NOTE: EXCLUDING the 0th term -> for real geometric series add `np.eye(N)`
    Parameters
    ----------
    S : numpy.ndarray
        matrix must have shape ``(N, N)`` with `N` being the number of nodes.

    Returns
    -------
    cut_geometric_series : numpy.ndarray
        'math':`\sum_{n=1}^{\infty} matrix^n`
    """
    N = S.shape[0]
    return np.linalg.inv(np.eye(N) - S) - np.eye(N)


def get_import_probabilities(
    limit_arrival_matrix, exit_probabilities, outbreak_location
):
    """
    Given the infinite-step arrival matrix and the conditional exit probabilities,
    compute the import probabilities of all target nodes.

    Parameters
    ----------
    limit_arrival_matrix : numpy.ndarray
        The limit arrival matrix where the entry ``arrival_matrix[j,i]`` gives
        the conditional probability that a walker jumped from node `i` to node `j`
        given that it did not exit at node `j` before, for any time step `t`.
        Must have shape ``(N, N)`` with `N` being the number of nodes.
    exit_probabilities : np.ndarray
        The entry ``exit_probabilities[i]`` will contain the probability that
        a passenger exits at node `i` given they started their walk at the defined
        outbreak location.

    Returns
    -------
    import_probabilities : numpy.ndarray
        The entry ``import_probabilities[i]`` will contain the ratio of people that
        will arrive at node `i` as their final destination given they entered a plane
        at the defined outbreak location.
    """

    return exit_probabilities * limit_arrival_matrix[:, outbreak_location]


def aggregate_passenger_steps(df, weight_arpts, name='mean_passenger_steps'):
    """
    computes the mean steps for a passenger to reach its final destination
    in the respective region (country)
    Formula:
        mean_passenger_steps_to_country = \sum_{arpts} (steps(arpt) * prob_to_reach_arpt_given_the_cntr * prob_to_start_in_source_arpt)
        thereby is prob_to_start_in_source_arpt=weight_arpts
        and the prob_to_reach_arpt_given_the_cntr needs to be computed
    INPUT:
        df pd.DataFrame
            columns = ['region', 'source', 'steps_from_source', 'import_prob', ...]
                with:
                    region = target region (e.g. target country)
                    source = source arpt
        weight_arpts np.ndarray
            the source weight, i.e. :math:`outflux_i / (\sum_j outflux_j)`
            with the sum over all sources
    OUTPUT:
        df pd.DataFrame
            columns = import_prob
            index = region
    """
    c_needed = ['region', 'source', 'steps_from_source', 'import_prob']
    assert np.all(np.in1d(c_needed, df.columns)), 'at least 1 essential column missing'
    df = df[c_needed].copy()
    df = df.merge(df.groupby((c_grp := ['region', 'source'])
                            ).sum()['import_prob'].reset_index(),
                  how='left', on=c_grp, suffixes=('', '_grp'))\
           .eval('prob_to_reach_arpt_given_the_cntr = import_prob / import_prob_grp')\
           .eval('steps_weighted = steps_from_source * prob_to_reach_arpt_given_the_cntr')
    df = aggregate_weighted_sum(df.region.values,
                                df.steps_weighted.values,
                                weight_arpts, name=name)
    return df


def aggregate_weighted_sum(region, import_prob, source_weight, name='import_prob'):
    """
    INPUT:
        source_weight np.ndarray
            the source weight, i.e. :math:`outflux_i / (\sum_j outflux_j)`
            with the sum over all sources
        region list
            region the respective airport is in
            e.g. country in which the respective airport is in
        import_prob np.ndarray
            import probability for each targest from each source
    OUTPUT:
        df pd.DataFrame
            columns = import_prob
            index = region
    """
    df = (
        pd.DataFrame(
            data=np.vstack([region, source_weight * import_prob]).T,
            columns=["region", name],
        )
        .groupby("region")
        .sum()
    )
    return df


def aggregate_min_distance(region, distance, name="dist_min"):
    """
    takes the minimal distance between all source
    nodes and a specific region
    INPUT:
        region list
            region the respective airport is in
            e.g. country in which the respective airport is in
        distance np.ndarray
            import probability for each targest from each source
    OUTPUT:
        df pd.DataFrame
            columns = import_prob
            index = region
    """
    df = (
        pd.DataFrame(data=np.vstack([region, distance]).T, columns=["region", name])
        .groupby("region")
        .min()
    )
    return df


def aggregate_multipath_distance(region, distance, name="dist_mp"):
    """
    takes the multi path distance between all source
    nodes and a specific region
    The arrival time is shortened if the walker
    can reach the target by 2 paths (a and b):
        :math:`e^{-D} = e^{-D_a} + e^{-D_b}`
        following [Gatreau et. al, Journal of Theoretical Biology (2008)]
    -> multiple path reduce the effective distance as follows
        :math:`D = -ln( \sum_a e^{-D_a} )`
    INPUT:
        region list
            region the respective airport is in
            e.g. country in which the respective airport is in
        distance np.ndarray
            import probability for each targest from each source
    OUTPUT:
        df pd.DataFrame
            columns = import_prob
            index = region
    """
    dist_exp = np.exp(-distance)
    df = (
        pd.DataFrame(data=np.vstack([region, dist_exp]).T, columns=["region", "foo"])
        .groupby("region")
        .sum()
    )
    df[name] = -np.log(df.foo)
    return df.drop(columns=["foo"])


class ImportRiskAnalyzer:

    outbreak_location = None
    effective_distance_offset = 1.0
    flux_scaling_exponent = 1.0
    estimate_population_with_outflux = True
    exit_prob_weighted_by = False

    flux_matrix = None
    transition_matrix = None
    arrival_matrix = None
    limit_arrival_matrix = None
    import_probabilities = None
    exit_probabilities = None
    shortest_path_tree = None

    def __init__(
        self,
        network,
        outbreak_location,
        effective_distance_offset=1.0,
        flux_scaling_exponent=1.0,
        estimate_population_with_outflux=True,
        global_exit_probabilities=None,
        exit_prob_weighted_by=False,
        dist_eff=True,
        dist_import=True,
        dist_exit_fractions=True,
    ):
        np.seterr(divide='ignore')  # to ignore warnings with zero division (happens for source country)
        self.network = network
        self.node_attributes = list( list( self.network.nodes(data=True))[0][1].keys())
        self.flux_matrix = get_flux_matrix(self.network)
        self.influx = get_node_influx(self.flux_matrix)
        self.outflux = get_node_outflux(self.flux_matrix)
        nodes = np.array(self.network.nodes())
        self.nodes_to_int = {node: i for i, node in enumerate(nodes)}
        self.int_to_nodes = {i: node for i, node in enumerate(nodes)}
        self.transition_matrix = get_transition_matrix(self.flux_matrix)
        self.geo_dist = self.compute_dist_geo()

        self.has_global_exit_probabilites = global_exit_probabilities is not None
        if self.has_global_exit_probabilites:
            self.exit_probabilities = global_exit_probabilities

        self.outbreak_location = self.nodes_to_int[outbreak_location]
        self.effective_distance_offset = effective_distance_offset
        self.flux_scaling_exponent = flux_scaling_exponent
        self.estimate_population_with_outflux = estimate_population_with_outflux
        self.exit_prob_weighted_by = exit_prob_weighted_by
        assert ((c := self.exit_prob_weighted_by) in
                (cc := [None, False, 'eff_dist', 'geo_dist'])), f'{c} not in {cc}'
        # distance flags
        self.dist_eff = dist_eff
        self.dist_import = dist_import
        self.dist_exit_fractions = dist_exit_fractions

    def compute(
        self,
        node_population=None,
    ):

        if not self.has_global_exit_probabilites:
            self.shortest_path_tree = self.compute_shortest_path_tree()

        self.exit_probabilities = self.compute_exit_probabilities(population=node_population)

        self.arrival_matrix = get_arrival_matrix(
            self.transition_matrix, self.exit_probabilities
        )
        self.limit_arrival_matrix = get_limit_arrival_matrix(
            self.transition_matrix, self.exit_probabilities
        )
        self.import_probabilities = get_import_probabilities(
            self.limit_arrival_matrix, self.exit_probabilities, self.outbreak_location
        )

        return self.import_probabilities


    def compute_exit_probabilities(self, population=None):
        # if glaobl_exit_probabilities are set -> no computing needed
        if self.has_global_exit_probabilites:
            self.exit_probabilities = np.array(self.global_exit_probabilities)
            self.exit_probabilities[self.outbreak_location] = 0
            return  self.exit_probabilities

        # given flags what is the scaling exponent and population
        scaling = self.flux_scaling_exponent
        if population is not None:
            scaling = 1.0
        elif self.estimate_population_with_outflux:
            population = self.outflux
        else:
            population = self.influx

        # add weight to population pop = N_i**scaling / d_{i,n0}
        if (by := self.exit_prob_weighted_by) not in [False, None]:
            if by == 'geo_dist':
                if self.geo_dist is None:
                    raise ValueError("if exit_prob_weighted_by_geo_dist: geo_dist must be computable")
                weights = 1 / self.geo_dist[self.outbreak_location]  # 1/d_{i,n0}
            elif by == 'eff_dist':
                T = self.shortest_path_tree
                df_ = pd.DataFrame(dict(T.nodes(data=True))).T
                weights = 1 / df_.distance.values     # 1 / d_{i,n0}
            weights[self.outbreak_location] = 1
            weights /= np.mean(weights)     # keeps numbers normal
            population = population**scaling * weights
            scaling = 1     # because scaling already applied to population

        # compute exit_probability
        self.exit_probabilities = get_exit_probabilities(
            self.shortest_path_tree, self.outbreak_location, population, scaling
        )
        return self.exit_probabilities

    def compute_dist_geo(self):
        self.geo_dist = None
        if np.all(np.in1d(['lat', 'lon'],self.node_attributes)):
            # return a pd.DataFrame with columns and index = iatas
            geo_dist = wan_2_geo_distance(self.network)
            # ensure correct order so np.ndarray can be used
            arpts_right_order = [self.int_to_nodes[i] for i in range(len(self.int_to_nodes))]
            self.geo_dist = geo_dist[arpts_right_order].T[arpts_right_order].T.values
        return self.geo_dist

    def compute_shortest_path_tree(self):
        self.shortest_path_tree = get_shortest_path_tree(
            self.transition_matrix,
            self.outbreak_location,
            self.effective_distance_offset,
            geo_distances=self.geo_dist
            )
        return self.shortest_path_tree

    def get_SPTree_measure(self, measure):
        assert measure in ['step_distance',
                           'geo_path_distance',
                           'distance'], f'{measure} not derived from SPT'
        if self.shortest_path_tree is None:
            self.shortest_path_tree = self.compute_shortest_path_tree()
        T = self.shortest_path_tree
        # columns: ['distance', 'step_distance', 'descendant_population']
        df = pd.DataFrame(dict(T.nodes(data=True))).T
        return df[measure].values

    def compute_effective_distance(self):
        return self.get_SPTree_measure('distance')

    def compute_step_distance(self):
        return self.get_SPTree_measure('step_distance')

    def compute_geo_path_distance(self):
        return self.get_SPTree_measure('geo_path_distance')

    def compute_exit_and_transit_fraction_spt(self):
        '''
        see doc from "get_exit_and_transit_fraction_spt"
        '''
        if not self.has_global_exit_probabilites:
            if self.shortest_path_tree is None:
                self.shortest_path_tree = self.compute_shortest_path_tree()
            if self.exit_probabilities is None:
                self.exit_probabilities = self.compute_exit_probabilities()
        df = get_exit_and_transit_fraction_spt(
            self.shortest_path_tree, self.flux_matrix.T,
            self.exit_probabilities, self.outbreak_location
            )
        return df


    def compute_for_multiple_sources(
        self,
        source_weights,
        node_population=None,
    ):
        """
        Returns exit probabilities for multiple sources.

        Parameters
        ----------
        source_weights : :obj:`dict`
            A dictionary mapping source nodes' IATA code to the
            weight with which it should be considered in the analysis.
            Weights could be number of people departing from the source
            node.

        Returns
        -------
        mean_import_probabilities : numpy.ndarray
            The entry ``import_probabilities[i]`` will contain the ratio of people that
            will arrive at node `i` as their final destination. Average weighted
            by the total number of people that depart from the given source nodes.

            Map back to IATA codes by using ``Risk.int_to_nodes[i]``.
        """

        N_source = len(source_weights)
        N = self.flux_matrix.shape[0]
        sources = ["" for _ in range(N_source)]
        weights = np.zeros((N_source,), dtype=float)

        for i, (k, v) in enumerate(source_weights.items()):
            sources[i] = k
            weights[i] = v

        weights /= weights.sum()
        p = np.zeros((N_source, N), dtype=float)

        for i, source in enumerate(sources):
            print(f"Computing exit for source {source} ({i+1}/{N_source}) ", end="\r")
            self.set_outbreak_location(source)
            _p = self.compute( node_population=node_population)
            p[i, :] = _p

        return weights.dot(p)

    def set_all_distance_flags(self, boo):
        self.dist_eff = boo
        self.dist_import = boo
        return self

    def set_global_exit_probabilities(self, global_exit_probabilities):
        if type(global_exit_probabilities) == dict:
            assert len(global_exit_probabilities) == len(self.influx)
            q = np.zeros_like(self.influx, dtype=float)
            for node, _q in global_exit_probabilities.items():
                q[self.nodes_to_int[node]] = _q
        else:
            q = global_exit_probabilities
        self.global_exit_probabilities = q
        self.has_global_exit_probabilites = True
        return self

    def set_outbreak_location(self, outbreak_location):
        '''
        if the outbreak location is changed
            -> the shortest path tree measures needs recomputation
        '''
        new = self.nodes_to_int[outbreak_location]
        if self.outbreak_location != new:
            self.outbreak_location = new
            self.shortest_path_tree = self.compute_shortest_path_tree()
        return self

    def set_effective_distance_offset(self, effective_distance_offset):
        if self.effective_distance_offset != effective_distance_offset:
            self.effective_distance_offset = effective_distance_offset
            self.shortest_path_tree = self.compute_shortest_path_tree()
        return self

    def set_flux_scaling_exponent(self, flux_scaling_exponent):
        if self.flux_scaling_exponent != flux_scaling_exponent:
            self.flux_scaling_exponent = flux_scaling_exponent
        return self

    def get_df(self, import_risk, name="import_prob"):
        # create dataframe with iata-target-id in target column
        result = [[self.int_to_nodes[i], ir] for i, ir in enumerate(import_risk)]
        df = pd.DataFrame(result, columns=["target", name])
        return df

    def get_df_distances(self):
        """
        join import risk and effective distance of a specific airport
        OUTPUT:
            df pd.DataFrame
                columns = ['target', 'import_prob', 'distance']
        """
        assert np.any(
            [self.dist_eff, self.dist_import]
        )
        out = {}
        if self.dist_import:
            out["import_prob"] = self.compute()
            out["dist_import_prob"] = -np.log(out["import_prob"])
        if self.dist_eff:   # all below measures are computed together (via SPT)
            out["distance"] = self.compute_effective_distance()
            out["steps_from_source"] = self.compute_step_distance()
            out["geo_path_dist"] = self.compute_geo_path_distance()
        if self.dist_exit_fractions:
            df_ = self.compute_exit_and_transit_fraction_spt()
            out[c] = df_[(c := 'transit_fraction')]
            out[c] = df_[(c := 'exit_fraction')]
            out['dist_exit_fraction'] = -np.log(out['exit_fraction'])
        for i, (k, v) in enumerate(out.items()):  # that could be simpler...
            df_ = self.get_df(v, name=k)
            df = df.merge(df_, how="left", on="target") if i else df_
        return df

    def collect_compute_multiple_sources(self, outbreak_airports):
        """
        for each outbreak airport saves the effective dist and import risk
        to be able to compare the import risk with effective distance
        and to compute country-level aggregation
        Note: the difference to "compute_for_multiple_sources" is
              that it does NOT aggregate
        INPUT:
            outbreak_airports dict
                keys=iata of airport
                values=weights of airports (as outflux)
        OUTPUT
            df_ pd.DataFrame
                columns = ['target',
                           'import_prob', 'distance',
                           'source', 'source_outflux']
        """
        for i, (otbrk_arpt, flux) in enumerate(outbreak_airports.items()):
            self.set_outbreak_location(otbrk_arpt)
            # cols=['target', 'import_prob', 'distance', ...]
            df = self.get_df_distances()
            df["source"] = otbrk_arpt
            df["source_outflux"] = flux
            df_ = pd.concat([df_, df]) if i else df
        return df_

    def aggregate_results_by_country(self, outbreak_airports, arpt_2_region):
        """
        IMPORT RISK aggregation per country is straight forward
            (use the outflux of the source airports as weight)
        EFFECTIVE DISTANCE is more complicated because it is
            proportional to the arrival time which scales differently with
            the number of passengers and the number of differnt paths
            available (see Iannelli et al., 2017 PRE)
        INPUT:
            outbreak_airports dict
                keys=iata of airport
                values=weights of airports (as outflux)
            arpt_2_region dict
                keys=iata of airport
                values=region (e.g. "country_id" or "region_id" or "continent")
        """

        def pd_2_np(series_list):
            return [s.values for s in series_list]

        if self.flux_scaling_exponent != 1:
            outbreak_airports = {
                arpt: weight**self.flux_scaling_exponent
                for (arpt, weight) in outbreak_airports.items()
            }
        df = self.collect_compute_multiple_sources(outbreak_airports)
        total_outflux = np.sum(list(outbreak_airports.values()))
        # join regions for aggregation
        region = pd.DataFrame.from_dict(
            arpt_2_region, orient="index", columns=["region"]
        )
        df = df.merge(region,
                      how="left", left_on="target", right_index=True)

        weight = df.source_outflux / total_outflux
        dfs_dists = []
        # IMPORT RISK aggregate
        if (c := "import_prob") in df.columns:
            args = pd_2_np([df.region, df.import_prob, weight])
            dfs_dists.append(aggregate_weighted_sum(*args))
            if 'steps_from_source' in df.columns:
                dfs_dists.append(aggregate_passenger_steps(df, weight))
        if (c := 'steps_from_source') in df.columns:
            args = pd_2_np([df.region, df[c]])
            dfs_dists.append(aggregate_min_distance(*args, name='min_steps'))
        if (c := "exit_fraction") in df.columns:
            args = pd_2_np([df.region, df[c], weight])
            dfs_dists.append(aggregate_weighted_sum(*args, name='import_prob_spt'))

        # DISTANCES aggregate
        c = 'distance' # = effective distance
        args = pd_2_np([df.region, df[c]])
        df_ = aggregate_min_distance(*args)
        dfs_dists.append(df_)
        df_ = aggregate_multipath_distance(*args)
        dfs_dists.append(df_)
        # aggregates + accounting for source_outflux
        args = pd_2_np([df.region, df[c] - np.log(weight)])
        df_ = aggregate_multipath_distance(*args, name='dist_Nmp')
        dfs_dists.append(df_)

        # JOIN MEASURES
        df_cntr = pd.concat(dfs_dists, axis=1)
        return df_cntr

    def compute_eff_distance_all_sources(self):
        '''
        returns the distance measures for all sources
        '''
        for i, node in self.int_to_nodes.items():
            self.outbreak_location = i
            dist = -np.log(self.compute())
            df_ = self.get_df(dist, name=node)
            df_dir = df_dir.merge(df_, how="left", on="target") if i else df_
            # now the same for eff-distance
            dist = self.compute_effective_distance()
            df_ = self.get_df(dist, name=node)
            df_def = df_def.merge(df_, how="left", on="target") if i else df_
        return df_dir, df_def


def normalize_probs_4_subset(df, cntrs):
    '''
    given the subset cntrs-> normalize
    INPUT:
        df pd.DataFrame
            columns = [source iso2-countries]
            index = [target iso2-countries]
            values = import probabilities
        cntrs list or None
            subset of countries
    OUTPUT:
        df pd.DataFrame
            as input but normalized subset
    '''
    if cntrs is None:
        cntrs = np.array(df.columns)
    cols = np.array(df.columns)
    cols = np.sort(cols[np.in1d(cols, cntrs)])
    df = df[cols].T[cols].T.copy().fillna(0)
    dat = df.values
    outflux = dat.sum(axis=0)
    df = pd.DataFrame(data=(dat / outflux),
                      columns=cols, index=cols)
    return df


def import_prob_2_flux(df_ip, df_outflux, return_numpy_input=False):
    '''
    INPUT:
        df_ip pd.DataFrame
            cols  =['AE', ..., 'DE', ..., 'US'] ('source_iso2')
            index =['AE', ..., 'DE', ..., 'US'] ('target_iso2')
            values='import_prob'
        df_outflux
            cols  =['AE', ..., 'DE', ..., 'US'] ('source_iso2')
            index =['outflux']
    OUTPUT:
        pd.DataFrame
            cols  =['AE', ..., 'DE', ..., 'US'] ('source_iso2')
            index =['AE', ..., 'DE', ..., 'US'] ('target_iso2')
    '''
    # ensure same columns
    cols = df_ip.columns
    df_outflux = df_outflux[cols]
    df_ip = df_ip[cols].T[cols].T
    if return_numpy_input:
        return df_ip.values, df_outflux.values.flatten(), cols
    flux = df_ip.values * df_outflux.values
    return pd.DataFrame(data=flux,
                        columns=cols,
                        index=cols.rename(df_ip.index.name))


def symmetrize_import_prob(df_ip, df_outflux, N=3):
    '''
    corrects import probabilities by
    1. convert probability to fluxes
    2. symmetrize fluxes
    3. convert back to probabilites
    4. repeat steps 1-3 N times
    INPUT:
        df_ip pd.DataFrame
            cols  =['AE', ..., 'DE', ..., 'US'] ('source_iso2')
            index =['AE', ..., 'DE', ..., 'US'] ('target_iso2')
            values='import_prob'
        df_outflux
            cols  =['AE', ..., 'DE', ..., 'US'] ('source_iso2')
            index =['outflux']
    '''
    prob, outflux, cols = import_prob_2_flux(df_ip, df_outflux,
                                             return_numpy_input=True)
    for i in range(N):
        flux = prob * outflux   # numpy: multiply with last dimension = source
        flux = flux + flux.T
        prob = flux / flux.sum(axis=0)
    df_ip = pd.DataFrame(data=prob,
                         columns=cols,
                         index=cols.rename(df_ip.index.name))
    return df_ip


def import_risk_pivot_and_normalize(df, cntrs, tolerance=1e-5,
                                    c_i='target_iso2',
                                    c_c='source_iso2',
                                    c_v='import_prob'):
    '''
    what it does:
        - 0. pivot the data such that col=source, indx=target
        - 1. ensures that NaN -> NA "Namibia"
        - 2. ensures that total prob. <= 1
        - 3. set self-loops = 0 AND normalize probs (\sum p = 1)
    INPUT:
        df pandas.DataFrame
            either
                cols = ['target_iso2', 'source_iso2', 'import_prob']
            or
                cols  =['AE', ..., 'DE', ..., 'US'] ('source_iso2')
                index =['AE', ..., 'DE', ..., 'US'] ('target_iso2')
                values='import_prob'
        cntrs list or None
            iso2 code for countries for which reference data exists
            if None: all countries are considered
    OUTPUT:
        df pandas.DataFrame
            cols  =['AE', ..., 'DE', ..., 'US'] ('source_iso2')
            index =['AE', ..., 'DE', ..., 'US'] ('target_iso2')
            values='import_prob'
            BUT only with countries that are also in "cntrs"
    '''
    # if the values are not only floats -> pivot needed
    if df.values.dtype != float:
        df = df.reset_index()   # only needed in case pivot needs to be done
        needed_cols = [c_i, c_c, c_v]
        assert np.all(np.in1d(needed_cols, df.columns)), f'missing columns (columns={df.columns}, needed_columns= {needed_cols})'
        df.fillna({c_i: 'NA', c_c: 'NA'}, inplace=True)
        df = df.pivot(columns=c_c, index=c_i, values=c_v)
    df.rename(columns={np.nan: 'NA'}, index={np.nan: 'NA'}, inplace=True)
    assert np.all(df.sum(axis=0) <= 1 + tolerance), 'the total probability is larger 1 for at least 1 source'
    assert np.all(df.columns == df.index), 'columns and index are not the same'
    # set import to source to 0 (country aggregated: unequal zero if multiple airports in country)
    np.fill_diagonal(df.values, 0)
    df = normalize_probs_4_subset(df, cntrs)
    return df


def import_risk_cntr_aggregate(G, arpt_2_region, cntr,
                               paras_risk=None):
    '''
    computes the import risk aggregated for 1 country as source
    '''
    if paras_risk is None:
        paras_risk = {'flux_scaling_exponent': 1}
    dummy_outbreak = list(G.nodes())[0]  # anyway changed in the function
    Risk = ImportRiskAnalyzer(G, dummy_outbreak, **paras_risk)
    source_weights = get_airports_of_country_and_outflux(cntr, G)
    df_out = Risk.aggregate_results_by_country(source_weights, arpt_2_region)
    df_out['source_iso2'] = cntr
    df_out['source_outflux'] = np.sum(list(source_weights.values()))
    df_out.index.name = 'target_iso2'
    return df_out


def import_risk_all_cntr_aggregate(G, symmetrize=True, paras_risk=None,
    parallel=True, worker_fraction=2/3, return_all=False):
    '''
    computes the import risk aggregated on country-level
    INPUT:
        G networkx.Graph
            WAN network with
                airports = nodes
                max passenger capaciity = links
        paras_risk dict or None
            dictionary with parameters for the ImportRiskAnalyzer object
            None: is set to default
        symmetrize bool
            flag to symmetrize on flux-level
            i.e.:
                flux_matrix = p_import_matrix * ourflux_vector[None, :]
                flux_matrix += flux_matrix.T
                p_import_matrix = flux_matrix / flux_matrix.sum(axis=0)[None, :]
        parallel bool
            if parallel computing is used or not
        worker_fraction float < 1
            faction of workers used in multiprocessing
            (already some multiprocessing at use from networkx
             -> using all worker is slower than a smaller fraction)
    OUTPUT:
        df_out pd.DataFrame
            columns= source_iso2
            index  = target_iso2
            values = import_risk
    '''
    # line below only works if nodes have country-information attributes
    df_arptInfo = pd.DataFrame([dic for iata, dic
                                in dict(G.nodes(data=True)).items()])
    arpt_2_region = df_arptInfo.set_index('iata')['country_id'].to_dict()
    # sort countries by airport-size
    #   (process largest first so 1 worker will not handle largest country at the end)
    cntrs = df_arptInfo.groupby('country_id').count()\
                       .sort_values('id', ascending=True).index
    # prepare computation
    aggregate_cntr_partial = partial(import_risk_cntr_aggregate,
                                     G,
                                     arpt_2_region,
                                     paras_risk=paras_risk)
    # COMPUTE
    if parallel:
        pool = mp.Pool(int(mp.cpu_count() * worker_fraction))
        out = pool.map(aggregate_cntr_partial, cntrs)
    else:
        out = []
        for cntr in cntrs:
            out.append(aggregate_cntr_partial(cntr))
    # join the dataframes
    df_out = pd.concat(out)
    # pivot the dataframe: cols=source index=target
    df_out = import_risk_pivot_and_normalize(df_out, None)
    # symmetrize
    if symmetrize:
        df_outflux = get_outflux_of_all_cntrs(G,
                                              out_of_outbreak_set=True)
        df_outflux = df_outflux.set_index('source_iso2').T
        df_out_sym = symmetrize_import_prob(df_out, df_outflux)
        if return_all:
            return df_out, df_out_sym
    return df_out


# local variables
d_data_ir = Path(__file__).absolute().parent.parent.parent / 'data/stage2/import_risks'


def round_date_2_month(date, fraction=2/3):
    lim = 30 * fraction
    d = pd.to_datetime(date)
    return f'{d.year}-{d.month + int(d.day > lim)}-01'


def get_file_names_precomputed(date, paras_risk=None,
                               top_traffic=0.99,
                               d_data=d_data_ir):
    '''
    INPUT:
        date
        paras_risk dict
            dictionary of optional parameters for ImportRiskAnalyzer
        d_data pathlib.Path object
            directory to where the precomputed data is located
    '''
    if paras_risk is None:
        paras_risk = {'exit_prob_weighted_by': 'geo_dist'}
    f_id = f'top_traffic_{top_traffic}__' + '__'.join([f'{k}_{v}' for k, v in paras_risk.items()])
    d = pd.to_datetime(round_date_2_month(date))
    dd = f'{d.year}-{d.month}'
    assert d_data.exists(), f'folder {d_data} does not exist'
    f_name_ir      = d_data / (f_id + f'__import_risk_{dd}.csv')
    f_name_outflux = d_data / (f_id + f'__outflux_{dd}.csv')
    return f_name_ir, f_name_outflux, d



def precomputed_import_risk_cntr_agg_all(date):
    '''
    load import risk that was computed previously via   "import_risk_all_cntr_aggregate"
    (takes between 10-20 minutes to compute 1 month)
    INPUT:
        date str
            e.g.: '2022-01-01' or '2022-01-22'
    OUTPUT:
        df_ir pd.DataFrame
            columns = source_iso2
            index   = target_iso2
        df_outflux pd.
    '''
    # standard parameters used in computation
    top_traffic = 0.99
    paras_risk = {'exit_prob_weighted_by': 'geo_dist'}
    f_name_ir, f_name_outflux, d = get_file_names_precomputed(date, paras_risk=paras_risk,
                                                          top_traffic=top_traffic,
                                                          d_data=d_data_ir)
    # load import risk
    f_name = f_name_ir
    assert f_name.exists(), f'{f_name} does not exists! Computation needed via "import_risk_all_cntr_aggregate"'
    df_ir = pd.read_csv(f_name)
    df_ir = df_ir.rename(columns={df_ir.columns[0]: (ci:= 'target_iso2')})\
                 .set_index(ci)
    # load outflux matrix
    f_name = f_name_outflux
    assert f_name.exists(), f'{f_name} does not exist! this should not happen since outflux is always computed with import risk'
    df_outflux = pd.read_csv(f_name)
    df_outflux = df_outflux.drop(columns=[df_outflux.columns[0]])\
                           .set_index('source_iso2')
    return df_ir, df_outflux


def precomputed_import_risk_cntr_agg(date, cntr):
    df_ir, df_outflux = precomputed_import_risk_cntr_agg_all(date)
    assert cntr in df_ir.columns, f'{cntr} is wrong iso2 code for country (or not in network)'
    return df_ir[cntr], df_outflux['outflux'][cntr]


if __name__ == "__main__":

    # create example test graph
    G = nx.empty_graph(4, nx.DiGraph)
    G.add_edge(0, 1, flux=0.1)
    G.add_edge(0, 2, flux=0.2)
    G.add_edge(1, 2, flux=0.3)
    G.add_edge(2, 3, flux=1.0)

    # define outbreak location, d0, and xi
    n0 = 0
    d0 = 0.5
    xi = 1.0

    F = get_flux_matrix(G)
    P = get_transition_matrix(F)
    f_in = get_node_influx(F)
    T = get_shortest_path_tree(P, n0, d0)
    q = get_exit_probabilities(T, n0, f_in, xi)
    S_inf = get_limit_arrival_matrix(P, q)
    p = get_import_probabilities(S_inf, q, n0)
    print("got:", p)
    print("expected:", [0, 1 / 3, 2 / 9, 4 / 9])

    with_outflux = False # example is designed for influx
    IRA = ImportRiskAnalyzer(G, n0, d0, xi,
                             estimate_population_with_outflux=with_outflux)
    p = IRA.compute()
    print("got:", p)
    print("expected:", [0, 1 / 3, 2 / 9, 4 / 9])

    p = IRA.compute_for_multiple_sources({0: 1})
    print("got:", p)
    print("expected:", [0, 1 / 3, 2 / 9, 4 / 9])

    # try out example WAN (30% nodes + shuffled and 20%-noise-added edges)
    # this code: 1. computes exit probability between all airports
    #            2. aggregates it on country level
    from importrisk_helper import load_wan_example
    G = load_wan_example()
    paras_risk = {'exit_prob_weighted_by': 'geo_dist'}
    df = import_risk_all_cntr_aggregate(G, paras_risk=paras_risk)
    # display results
    N = 10
    print('Given a subgraph of 30% randomly selected nodes with their edges shuffled')
    print(f'The {N} countries with the highest mean import probability as a target are:')
    print((df.mean(axis=1).sort_values() * 100).tail(N))

    # test loading of precompiled import_risk files
    date = '2021-01-01'
    cntr = 'DE'
    N = 10
    df_ir, df_outflux = precomputed_import_risk_cntr_agg(date, cntr)
    print(f'the {N} countries with the highest precomputed import risk from {cntr} at {date} are')
    print(df_ir.sort_values().tail(N))
