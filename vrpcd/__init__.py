# Copyright 2015 Thodoris Sotiropoulos
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
==================================================
Vehicle Routing Problem with Cross-Docking (VRPCD)
==================================================

Implementation of a Tabu Search algorithm for solving Vehicle Routing
Problem with Cross- Docking (VRPCD).

For information about Tabu Search algorithm see below:
https://en.wikipedia.org/wiki/Tabu_search

For information about VRPCD see below:
http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.87.1498

"""

import math

__author__ = 'Thodoris Sotiropoulos'


def tabu_search_vrpcd(G, cross_dock, Q, T, load, px='px', py='py',
                      node_type='type', quantity='quantity', pair='pair',
                      cons='consolidation', tolerance=100, L=12, k=2):
    """
    Implementation of a Tabu Search algorithm for solving VRPCD.
    Algorithm uses a 2.5 opt move to generate neighbor solution at every
    iteration.

    :param G:, graph, a graph with nodes which represents suppliers and customers.
    :param cross_dock: label of node representing cross-dock.
    :param Q: int, maximum capacity for every vehicle.
    :param T: int, maximum route duration.
    :param load: int, constant time of load's preparation.
    :param px: string, optional (default='px').
    Node data key corresponding to the position of node to x-axis.
    :param py: string, optional (default='py').
    Node data key corresponding to the position of node to y-axis.
    :param node_type: string, optional (default='type').
    Node data key corresponding to the type of node (supplier or customer).
    :param quantity: string, optional (default='quantity').
    Node data key corresponding to the quantity of node.
    :param pair: string, optional (default='pair').
    Node data key corresponding to the label of node's pair.
    :param cons: string, optional (default='consolidation')
    Node data key corresponding to the number of loads/unloads every vehicle
    does at the cross-dock.
    :param tolerance: int, optional (default=20)
    Number of consecutive iterations of algorithm loop that cost of best
    solution is not decreased.
    :param L: int, optional (default=12)
    Maximum length of tabu list.
    :param k:, int, optional (default=2)
    Constant representing the significance to use less vehicles.
    :return: G, graph representing the best solution, float, cost of best
    solution.
    """
    dist = lambda u, v: (math.sqrt((G.node[u][px] - G.node[v][px]) ** 2 + (
        G.node[u][py] - G.node[v][py]) ** 2.0))

    if G.number_of_edges() == 0:
        vehicle_id = 0
        capacity = {}
        duration = {}
        for u in G:
            if G.node[u][node_type] == 'pickup':
                capacity[vehicle_id] = [G.node[u][quantity],
                                        G.node[u][quantity]]
                duration[vehicle_id] = 2 * dist(cross_dock, u) + 2 * dist(
                    cross_dock, G.node[u][pair])
                vehicle_id += 1

        # Construct an initial-primitive solution.
        G = construct_primitive_solution(G, cross_dock, node_type, pair, dist)
        G = clarke_wright(G, cross_dock, dist, duration, capacity, node_type,
                        quantity, pair, Q, T)
        used_vehicles = sum(1 for vehicle in capacity.keys() if
                            capacity[vehicle][0] != 0 and capacity[vehicle][1] != 0)
    cost = sum(duration.values()) + used_vehicles ** k
    best_cost = float('inf')
    best_sol = None
    tabu_list = []
    max_iter = 0

    # Tabu algorithm's main loop.
    while max_iter < tolerance:
        try:
            cost = _apply_move(G, cross_dock, cost, best_cost, dist,
                               tabu_list, capacity, duration, load, Q, T, L,
                               node_type, quantity, cons, pair, used_vehicles,
                               k)
        except IOError:
            break
        if cost < best_cost:
            max_iter = 0
            best_cost = cost
            best_sol = G.copy()
        else:
            max_iter += 1

    return best_sol, best_cost, capacity, duration


def clarke_wright(G, cross_dock, dist, duration, capacity, node_type,
                  quantity, pair, Q, T):
    """
    Implementation of Clarke-Wright algorithm to construct an initial solution
    for Tabu Search local search algorithm.

    :param G:, graph, a graph representing the initial-primitive solution needed
    for this algorithm.
    :param cross_dock: label of node representing cross-dock.
    :param dist: array of distance between every pair of nodes.
    :param duration, dict, duration of every vehicle's route.
    :param capacity, array of used capacity of every vehicle.
    :param node_type: string, optional (default='type').
    Node data key corresponding to the type of node (supplier or customer).
    :param quantity: string, optional (default='quantity').
    Node data key corresponding to the quantity of node.
    :param pair: string, optional (default='pair').
    :param Q: int, maximum capacity for every vehicle.
    :param T: int, maximum route duration.
    :return: G, graph representing the initial solution for Tabu Search.
    """
    savings = {}
    for u in G:
        for v in G:
            if u != cross_dock and v != cross_dock and u != v \
                    and G.node[u][node_type] == G.node[v][node_type]:
                pair_u = G.node[u][pair]
                pair_v = G.node[v][pair]
                savings[(u, v)] = dist(cross_dock, u) + dist(cross_dock, v) \
                                  - dist(u, v) + dist(cross_dock, pair_u) \
                                  + dist(cross_dock, pair_v) - dist(pair_u,
                                                                    pair_v)

    import operator

    sorted_savings = sorted(savings.items(), key=operator.itemgetter(1),
                            reverse=True)

    for (u, v), savings in sorted_savings:
        pair_u = G.node[u][pair]
        pair_v = G.node[v][pair]
        if (G.has_edge(cross_dock, u)) and (G.has_edge(v, cross_dock))\
                and (not G.has_edge(u, v)):
            index, pair_index = (0, 1) if G.node[v][node_type] == 'pickup' else (1, 0)

            # Check if vehicle's capacity is sufficient to service new node.
            if capacity[G.node[v]['vehicle']][index] + G.node[u][quantity] > Q \
                    and capacity[G.node[pair_v]['vehicle']][pair_index] + \
                            G.node[u][quantity] > Q:
                continue

            # Initialize object to check the feasibility of solution
            # based on its duration.
            if duration[G.node[v]['vehicle']] + (
                    dist(pair_u, pair_v) + dist(u, v) - dist(cross_dock, v)) - (
                    dist(cross_dock, pair_v) + dist(cross_dock, u)
                    + dist(cross_dock, pair_u)) > T:
                continue

            if not G.has_edge(u, cross_dock) or not G.has_edge(cross_dock, v) \
                    or not G.has_edge(pair_u, cross_dock) or not G.has_edge(
                    cross_dock, pair_v):
                continue

            # Update capacity of vehicles.
            capacity[G.node[v]['vehicle']][index] += G.node[u][quantity]
            capacity[G.node[v]['vehicle']][pair_index] += G.node[u][quantity]
            capacity[G.node[u]['vehicle']][index] -= G.node[u][quantity]
            capacity[G.node[u]['vehicle']][pair_index] -= G.node[u][quantity]

            # Update duration of vehicles.
            duration[G.node[v]['vehicle']] += (
                dist(pair_u, pair_v) + dist(u, v) - dist(cross_dock, v) - dist(
                cross_dock, pair_v) + dist(cross_dock, u) + dist(cross_dock, pair_u))
            duration[G.node[u]['vehicle']] -= (
                dist(u, cross_dock) + dist(cross_dock, u)) + (
                dist(pair_u, cross_dock) + dist(cross_dock, pair_u))

            # Update solution.
            G.node[u]['vehicle'] = G.node[v]['vehicle']
            G.node[pair_u]['vehicle'] = G.node[v]['vehicle']
            G.remove_edge(u, cross_dock)
            G.remove_edge(cross_dock, v)
            G.edge[cross_dock][u]['vehicle'] = G.node[v]['vehicle']
            G.edge[cross_dock][pair_u]['vehicle'] = G.node[v]['vehicle']
            G.add_edge(u, v, weight=dist(u, v), type=G.node[v][node_type],
                       vehicle=G.node[v]['vehicle'])
            G.remove_edge(pair_u, cross_dock)
            G.remove_edge(cross_dock, pair_v)
            G.add_edge(pair_u, pair_v, weight=dist(pair_u, pair_v),
                       type=G.node[pair_v][node_type], vehicle=G.node[v]['vehicle'])
    return G


def construct_primitive_solution(G, cross_dock, node_type, pair, dist):
    """
    Constructs a fist solution for Clarke-Wright algorithm.

    This solution uses for every pair of nodes (supplier-customer) a different
    vehicle to serve them.

    :param G: graph, a graph with nodes which represents suppliers and customers.
    :param cross_dock: label of node which represents cross-dock.
    :param node_type: string , node data key corresponding to the type of node
    (supplier or customer).
    :param pair: string, node data key corresponding to the label of node's pair.
    :param dist: array of distance between every pair of nodes.
    """
    vehicle_id = 0
    for u in G:
        if G.node[u][node_type] == 'pickup':
            pair_node = G.node[u][pair]
            G.node[u]['vehicle'] = vehicle_id
            G.node[pair_node]['vehicle'] = vehicle_id
            G.add_edge(cross_dock, u, weight=dist(cross_dock, u),
                       type=G.node[u][node_type], vehicle=vehicle_id)
            G.add_edge(u, cross_dock, weight=dist(u, cross_dock),
                       type=G.node[u][node_type], vehicle=vehicle_id)
            G.add_edge(cross_dock, pair_node,
                       weight=dist(cross_dock, pair_node),
                       type=G.node[pair_node][node_type], vehicle=vehicle_id)
            G.add_edge(pair_node, cross_dock,
                       weight=dist(pair_node, cross_dock),
                       type=G.node[pair_node][node_type], vehicle=vehicle_id)
            vehicle_id += 1
    return G


def _apply_move(G, cross_dock, cost, min_cost, dist, tabu_list, capacity,
                duration, load, Q, T, L, node_type, quantity, cons, pair, used,
                k):
    """
    Apply a move to a current solution to generate neighbor solutions.

    Use a 1-0 exchange move. It moves an element of solution in a new position.
    For example if we apply 1-0 exchange in the solution  = [3, 2, 1, 4, 3]
    we can get the following the transfer of 4 element to the second position.
    A' = [3, 4, 2, 1, 3].

    Then keep the best feasible neighbor solution to be current solution to the
    next iteration of algorithm's loop.

    :param G, graph which represents current solution.
    :param min_cost, float, overall minimum cost.
    :param dist, array of distance between any pair of nodes.
    :param tabu_list, list of characteristics which are tabu in order to be
    excluded.
    :param capacity, array of used capacity of every vehicle.
    :param duration, dict, duration of every vehicle's route.
    :param load, int, constant time of load's preparation.
    :param Q, int, maximum vehicle capacity.
    :param T, int maximum route duration for every vehicle.
    :param L, int, maximum length of tabu list.
    :param node_type: string, node data key corresponding to the type of node
    (supplier or customer).
    :param quantity: string, node data key corresponding to the quantity of node.
    :param pair: string, node data key corresponding to the label of node's pair.
    :param cons: string, node data key corresponding to the number of
    loads/unloads every vehicle does at the cross-dock.
    :param used: int, number of used vehicles on current solution.
    :param k: int, constant representing the significance to use less vehicles.
    :return Cost of best neighbor solution (which is now current solution for
    next iteration of the algorithm.
    :raise IOError when no feasible neighbor solution can be found.
    """
    min_neighbor_cost = float('inf')
    min_neighbor_move = ()
    for n, node_data in G.nodes(data=True):

        # Ignore node which represents cross-dock.
        if n == cross_dock:
            continue

        for u, v, data in G.edges(data=True):
            if data[node_type] != node_data[node_type]:
                continue

            if n == u or n == v:
                continue

            index = 0 if data[node_type] == 'pickup' else 1

            # Check if vehicle's capacity is sufficient to service new node.
            if capacity[data['vehicle']][index] + node_data[quantity] > Q:
                continue

            # Initialize object to check the feasibility of solution
            # based on its duration.
            checker = TimeChecker(G, cross_dock, data['vehicle'],
                                  node_data['vehicle'],
                                  G.node[node_data[pair]]['vehicle'], n,
                                  duration, node_type, cons)

            # Check if vehicle's total duration does not surpress T.
            if not checker.is_feasible(dist, u, v, node_data[quantity],
                                       load, T):
                continue
            withdrawn = True if sum(capacity[node_data['vehicle']]) \
                                - node_data[quantity] == 0 else False
            used_vehicles = used - 1 if withdrawn else used
            # Calculate overall cost of neighbor's solution.
            neighbor_cost = (checker.new_vehicle_dur - duration[
                checker.vehicle]) + (
                checker.previous_vehicle_dur - duration[
                    checker.previous_vehicle]) + (
                                checker.pair_vehicle_dur - duration[
                                    checker.pair_vehicle]) + used_vehicles ** k

            if neighbor_cost < min_neighbor_cost:
                # Check if best (so far neighbor solution) is in the tabu list.
                if any((u, (n, v)) == tabu[0] or (v, (u, n)) == tabu[1] or (
                        n, (u, v)) == tabu[2]
                       for tabu in tabu_list) and (
                                neighbor_cost + cost >= min_cost):
                    continue

                min_neighbor_move = (n, u, v, node_data[node_type],
                                     data['vehicle'])
                min_neighbor_cost = neighbor_cost
                best = index
                vehicles = checker.vehicle, checker.previous_vehicle, \
                           checker.pair_vehicle
                consolidations = checker.consolidation_vehicle, \
                                 checker.consolidation_previous, \
                                 checker.consolidation_pair
                new_durations = checker.new_vehicle_dur, \
                                checker.previous_vehicle_dur, \
                                checker.pair_vehicle_dur
                withdrawn_vehicle = withdrawn

    if min_neighbor_move == ():
        raise IOError('No any feasible neighbor solution found.')

    # Adapt graph based on the current solution which is the best neighbor
    # solution.
    update_edges(G, min_neighbor_move, dist)
    update_node_data(G, cross_dock, min_neighbor_move[0], vehicles,
                     new_durations,
                     consolidations, capacity, duration, best, node_type,
                     quantity, cons)
    update_tabu_list(min_neighbor_move, tabu_list, L)
    used -= 1 if withdrawn_vehicle else 0
    return sum(duration.values()) + used ** k


def update_tabu_list(neighbor_move, tabu_list, L):
    """
    Adds characteristics of best neighbor solution to the tabu list.

    :param neighbor_move: tuple, characteristics of neighbor solution.
    :param tabu_list: list, current tabu list
    :param L: int, maximum length of tabu list.
    """
    n, u, v, node_type, vehicle = neighbor_move
    if len(tabu_list) >= L:
        tabu_list.pop(0)
    tabu_list.append(((u, (n, v)), (v, (u, n)), (n, (v, u))))


def update_node_data(G, cross_dock, n, vehicles, new_durations,
                     new_consolidations, capacity, duration, best, node_type,
                     quantity, cons):
    """
    Updates data of node to be compatible with the best neighbor solution
    (current solution for next iteration).

    :param G: graph, current solution of VRPCD.
    :param n: int, node's label which was moved.
    :param vehicles: tuple, (current vehicle which serves moved node,
    previous vehicle which served moved node, vehicle which serves pair of moved
    node).
    :param new_durations: tuple, new duration's times for current, previous and
    pair vehicles.
    :param new_consolidations: tuple, number of unloads/reloads for current,
    previous, and par vehicles.
    :param capacity: dict, used capacity for every vehicle.
    :param duration: dict, duration of every vehicle's route.
    :param best: int, index which defines if there was an addition/deletion to
    the number of unloads or loads which is defined by the best neighbor
    solution.
    :param node_type: string, node data key corresponding to the type of node
    (supplier or customer).
    :param quantity: string, node data key corresponding to the quantity of node.
    :param cons: string, node data key corresponding to the number of
    loads/unloads every vehicle does at the cross-dock.
    """
    if G.node[n]['vehicle'] != vehicles[0]:
        index = 0 if G.node[n][node_type] == 'pickup' else 1
        capacity[G.node[n]['vehicle']][index] -= G.node[n][quantity]
        G.node[n]['vehicle'] = vehicles[0]
        capacity[vehicles[0]][index] += G.node[n][quantity]
    for i in range(3):
        duration[vehicles[i]] = new_durations[i]
        pair_index = 0 if best == 1 else 1
        G.node[cross_dock][cons][vehicles[0]][best] = new_consolidations[0]
        G.node[cross_dock][cons][vehicles[1]][best] = new_consolidations[1]
        G.node[cross_dock][cons][vehicles[2]][pair_index] = new_consolidations[
            2]


def update_edges(G, best_neighbor, dist):
    """
    Updates edges (edges removal/addition) to be compatible with the best
     neighbor solution (current solution for next iteration).

    :param G: graph, current solution of VRPCD.
    :param neighbor_move: tuple, characteristics of neighbor solution.
    :param dist: array of distance between any pair of nodes.
    """
    n, u, v, node_type, vehicle = best_neighbor
    pred = G.predecessors(n)[0]
    succ = G.successors(n)[0]
    G.remove_edge(u, v)
    G.remove_edge(n, succ)
    G.remove_edge(pred, n)
    if pred != succ:
        G.add_edge(pred, succ, weight=dist(pred, succ), type=node_type,
                   vehicle=G.node[n]['vehicle'])
    G.add_edge(u, n, weight=dist(u, n), type=node_type, vehicle=vehicle)
    G.add_edge(n, v, weight=dist(n, v), type=node_type, vehicle=vehicle)


class TimeChecker:
    """
    This class is a checker of feasibility of one solution based on the duration
    of every vehcile
    """

    def __init__(self, G, cross_dock, vehicle, previous_vehicle, pair_vehicle,
                 node, duration, node_type, cons):
        self.G = G
        self.vehicle = vehicle
        self.node = node
        self.duration = duration
        self.pair_vehicle = pair_vehicle
        self.previous_vehicle = previous_vehicle
        self.new_vehicle_dur = 0
        self.previous_vehicle_dur = 0
        self.pair_vehicle_dur = 0
        self.type = G.node[self.node][node_type]
        self.consolidation = self.G.node[cross_dock][cons]

    def is_feasible(self, dist, u, v, quantity, load, T):
        """
        Checks if solution is feasible by checking that no vehicle's duration
         surpresses the maximum route's duration.

        :param dist, array of distance between any pair of nodes.
        :param u, label of source node of edge that node is moving.
        :param v, label of target node of edge that node is moving.
        :param quantity, int, quantity that has to be picked up from a supplier
         or delivered to a customer.
        :param load, int, constant time of load's preparation.
        :param T, int, maximum route duration for every vehicle.
        :return True if solution is feasible; False otherwise.
        """
        self.calculate_new_durations(dist, u, v,
                                     self.calculate_consolidation(load,
                                                                  quantity))
        if self.new_vehicle_dur > T:
            return False
        if self.previous_vehicle_dur > T:
            return False
        if self.pair_vehicle_dur > T:
            return False
        return True

    def calculate_new_durations(self, dist, u, v, cd_time):
        """
        This method sets new total duration of every vehicle's trip (vehicle
        which now serves node, vehicle which previously served node, vehicle
        which serves pair of node).

        :param dist, array, distance from every pair of nodes.
        :param u, label of source node of edge that node is moving.
        :param v, label of target node of edge that node is moving.
        :param quantity, int, quantity that has to be picked up from a supplier
         or delivered to a customer.
        :param cd_time Time spending at cross-dock for every vehicle.
        """
        pred = self.G.predecessors(self.node)[0]
        succ = self.G.successors(self.node)[0]
        v_cd, previous_cd, pair_cd = cd_time
        self.new_vehicle_dur = self.duration[self.vehicle] + dist(u, self.node) \
                               + dist(self.node, v) - dist(u, v) + v_cd
        self.previous_vehicle_dur = self.duration[self.previous_vehicle] \
                                    + dist(pred, succ) - dist(pred, self.node) \
                                    - dist(self.node, succ) + previous_cd
        self.pair_vehicle_dur = self.duration[self.pair_vehicle] + pair_cd
        if self.vehicle == self.previous_vehicle and self.vehicle != self.pair_vehicle:
            self.new_vehicle_dur = self.duration[self.vehicle] + dist(
                u, self.node) + dist(self.node, v) - dist(u, v) + v_cd \
                                   - dist(pred, self.node) - dist(
                self.node, succ) + dist(pred, succ)
            self.previous_vehicle_dur = self.new_vehicle_dur
        elif self.previous_vehicle == self.pair_vehicle \
                and self.previous_vehicle != self.vehicle:
            self.pair_vehicle_dur = self.previous_vehicle_dur
        elif self.vehicle == self.pair_vehicle \
                and self.vehicle != self.previous_vehicle:
            self.pair_vehicle_dur = self.new_vehicle_dur
        elif self.vehicle == self.pair_vehicle == self.previous_vehicle:
            self.pair_vehicle_dur = self.new_vehicle_dur
            self.previous_vehicle_dur = self.new_vehicle_dur

    def calculate_consolidation(self, load, quantity):
        """
        This method calculates new times spending at cross-dock for
        every vehicle (vehicle which now serves node, vehicle which previously
        served node, vehicle which serves pair of node).

        :param load, int, constant time of load's preparation.
        :param quantity, int, quantity that has to be picked up from a supplier
         or delivered to a customer.
        :return tuple, new times spending at cross-dock for every
        vehicle.
        """
        self.init_cd_consolidations()
        (cd_vehicle, n_cd_vehicle), (cd_previous, n_cd_previous), cd_pair \
            = self.init_cd_diff_times(load, quantity)
        vehicle_at_cd = 0
        pair_vehicle_at_cd = 0
        previous_vehicle_at_cd = 0
        if self.previous_vehicle != self.vehicle:
            if self.previous_vehicle == self.pair_vehicle:
                previous_vehicle_at_cd = cd_previous
                self.consolidation_pair += 1
            else:
                self.consolidation_previous -= 1
                previous_vehicle_at_cd = n_cd_previous
        if self.vehicle == self.pair_vehicle:
            if self.vehicle != self.previous_vehicle:
                self.consolidation_pair -= 1
                vehicle_at_cd = n_cd_vehicle
        else:
            if self.vehicle != self.previous_vehicle:
                self.consolidation_vehicle += 1
                vehicle_at_cd = cd_vehicle
            if self.pair_vehicle == self.previous_vehicle:
                pair_vehicle_at_cd = cd_pair
        return vehicle_at_cd, previous_vehicle_at_cd, pair_vehicle_at_cd

    def init_cd_consolidations(self):
        """Initialize number of loads/unloads at the cross-dock for every
        vehicle."""
        index, pair_index = (0, 1) if self.type == 'pickup' else (1, 0)
        self.consolidation_vehicle = self.consolidation[self.vehicle][index]
        self.consolidation_pair = self.consolidation[self.pair_vehicle][
            pair_index]
        self.consolidation_previous = self.consolidation[self.previous_vehicle][
            index]

    def init_cd_diff_times(self, load, quantity):
        """
        This method calculates the variation of time spending at cross-dock for
        every vehicle (vehicle which now serves node, vehicle which previously
        served node, vehicle which serves pair of node).

        :param load, int, constant time of load's preparation.
        :param quantity, int, quantity that has to be picked up from a supplier
         or delivered to a customer.
        :return tuple, variation of time spending at cross-dock for every
        vehicle.
        """
        diff1 = load + quantity
        diff2 = quantity
        cd_vehicle = diff1 if self.consolidation_vehicle == 0 else diff2
        negative_cd_vehicle = -diff1 if self.consolidation_vehicle == 1 else -diff2
        cd_previous = diff1 if self.consolidation_previous == 0 else diff2
        negative_cd_previous = -diff1 if self.consolidation_previous == 1 else -diff2
        cd_pair = diff1 if self.consolidation_pair == 0 else diff2
        return (cd_vehicle, negative_cd_vehicle), \
               (cd_previous, negative_cd_previous), cd_pair
