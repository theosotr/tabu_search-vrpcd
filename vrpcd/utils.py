import csv
import os
import networkx as nx


def _positive(value):
    return abs(int(value))

NODE_ID = 'ID'
X_AXIS = 'x_axis'
Y_AXIS = 'y_axis'
QUANTITY = 'quantity'
PAIR = 'pair'
NODE_TYPE = 'type'
CONSOLIDATION = 'consolidation'
VEHICLE = 'vehicle'

COLUMN_NAMES = (NODE_ID, X_AXIS, Y_AXIS, QUANTITY, PAIR)
COLUMN_TYPES = {
    NODE_ID: str,
    X_AXIS: int,
    Y_AXIS: int,
    QUANTITY: _positive,
    PAIR: str,
}
CROSS_DOCK = 0
SUPPLIER = 1
CUSTOMER = 2


def get_consolidation(nu_pairs):
    return {i: [0, 0] for i in range(nu_pairs)}


def get_node_type(quantity):
    quantity = int(quantity)
    if quantity == 0:
        return CROSS_DOCK
    elif quantity > 0:
        return SUPPLIER
    else:
        return CUSTOMER


def validate_data(data):
    parsed_data = []
    for i, column_name in enumerate(COLUMN_NAMES):
        column_value = data[i]
        column_type = COLUMN_TYPES.get(column_name)
        try:
            parsed_data.append(column_type(column_value))
        except TypeError:
            msg = 'Value of {!r} is not {!r} ({!r}).'
            raise TypeError(msg.format(column_value, column_type,
                                       column_value))
    return parsed_data


def _get_node_kwargs(node_type, node_data, nu_pairs):
    kwargs = {COLUMN_NAMES[i + 1]: node_attr
              for i, node_attr in enumerate(node_data)}
    kwargs.update({NODE_TYPE: node_type})
    if node_type == CROSS_DOCK:
        kwargs.update(
            {CONSOLIDATION: get_consolidation(nu_pairs)})
    return kwargs


def load_data(path, delimiter=' '):
    if not os.path.isfile(path):
        msg = 'Given path {!r} is not a file'
        raise IOError(msg.format(path))
    G = nx.DiGraph()
    with open(path, 'r') as fil:
        nodes = list(csv.reader(fil, delimiter=delimiter))
        assert len(nodes) % 2, (
            'The number of nodes must be an odd number. There must be a pair'
            ' of suppliers and customers and a cross dock.'
            ' Found: {!s}'.format(len(nodes)))
        for i, row in enumerate(nodes):
            assert len(row) == len(COLUMN_NAMES), (
                'Row {!r} is incomplete. Missing node data.'.format(i + 1))
            data = validate_data(row)
            node, node_data = data[0], data[1:]
            node_type = get_node_type(row[-2])
            nu_pairs = (len(nodes) - 1) / 2
            node_kwargs = _get_node_kwargs(node_type, node_data, nu_pairs)
            G.add_node(node, **node_kwargs)
    return G
