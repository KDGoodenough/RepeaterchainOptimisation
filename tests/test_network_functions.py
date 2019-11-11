import numpy as np
import networkx as nx
import global_file

from scheme_class import *
import collections
import copy
import operator
import warnings
from network_functions import *


import unittest
assertions = unittest.TestCase('__init__')


class MyTest(unittest.TestCase):

    # network tests

    def test_create_repeater_chain(self):
        G = create_repeater_chain(50, 3)
        self.assertEqual(set(G.nodes()), set(('A', 'AB1', 'AB2', 'AB3', 'B')))
        for u, v in G.nodes(data=True):
            if global_file.pert == 0:
                self.assertTrue(v['pos'][1] in [0, 12.5, 25.0, 37.5, 50])
            else:
                print("Perturbation of nodes is non-zero, pert = ", pert)
                self.assertTrue(True)

    def test_get_ordered_nodes(self):
        G = create_repeater_chain(5, 4)

        nodes = get_ordered_nodes(G, ('A', 'AB1', 'AB2', 'AB3', 'AB4', 'B'),
                                  ('A', 'AB1', 'AB2'), ('AB3', 'AB4', 'B'))
        positions = [node['pos'][1] for node in nodes]
        self.assertEqual(positions, [0, 2, 3, 5])

        nodes = get_ordered_nodes(G, ('A', 'AB1', 'AB2', 'AB3', 'AB4', 'B'),
                                  ('AB3', 'AB4', 'B'), ('A', 'AB1', 'AB2'))
        positions = [node['pos'][1] for node in nodes]
        self.assertEqual(positions, [5, 3, 2, 0])

        nodes = get_ordered_nodes(G, ('A', 'AB1', 'AB2', 'AB3', 'AB4'),
                                  ('A', 'AB1', 'AB2'), ('AB2', 'AB3', 'AB4'))
        positions = [node['pos'][1] for node in nodes]
        self.assertEqual(positions, [0, 2, 4])

        nodes = get_ordered_nodes(G, ('A', 'AB1', 'AB2', 'AB3', 'AB4'),
                                  ('AB2', 'AB3', 'AB4'), ('A', 'AB1', 'AB2'))
        positions = [node['pos'][1] for node in nodes]
        self.assertEqual(positions, [4, 2, 0])
        #           self.assertEqual(node['pos'], (0, i*10))


if __name__ == '__main__':
    unittest.main()
