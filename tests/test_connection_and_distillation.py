import unittest
from ..connection_and_distillation import get_bipartitions, iterate_over_paths_of_length_k
from ..network_functions import create_repeater_chain, find_candidate_paths
from ..scheme_class import Scheme
assertions = unittest.TestCase('__init__')


class MyTest(unittest.TestLoader):

    def test_iterate_over_schemes_of_subpaths(self):
        scheme1 = Scheme(('A'), None, None, None, [0.9, 0, 0, 0.1], 1, 1)
        scheme2 = Scheme(('A'), None, None, None, [1, 0, 0, 0], 1.1, 1)
        scheme3 = Scheme(('B'), None, None, None, [0.7, 0, 0, 0.1], 2, 1)
        scheme4 = Scheme(('B'), None, None, None, [0.8, 0, 0, 0], 6.5, 1)
        scheme5 = Scheme(('B'), None, None, None, [0.85, 0, 0, 0], 7.5, 1)

        add_schemes([scheme1, scheme2, scheme3, scheme4, scheme5])

        self.assertTrue(len(list(iterate_over_schemes_of_subpaths(('A'), ('B')))) == 6)

    def test_iterate_distillation_schemes(self):
        scheme1 = Scheme(('A'), None, None, None, [0.9, 0, 0, 0.1], 1, 1)
        scheme2 = Scheme(('A'), None, None, None, [0.95, 0, 0, 0], 1.1, 1)

        add_schemes([scheme1, scheme2])
        iterator = iterate_distillation_subschemes(
            ('A'), general_distillation=True)
        self.assertEqual(set(iterator), {[(scheme1, scheme1), (scheme1, scheme2),
                                          (scheme2, scheme1), (scheme2, scheme2)]})

        iterator = iterate_distillation_subschemes(
            ('A'), general_distillation=False)
        self.assertEqual(set(iterator), {[(scheme1, scheme1),
                                          (scheme2, scheme2)]})


if __name__ == '__main__':
    unittest.main()
