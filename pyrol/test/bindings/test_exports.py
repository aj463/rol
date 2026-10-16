import unittest

import pyrol
import pyrol.unsupported


class TestExports(unittest.TestCase):

    def test_supported_classes(self):
        for name in ["Bounds", "Constraint", "LinearConstraint", "LinearOperator",
                     "Objective", "Problem", "Solver", "Vector"]:
            with self.subTest(name=name):
                cls = getattr(pyrol, name)
                forScalar = getattr(pyrol, name + "_forScalar")
                self.assertIs(forScalar(), cls)
                self.assertIs(forScalar("double"), cls)

    def test_supported_functions(self):
        for name in ["getCout", "getParametersFromXmlFile", "getParametersFromYamlFile"]:
            with self.subTest(name=name):
                self.assertTrue(callable(getattr(pyrol, name)))

    def test_unsupported_classes(self):
        for name in ["Constraint_SimOpt", "Objective_SimOpt", "PartitionedVector",
                     "Reduced_Objective_SimOpt", "Vector_SimOpt"]:
            with self.subTest(name=name):
                self.assertTrue(hasattr(pyrol.unsupported, name))
                self.assertFalse(hasattr(pyrol, name))

    def test_vectors_star_import(self):
        namespace = {}
        exec("from pyrol.vectors import *", namespace)
        self.assertIn("NumPyVector", namespace)

    def test_version(self):
        self.assertTrue(pyrol.version().startswith("PyROL version"))


if __name__ == "__main__":
    unittest.main()
