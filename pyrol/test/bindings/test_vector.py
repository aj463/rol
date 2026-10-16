import unittest

import numpy as np

import pyrol
from pyrol import Vector


class CountingVector(Vector):
    """Minimal vector that relies on the default ROL::Vector::axpy."""

    clones = 0

    def __init__(self, values):
        self.values = np.array(values, dtype=float)
        super().__init__()

    def plus(self, x):
        self.values += x.values

    def scale(self, alpha):
        self.values *= alpha

    def dot(self, x):
        return float(np.dot(self.values, x.values))

    def norm(self):
        return float(np.sqrt(self.dot(self)))

    def clone(self):
        CountingVector.clones += 1
        return CountingVector(np.zeros_like(self.values))


class TestVectorClone(unittest.TestCase):

    def test_cpp_clone_calls_python_clone_once(self):
        y = CountingVector([1.0, 2.0])
        x = CountingVector([3.0, 4.0])
        CountingVector.clones = 0
        y.axpy(2.0, x)  # ROL::Vector::axpy clones x exactly once
        self.assertEqual(CountingVector.clones, 1)
        np.testing.assert_allclose(y.values, [7.0, 10.0])


class TestVectorExports(unittest.TestCase):

    def test_star_import(self):
        namespace = {}
        exec("from pyrol.vectors import *", namespace)
        self.assertIn("NumPyVector", namespace)

    def test_forScalar(self):
        self.assertIs(pyrol.Vector_forScalar(), pyrol.Vector)
        self.assertIs(pyrol.Vector_forScalar("double"), pyrol.Vector)


if __name__ == "__main__":
    unittest.main()
