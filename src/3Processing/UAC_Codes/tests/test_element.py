#!/usr/bin/python3
#
import unittest
from uac.element import Element


class TestElement(unittest.TestCase):
    def test_init(self):
        e = Element(1, ["1", "2", 3.0, "4", 5], "11")
        self.assertEqual(e.type, 1)
        self.assertEqual(e.n, [1, 2, 3, 4, 5])
        self.assertEqual(e.reg, 11)

    def test_eq(self):
        e1 = Element(0, [0], 0)
        e2 = e1
        self.assertNotEqual(e1, Element(0, [0], 0))
        self.assertEqual(e1, e2)

    def test_region(self):
        e = Element(0, [0], 0)

        e.region(5)
        self.assertEqual(e.reg, 5)

        e.region(11.5)
        self.assertEqual(e.reg, 11)

        e.region("0.5")
        self.assertEqual(e.reg, 0)
