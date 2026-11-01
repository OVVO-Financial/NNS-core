import importlib.util
import math
import unittest
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "run_fidelity", Path(__file__).resolve().parents[2] / "scripts" / "run_fidelity.py")
run_fidelity = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(run_fidelity)
within_tol = run_fidelity.within_tol
compare = run_fidelity.compare

INF = math.inf
NAN = math.nan


class WithinTolTest(unittest.TestCase):
    def test_close_values(self):
        self.assertTrue(within_tol(1.0, 1.0 + 1e-12, 1e-9))
        self.assertTrue(within_tol(1e6, 1e6 * (1 + 1e-10), 1e-9))
        self.assertFalse(within_tol(1.0, 1.1, 1e-9))

    def test_nan_matches_only_nan(self):
        self.assertTrue(within_tol(NAN, NAN, 1e-9))
        self.assertFalse(within_tol(NAN, 1.0, 1e-9))
        self.assertFalse(within_tol(1.0, NAN, 1e-9))

    def test_infinity_matches_only_the_same_infinity(self):
        self.assertTrue(within_tol(INF, INF, 1e-9))
        self.assertTrue(within_tol(-INF, -INF, 1e-9))
        self.assertFalse(within_tol(INF, 1.0, 1e-9))
        self.assertFalse(within_tol(1.0, INF, 1e-9))
        self.assertFalse(within_tol(INF, -INF, 1e-9))
        self.assertFalse(within_tol(-INF, 1e308, 1e-9))

    def test_compare_rejects_any_infinite_mismatch(self):
        want = [0.5, 1.0, 2.0]
        self.assertTrue(compare(list(want), want, 1e-9)[0])
        self.assertFalse(compare([INF] * 3, want, 1e-9)[0])
        self.assertFalse(compare([0.5, 1.0], want, 1e-9)[0])


if __name__ == "__main__":
    unittest.main()
