import importlib.util
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

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

    def test_compare_per_element_tolerance(self):
        want = [4e6, 0.5]
        self.assertTrue(compare([4e6, 0.5 + 1e-9], want, [1e-12, 1e-8])[0])
        self.assertFalse(compare([4e6, 0.5 + 1e-9], want, [1e-12, 1e-12])[0])
        self.assertFalse(compare(list(want), want, [1e-12])[0])


class XfailGateTest(unittest.TestCase):
    """--xfail-r-fixes only excuses a value mismatch on an r_fix case."""

    def gate(self, run_case_result):
        case = {"name": "c", "func": "f", "expected": [1.0, 2.0], "r_fix": "X"}
        with tempfile.TemporaryDirectory() as tmp:
            truth = Path(tmp) / "truth.json"
            truth.write_text(json.dumps({"cases": [case]}), encoding="utf-8")
            argv = ["run_fidelity.py", "--truth", str(truth), "--eval-bin", str(truth),
                    "--report", str(Path(tmp) / "report.md"), "--xfail-r-fixes"]
            with mock.patch.object(sys, "argv", argv), \
                    mock.patch.object(run_fidelity, "run_case", return_value=run_case_result), \
                    mock.patch("builtins.print"):
                return run_fidelity.main()

    def test_value_mismatch_is_xfail(self):
        self.assertEqual(self.gate(([1.0, 3.0], "")), 0)

    def test_evaluator_failure_fails(self):
        self.assertEqual(self.gate((None, "terminate called")), 1)

    def test_length_mismatch_is_xfail(self):
        self.assertEqual(self.gate(([1.0], "")), 0)

    def test_unparsable_output_fails(self):
        with mock.patch.object(run_fidelity.subprocess, "run") as run:
            run.return_value.returncode = 0
            run.return_value.stdout = "1,oops"
            got, err = run_fidelity.run_case("nns_eval", {"func": "f"})
        self.assertIsNone(got)
        self.assertIn("unparsable", err)
        self.assertEqual(self.gate((got, err)), 1)


class CliGateTest(unittest.TestCase):
    """The gate run end to end, with a Python script standing in for nns_eval."""

    def run_gate(self, evaluator):
        with tempfile.TemporaryDirectory() as tmp:
            script = Path(tmp) / "eval.py"
            script.write_text(evaluator, encoding="utf-8")
            # run_case calls [eval_bin, func, "--stdin"], so func is the script.
            case = {"name": "c", "func": str(script), "expected": [1.0, 2.0], "r_fix": "X"}
            truth = Path(tmp) / "truth.json"
            truth.write_text(json.dumps({"cases": [case]}), encoding="utf-8")
            proc = subprocess.run(
                [sys.executable, str(_spec.origin), "--truth", str(truth),
                 "--eval-bin", sys.executable, "--report", str(Path(tmp) / "report.md"),
                 "--xfail-r-fixes"], capture_output=True, text=True)
            return proc.returncode, proc.stdout

    def test_value_mismatch_is_xfail(self):
        code, out = self.run_gate("print('1,3')\n")
        self.assertEqual(code, 0, out)
        self.assertIn("XFAIL", out)

    def test_empty_output_fails(self):
        code, out = self.run_gate("")
        self.assertEqual(code, 1, out)
        self.assertIn("no output", out)

    def test_nonzero_exit_fails(self):
        code, out = self.run_gate("raise SystemExit(3)\n")
        self.assertEqual(code, 1, out)


if __name__ == "__main__":
    unittest.main()
