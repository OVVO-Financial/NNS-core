"""Compare NNS-core C++ outputs against live R (Rcpp) truth.

Reads a truth file produced by tests/fidelity/generate_truth.R (which calls the
installed R NNS package), drives the portable C++ core via the `nns_eval` CLI on
the same inputs, and asserts every value matches within tolerance.

`expected` may be a scalar or a list (vectors and column-major matrices); every
element is compared. A case may carry its own `tol`, a scalar or one per
element. Cases tagged `r_fix` cover an R bug that NNS-core fixes ahead of R;
with --xfail-r-fixes a value mismatch on such a case is reported as XFAIL
(expected against an R build without the fix) and does not fail the run, and a
pass is reported as XPASS. An evaluator failure or unparsable output always
fails.

Exit code 0 only if every case passes. The sync workflow uses this gate to
decide whether a sync PR may be marked ready: no pass, no ready PR.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
from pathlib import Path

# Keys that are not function inputs and must not be forwarded as flags.
META_KEYS = {"name", "func", "expected", "tolerance", "tol", "r_fix"}


def repr_num(v) -> str:
    # Full-precision, plain decimal (avoid surprises from str() on ints/floats).
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, str):
        return v
    return repr(float(v))


def case_to_stdin(case: dict) -> str:
    lines = []
    for key, value in case.items():
        if key in META_KEYS:
            continue
        if isinstance(value, list):
            text = ",".join(repr_num(v) for v in value)
        else:
            text = repr_num(value)
        lines.append(f"{key}\t{text}")
    return "\n".join(lines) + "\n"


def as_list(v) -> list[float]:
    return [float(x) for x in v] if isinstance(v, list) else [float(v)]


def within_tol(got: float, want: float, tol: float) -> bool:
    if math.isnan(got) or math.isnan(want):
        return math.isnan(got) and math.isnan(want)
    if got == want:  # also covers matching infinities
        return True
    if math.isinf(got) or math.isinf(want):
        return False
    diff = abs(got - want)
    return diff <= tol or diff <= tol * max(abs(got), abs(want))


def compare(got: list[float], want: list[float],
            tol: float | list[float]) -> tuple[bool, float, str]:
    if len(got) != len(want):
        return False, float("nan"), f"length {len(got)} != {len(want)}"
    tols = tol if isinstance(tol, list) else [tol] * len(want)
    if len(tols) != len(want):
        return False, float("nan"), f"tol length {len(tols)} != {len(want)}"
    worst = 0.0
    ok = True
    for g, w, t in zip(got, want, tols):
        if not within_tol(g, w, t):
            ok = False
        if not (math.isnan(g) or math.isnan(w) or g == w):
            worst = max(worst, abs(g - w))
    return ok, worst, ""


def run_case(eval_bin: str, case: dict) -> tuple[list[float] | None, str]:
    proc = subprocess.run([eval_bin, case["func"], "--stdin"], input=case_to_stdin(case),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        return None, proc.stderr.strip()
    text = proc.stdout.strip()
    try:
        return ([float(t) for t in text.split(",")] if text else []), ""
    except ValueError:
        return None, f"unparsable output {text[:60]!r}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth", required=True)
    parser.add_argument("--eval-bin", required=True)
    parser.add_argument("--report", default="tests/fidelity/fidelity_report.md")
    parser.add_argument("--tol", type=float, default=None,
                        help="override tolerance from the truth file (and per-case tol)")
    parser.add_argument("--xfail-r-fixes", action="store_true",
                        help="report failures of r_fix-tagged cases as XFAIL (R build "
                             "without the fix) instead of failing")
    args = parser.parse_args()

    truth_path = Path(args.truth)
    if not truth_path.is_file():
        raise SystemExit(f"missing truth file: {truth_path}")
    if not Path(args.eval_bin).is_file():
        raise SystemExit(f"missing eval binary: {args.eval_bin}")

    truth = json.loads(truth_path.read_text(encoding="utf-8"))
    default_tol = float(truth.get("tolerance", 1e-9))
    cases = truth.get("cases", [])

    rows: list[tuple[str, str, int, float, str]] = []
    failures = xfail = xpass = 0
    for case in cases:
        tol = args.tol if args.tol is not None else case.get("tol", default_tol)
        tol = [float(t) for t in tol] if isinstance(tol, list) else float(tol)
        want = as_list(case["expected"])
        got, err = run_case(args.eval_bin, case)
        if got is None:
            ok, worst, note = False, float("nan"), err or "nns_eval failed"
        else:
            ok, worst, note = compare(got, want, tol)
        fix = case.get("r_fix")
        if ok:
            status = "XPASS" if (fix and args.xfail_r_fixes) else "ok"
            xpass += status == "XPASS"
        # A crash is never the missing R fix; any mismatch, even in length (a
        # fix can change how many seasonal periods come back), can be.
        elif fix and args.xfail_r_fixes and got is not None:
            status = "XFAIL"
            xfail += 1
        else:
            status = "FAIL"
            failures += 1
        if fix:
            note = (note + "; " if note else "") + f"r_fix {fix}"
        rows.append((case["name"], status, len(want), worst, note))

    lines: list[str] = []
    lines.append("# Live Rcpp fidelity report")
    lines.append("")
    lines.append(f"- R version: `{truth.get('r_version', '?')}`")
    lines.append(f"- R commit: `{truth.get('r_commit', '?')}`")
    lines.append(f"- default tolerance: `{default_tol:g}`")
    lines.append(f"- cases: {len(cases)}  |  failures: {failures}  |  "
                 f"xfail (R fix not in this R build): {xfail}  |  xpass: {xpass}")
    lines.append("")
    lines.append("| case | status | values | max abs diff | note |")
    lines.append("| --- | --- | --- | --- | --- |")
    for name, status, nvals, worst, note in rows:
        lines.append(f"| {name} | {status} | {nvals} | {worst:.3g} | {note} |")
    lines.append("")
    lines.append("## Result")
    lines.append("")
    if failures == 0 and cases:
        lines.append("All fidelity cases passed. Core matches live R within tolerance.")
        if xfail:
            lines.append("")
            lines.append(f"{xfail} case(s) tagged `r_fix` differ because this R build does "
                         "not contain the corresponding R fix yet (see "
                         "docs/r-source-sync.md).")
    elif not cases:
        lines.append("No fidelity cases found in the truth file.")
    else:
        lines.append(f"{failures} fidelity case(s) FAILED. This sync must not be "
                     "marked ready until the C++ core matches live R.")
    lines.append("")

    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")

    passed = len(cases) - failures - xfail
    print(f"fidelity: {passed}/{len(cases)} passed, {failures} failed, {xfail} xfail, "
          f"{xpass} xpass; wrote {report_path}")
    for name, status, _, worst, note in rows:
        if status in ("FAIL", "XFAIL"):
            print(f"  {status} {name}: max abs diff {worst:.3g} {note}")
    # No cases is treated as failure: we cannot claim fidelity without evidence.
    return 0 if (failures == 0 and cases) else 1


if __name__ == "__main__":
    raise SystemExit(main())
