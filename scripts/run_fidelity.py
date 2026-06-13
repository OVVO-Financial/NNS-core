"""Compare NNS-core C++ outputs against live R (Rcpp) truth.

Reads a truth file produced by tests/fidelity/generate_truth.R (which calls the
installed R NNS package), drives the portable C++ core via the `nns_eval` CLI on
the same inputs, and asserts every value matches within tolerance.

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
META_KEYS = {"name", "func", "expected", "tolerance"}


def case_to_argv(eval_bin: str, case: dict) -> list[str]:
    argv = [eval_bin, case["func"]]
    for key, value in case.items():
        if key in META_KEYS:
            continue
        if isinstance(value, list):
            argv += [f"--{key}", ",".join(repr_num(v) for v in value)]
        else:
            argv += [f"--{key}", repr_num(value)]
    return argv


def repr_num(v) -> str:
    # Full-precision, plain decimal (avoid surprises from str() on ints/floats).
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, int):
        return str(v)
    return repr(float(v))


def within_tol(got: float, want: float, tol: float) -> bool:
    if math.isnan(got) and math.isnan(want):
        return True
    diff = abs(got - want)
    return diff <= tol or diff <= tol * max(abs(got), abs(want))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth", required=True)
    parser.add_argument("--eval-bin", required=True)
    parser.add_argument("--report", default="tests/fidelity/fidelity_report.md")
    parser.add_argument("--tol", type=float, default=None,
                        help="override tolerance from the truth file")
    args = parser.parse_args()

    truth_path = Path(args.truth)
    if not truth_path.is_file():
        raise SystemExit(f"missing truth file: {truth_path}")
    if not Path(args.eval_bin).is_file():
        raise SystemExit(f"missing eval binary: {args.eval_bin}")

    truth = json.loads(truth_path.read_text(encoding="utf-8"))
    tol = args.tol if args.tol is not None else float(truth.get("tolerance", 1e-9))
    cases = truth.get("cases", [])

    rows: list[tuple[str, bool, float, float, float]] = []
    failures = 0
    for case in cases:
        argv = case_to_argv(args.eval_bin, case)
        proc = subprocess.run(argv, capture_output=True, text=True)
        if proc.returncode != 0:
            rows.append((case["name"], False, float("nan"),
                         float(case["expected"]), float("nan")))
            failures += 1
            continue
        got = float(proc.stdout.strip())
        want = float(case["expected"])
        ok = within_tol(got, want, tol)
        rows.append((case["name"], ok, got, want, abs(got - want)))
        if not ok:
            failures += 1

    lines: list[str] = []
    lines.append("# Live Rcpp fidelity report")
    lines.append("")
    lines.append(f"- R version: `{truth.get('r_version', '?')}`")
    lines.append(f"- R commit: `{truth.get('r_commit', '?')}`")
    lines.append(f"- tolerance: `{tol:g}`")
    lines.append(f"- cases: {len(cases)}  |  failures: {failures}")
    lines.append("")
    lines.append("| case | status | C++ | R | abs diff |")
    lines.append("| --- | --- | --- | --- | --- |")
    for name, ok, got, want, diff in rows:
        status = "ok" if ok else "FAIL"
        lines.append(f"| {name} | {status} | {got:.12g} | {want:.12g} | {diff:.3g} |")
    lines.append("")
    lines.append("## Result")
    lines.append("")
    if failures == 0 and cases:
        lines.append("All fidelity cases passed. Core matches live R within tolerance.")
    elif not cases:
        lines.append("No fidelity cases found in the truth file.")
    else:
        lines.append(f"{failures} fidelity case(s) FAILED. This sync must not be "
                     "marked ready until the C++ core matches live R.")
    lines.append("")

    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")

    print(f"fidelity: {len(cases) - failures}/{len(cases)} passed; wrote {report_path}")
    # No cases is treated as failure: we cannot claim fidelity without evidence.
    return 0 if (failures == 0 and cases) else 1


if __name__ == "__main__":
    raise SystemExit(main())
