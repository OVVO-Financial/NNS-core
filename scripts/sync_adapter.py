"""Conservative R->NNS-core sync adapter.

This adapter is deliberately timid. The ONLY transformation it applies to the
C++ core is rewriting a managed provenance comment block that records which R
`src/**` file (and which R commit) the core file is derived from. It never edits
C++ logic.

For each changed R `src` file it consults `sync/r_src_map.json` and:

  * verifies the incoming R source still contains the expected function
    signatures (markers). If a signature is missing it REFUSES and routes the
    file to manual review, because a changed public signature means a real port
    is required, not a provenance bump.
  * for `auto: true` mappings (only those covered by the live Rcpp fidelity
    suite), updates the provenance block in each mapped core file.
  * for `auto: false` mappings, records that a human must port and extend the
    fidelity suite -- it does not invent a port.

Changed R files that are not in the map and not in the `ignore` glue list are
reported as unmapped and routed to manual review. The adapter writes a human
report and a machine-readable result JSON, and never edits anything outside the
provenance comment block.
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass, field
from pathlib import Path

DEFAULT_MAP = Path("sync/r_src_map.json")
DEFAULT_REPORT = Path("sync/last_sync_report.md")
DEFAULT_RESULT = Path("sync/last_sync_result.json")

COMMIT_RE = re.compile(r"^[0-9a-fA-F]{7,40}$")


@dataclass
class Outcome:
    updated: list[dict] = field(default_factory=list)       # auto provenance bumps
    manual: list[dict] = field(default_factory=list)        # mapped but needs a human
    refused: list[dict] = field(default_factory=list)       # markers/signatures missing
    unmapped: list[str] = field(default_factory=list)       # changed, not mapped, not glue
    ignored: list[str] = field(default_factory=list)        # known glue/build files
    out_of_scope: list[str] = field(default_factory=list)   # changed outside source_root
    changed_core_files: list[str] = field(default_factory=list)


def load_changed_files(raw: str | None, source_root: str) -> list[str]:
    if not raw:
        return []
    text = raw
    candidate = Path(raw)
    if candidate.is_file():
        text = candidate.read_text(encoding="utf-8")
    text = text.strip()
    if not text:
        return []
    try:
        data = json.loads(text)
    except json.JSONDecodeError as exc:
        raise SystemExit(f"invalid --changed-files-json: {exc}")
    if not isinstance(data, list):
        raise SystemExit("--changed-files-json must contain a JSON array")
    return [str(item).strip() for item in data if str(item).strip()]


def update_provenance(core_path: Path, prov: dict, r_src: str, commit: str) -> str:
    """Ensure the managed provenance block records `r_src @ commit`.

    Returns one of: "updated", "added-line", "created-block", "unchanged".
    Creates the block (after the SPDX line, else at top) when absent.
    """
    begin = prov["begin_marker"]
    end = prov["end_marker"]
    src_line = prov["source_line_format"].format(r_src=r_src, commit=commit)
    # Match any commit for this r_src so we can rewrite it.
    line_re = re.compile(
        r"^// NNS-R-SOURCE: " + re.escape(r_src) + r" @ [0-9a-fA-F]{7,40}\s*$"
    )

    text = core_path.read_text(encoding="utf-8")
    lines = text.splitlines()

    try:
        b = next(i for i, ln in enumerate(lines) if ln.strip() == begin)
        e = next(i for i, ln in enumerate(lines) if ln.strip() == end)
    except StopIteration:
        b = e = None

    if b is not None and e is not None and b < e:
        block = lines[b + 1 : e]
        for i, ln in enumerate(block):
            if line_re.match(ln.strip()):
                if ln.strip() == src_line:
                    return "unchanged"
                block[i] = src_line
                lines[b + 1 : e] = block
                core_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
                return "updated"
        # block present but no line for this r_src -> append one
        block.append(src_line)
        lines[b + 1 : e] = block
        core_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return "added-line"

    # No block: create one after the SPDX line if present, else at the very top.
    insert_at = 0
    for i, ln in enumerate(lines):
        if "SPDX-License-Identifier" in ln:
            insert_at = i + 1
            break
    new_block = ["", begin, src_line, end]
    lines[insert_at:insert_at] = new_block
    core_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return "created-block"


def process(args) -> Outcome:
    mapping_doc = json.loads(Path(args.map).read_text(encoding="utf-8"))
    source_root = mapping_doc.get("source_root", "src")
    prov = mapping_doc["provenance"]
    ignore = set(mapping_doc.get("ignore", []))
    by_src = {m["r_src"]: m for m in mapping_doc["mappings"]}

    r_checkout = Path(args.r_checkout).resolve()
    commit = args.r_commit
    if not COMMIT_RE.match(commit):
        raise SystemExit(f"--r-commit does not look like a git sha: {commit!r}")

    if args.init:
        changed = sorted(by_src.keys())
    else:
        changed = load_changed_files(args.changed_files_json, source_root)

    out = Outcome()
    core_files_touched: set[str] = set()

    for rel in changed:
        prefix = source_root.rstrip("/") + "/"
        if not (rel == source_root or rel.startswith(prefix)):
            out.out_of_scope.append(rel)
            continue
        if rel in ignore:
            out.ignored.append(rel)
            continue
        mapping = by_src.get(rel)
        if mapping is None:
            out.unmapped.append(rel)
            continue

        # Verify the incoming R source still has the expected signatures.
        r_file = r_checkout / rel
        if not r_file.is_file():
            out.refused.append({"r_src": rel, "reason": "incoming R source file is missing"})
            out.manual.append({"r_src": rel, "reason": "incoming R source file is missing"})
            continue
        r_text = r_file.read_text(encoding="utf-8", errors="replace")
        missing = [s for s in mapping.get("required_r_signatures", []) if s not in r_text]
        if missing:
            reason = "missing required R signature(s): " + ", ".join(missing)
            out.refused.append({"r_src": rel, "reason": reason})
            out.manual.append({"r_src": rel, "reason": reason})
            continue

        if not mapping.get("auto", False):
            out.manual.append(
                {
                    "r_src": rel,
                    "reason": "mapped but auto=false (not covered by the fidelity "
                    "suite); requires a human port and a fidelity case",
                }
            )
            continue

        # auto + signatures intact -> rewrite provenance only.
        actions = []
        for core_rel in mapping["core_paths"]:
            core_path = Path(core_rel)
            if not core_path.is_file():
                actions.append({"core": core_rel, "action": "missing-core-file"})
                continue
            action = update_provenance(core_path, prov, rel, commit)
            actions.append({"core": core_rel, "action": action})
            if action != "unchanged":
                core_files_touched.add(core_rel)
        out.updated.append({"r_src": rel, "commit": commit, "actions": actions})

    out.changed_core_files = sorted(core_files_touched)
    return out


def write_report(out: Outcome, args, path: Path) -> None:
    L: list[str] = []
    L.append("# R src conservative sync report")
    L.append("")
    L.append(f"- R repo: `{args.r_repo}`")
    L.append(f"- R commit: `{args.r_commit}`")
    L.append(f"- mode: {'init/reseed' if args.init else 'dispatch'}")
    L.append("")

    L.append("## Auto-updated provenance (fidelity-covered)")
    if out.updated:
        for u in out.updated:
            L.append(f"- `{u['r_src']}` -> `{u['commit']}`")
            for a in u["actions"]:
                L.append(f"  - `{a['core']}`: {a['action']}")
    else:
        L.append("- none")
    L.append("")

    L.append("## Manual review required")
    if out.manual:
        for m in out.manual:
            L.append(f"- `{m['r_src']}`: {m['reason']}")
    else:
        L.append("- none")
    L.append("")

    L.append("## Unmapped changed R src (no port invented)")
    if out.unmapped:
        for f in out.unmapped:
            L.append(f"- `{f}`")
    else:
        L.append("- none")
    L.append("")

    if out.refused:
        L.append("## Refused (marker/signature drift)")
        for r in out.refused:
            L.append(f"- `{r['r_src']}`: {r['reason']}")
        L.append("")

    if out.ignored:
        L.append("## Ignored Rcpp glue / build files")
        for f in out.ignored:
            L.append(f"- `{f}`")
        L.append("")

    if out.out_of_scope:
        L.append("## Out of scope (outside source_root)")
        for f in out.out_of_scope:
            L.append(f"- `{f}`")
        L.append("")

    needs_manual = bool(out.manual or out.unmapped or out.refused)
    L.append("## Result")
    L.append("")
    if needs_manual:
        L.append(
            "Manual attention is required. This sync must NOT be marked ready "
            "for merge until the listed files are ported by a human and the live "
            "Rcpp fidelity suite passes."
        )
    elif out.updated:
        L.append(
            "Only provenance was updated for fidelity-covered files. Readiness "
            "still depends on the live Rcpp fidelity suite passing."
        )
    else:
        L.append("No in-scope R `src/**` changes required adapter action.")
    L.append("")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(L), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r-checkout", required=True)
    parser.add_argument("--r-commit", required=True)
    parser.add_argument("--r-repo", default="OVVO-Financial/NNS")
    parser.add_argument("--changed-files-json", default=None)
    parser.add_argument("--map", default=str(DEFAULT_MAP))
    parser.add_argument("--report", default=str(DEFAULT_REPORT))
    parser.add_argument("--result-json", default=str(DEFAULT_RESULT))
    parser.add_argument(
        "--init",
        action="store_true",
        help="process every mapping (seed/refresh provenance blocks)",
    )
    parser.add_argument("--github-output", default=None)
    args = parser.parse_args()

    out = process(args)
    write_report(out, args, Path(args.report))

    needs_manual = bool(out.manual or out.unmapped or out.refused)
    ready_candidate = (not needs_manual) and bool(out.updated)
    result = {
        "r_repo": args.r_repo,
        "r_commit": args.r_commit,
        "updated": out.updated,
        "manual": out.manual,
        "refused": out.refused,
        "unmapped": out.unmapped,
        "ignored": out.ignored,
        "out_of_scope": out.out_of_scope,
        "changed_core_files": out.changed_core_files,
        "needs_manual": needs_manual,
        "ready_candidate": ready_candidate,
    }
    result_path = Path(args.result_json)
    result_path.parent.mkdir(parents=True, exist_ok=True)
    result_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    gh_out = args.github_output
    if gh_out is None:
        import os

        gh_out = os.environ.get("GITHUB_OUTPUT")
    if gh_out:
        with open(gh_out, "a", encoding="utf-8") as fh:
            fh.write(f"needs_manual={'true' if needs_manual else 'false'}\n")
            fh.write(f"ready_candidate={'true' if ready_candidate else 'false'}\n")
            fh.write(f"updated_count={len(out.updated)}\n")
            fh.write(f"unmapped_count={len(out.unmapped)}\n")
            fh.write(f"manual_count={len(out.manual)}\n")

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
