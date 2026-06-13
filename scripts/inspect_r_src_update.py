from __future__ import annotations

import argparse
import json
from pathlib import Path

from hash_r_src import hash_tree

MANIFEST_PATH = Path("sync/nns_r_source.json")
INSPECTION_PATH = Path("sync/last_r_src_inspection.md")


def load_changed_files(raw: str | None) -> list[str]:
    """Parse the --changed-files-json argument.

    Accepts either a path to a JSON file or an inline JSON string. Returns an
    empty list when nothing usable is provided.
    """
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

    return [str(item) for item in data]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r-checkout", type=Path, required=True)
    parser.add_argument("--r-commit", required=True)
    parser.add_argument("--r-version", required=True)
    parser.add_argument("--changed-files-json", default=None)
    args = parser.parse_args()

    r_checkout = args.r_checkout.resolve()
    if not (r_checkout / "src").is_dir():
        raise SystemExit(f"missing src directory in R checkout: {r_checkout / 'src'}")

    if not MANIFEST_PATH.is_file():
        raise SystemExit(f"missing manifest: {MANIFEST_PATH}")

    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    previous_hash = manifest.get("r_src_tree_hash", "")
    r_repo = manifest.get("r_repo", "OVVO-Financial/NNS")

    incoming_hash = hash_tree(r_checkout)
    changed_files = load_changed_files(args.changed_files_json)
    changed = incoming_hash != previous_hash

    lines: list[str] = []
    lines.append("# R src sync inspection")
    lines.append("")
    lines.append(f"- R repo: {r_repo}")
    lines.append(f"- R commit: `{args.r_commit}`")
    lines.append(f"- R version: `{args.r_version}`")
    lines.append(f"- Previous src hash: `{previous_hash}`")
    lines.append(f"- Incoming src hash: `{incoming_hash}`")
    lines.append("- Changed files:")
    if changed_files:
        for f in changed_files:
            lines.append(f"  - `{f}`")
    else:
        lines.append("  - (not reported by dispatch payload)")
    lines.append("")
    lines.append("## Result")
    lines.append("")
    if changed:
        lines.append(
            "`src/**` changed. Review and port the affected R C++ files into "
            "NNS-core."
        )
        lines.append("")
        lines.append(
            "This inspection does not auto-port arbitrary C++ changes. After the "
            "core source has actually been updated and C++ build and tests pass, "
            "update `sync/nns_r_source.json` with the new `r_commit`, "
            "`r_version`, `r_src_tree_hash`, and `core_commit`."
        )
    else:
        lines.append(
            "`src/**` is unchanged relative to the recorded manifest hash. No "
            "NNS-core source update is required."
        )
    lines.append("")

    INSPECTION_PATH.parent.mkdir(parents=True, exist_ok=True)
    INSPECTION_PATH.write_text("\n".join(lines), encoding="utf-8")

    print(f"changed={'true' if changed else 'false'}")
    print(f"previous_hash={previous_hash}")
    print(f"incoming_hash={incoming_hash}")
    print(f"wrote {INSPECTION_PATH}")


if __name__ == "__main__":
    main()
