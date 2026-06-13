from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


def hash_tree(root: Path) -> str:
    src = root / "src"
    if not src.is_dir():
        raise SystemExit(f"missing src directory: {src}")

    digest = hashlib.sha256()
    for path in sorted(p for p in src.rglob("*") if p.is_file()):
        rel = path.relative_to(root).as_posix()
        digest.update(rel.encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("r_checkout", type=Path)
    args = parser.parse_args()
    print(hash_tree(args.r_checkout.resolve()))


if __name__ == "__main__":
    main()
