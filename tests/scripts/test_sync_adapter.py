import importlib.util
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ADAPTER = ROOT / "scripts" / "sync_adapter.py"
_spec = importlib.util.spec_from_file_location("sync_adapter", ADAPTER)
sync_adapter = importlib.util.module_from_spec(_spec)
sys.modules["sync_adapter"] = sync_adapter  # dataclasses look the module up
_spec.loader.exec_module(sync_adapter)
code_signature = sync_adapter.code_signature


class CodeSignatureTest(unittest.TestCase):
    def test_comments_and_whitespace_are_ignored(self):
        a = "int f(int x) {\n  return x + 1;  // add one\n}\n"
        b = "/* doc */ int f(int x)\n{\n\treturn x + 1;\n}"
        self.assertEqual(code_signature(a), code_signature(b))

    def test_code_changes_are_seen(self):
        self.assertNotEqual(code_signature("return x + 1;"), code_signature("return x + 2;"))

    def test_literals_are_kept_byte_for_byte(self):
        self.assertNotEqual(code_signature('sizeof("a b")'), code_signature('sizeof("a  b")'))
        self.assertNotEqual(code_signature("c == ' '"), code_signature("c == '\t'"))
        self.assertEqual(code_signature('f("a // b");'), 'f("a // b");')


class LiteralChangeRoutingTest(unittest.TestCase):
    def test_changed_literal_needs_manual_review(self):
        decl = '\nint literal_size() { return sizeof("a b"); }\n'
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            for rel in ("src/partial_moments.cpp", "include/nns/partial_moments.hpp",
                        "original_src/partial_moments.cpp", "original_src/partial_moments.h"):
                (work / rel).parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(ROOT / rel, work / rel)
            incoming = work / "incoming" / "src"
            incoming.mkdir(parents=True)
            for name in ("partial_moments.cpp", "partial_moments.h"):
                shutil.copyfile(ROOT / "original_src" / name, incoming / name)
            snapshot = work / "original_src" / "partial_moments.cpp"
            snapshot.write_text(snapshot.read_text(encoding="utf-8") + decl, encoding="utf-8")
            cpp = incoming / "partial_moments.cpp"
            cpp.write_text(cpp.read_text(encoding="utf-8") + decl.replace("a b", "a  b"),
                           encoding="utf-8")
            proc = subprocess.run(
                [sys.executable, str(ADAPTER), "--r-checkout", str(work / "incoming"),
                 "--r-commit", "1234567", "--map", str(ROOT / "sync" / "r_src_map.json"),
                 "--changed-files-json", '["src/partial_moments.cpp", "src/partial_moments.h"]'],
                cwd=work, capture_output=True, text=True, check=True)
            result = json.loads(proc.stdout)
            self.assertTrue(result["needs_manual"])
            self.assertFalse(result["ready_candidate"])


if __name__ == "__main__":
    unittest.main()
