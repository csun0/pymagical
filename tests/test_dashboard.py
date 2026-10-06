"""Small dashboard checks; runnable without inference dependencies."""

from contextlib import redirect_stderr, redirect_stdout
import io
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from pymagical.cli import main
from pymagical.dashboard import generate_dashboard, RESULT_COLUMNS


HEADER = "\t".join(RESULT_COLUMNS) + "\tTFs(prob)\n"
ROW = "GeneA\tchr1\t100\tchr1\t80\t120\t0.9\tTF1 (0.8),\n"


class DashboardTests(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)

    def result(self, name="sample", row=ROW):
        (self.directory / f"{name}.txt").write_text(HEADER + row, encoding="utf-8")

    def invoke(self, *args):
        with patch.object(sys, "argv", ["pymagical", "dashboard", *args]):
            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                main()

    def test_embeds_multiple_results_and_skips_sidecars(self):
        self.result("astrocytes")
        self.result("neurons")
        (self.directory / "sample_B_matrix.txt").write_text("1\t2\n")
        (self.directory / "sample_timing_stats.txt").write_text("time: 10\n")
        (self.directory / "nested").mkdir()
        (self.directory / "nested" / "unused.txt").write_text(HEADER + ROW)
        path = generate_dashboard(self.directory)
        html = path.read_text()
        self.assertEqual(path, self.directory / "dashboard.html")
        self.assertNotIn('src="./data.js"', html)
        payload = html.split("window.MAGICAL_RESULTS = ", 1)[1].split(";</script>", 1)[0]
        self.assertEqual(set(json.loads(payload)), {"astrocytes", "neurons"})

    def test_embedded_data_cannot_close_script_and_supports_bom(self):
        row = ROW.replace("GeneA", "</script><script>alert(1)</script>é")
        (self.directory / "sample.txt").write_text(HEADER + row, encoding="utf-8-sig")
        output = self.directory / "custom.html"
        self.assertEqual(generate_dashboard(self.directory, output), output)
        html = output.read_text()
        self.assertNotIn("</script><script>alert(1)", html)
        payload = html.split("window.MAGICAL_RESULTS = ", 1)[1].split(";</script>", 1)[0]
        self.assertEqual(json.loads(payload)["sample"], HEADER + row)

    def test_missing_and_empty_directories(self):
        with self.assertRaisesRegex(ValueError, "not a directory"):
            generate_dashboard(self.directory / "missing")
        with self.assertRaisesRegex(ValueError, "No MAGICAL result"):
            generate_dashboard(self.directory)

    def test_required_flag_and_cli_errors(self):
        for args in [(), ("--input_dir", str(self.directory)),
                     ("--input_dir", str(self.directory / "missing"))]:
            with self.subTest(args=args), self.assertRaises(SystemExit) as exc:
                self.invoke(*args)
            self.assertEqual(exc.exception.code, 2)

    def test_cli_browser_and_no_open_alias(self):
        self.result()
        with patch("pymagical.dashboard.webbrowser.open", return_value=True) as browser:
            self.invoke("--input_dir", str(self.directory))
            browser.assert_called_once_with((self.directory / "dashboard.html").as_uri())
            browser.reset_mock()
            self.invoke("--input-dir", str(self.directory), "--no-open")
            browser.assert_not_called()


if __name__ == "__main__":
    unittest.main()
