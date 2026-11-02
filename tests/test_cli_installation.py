"""Installed entry point diagnostics must match the module users actually import."""
import subprocess
import sys
import unittest
from pathlib import Path
from batch_primer_design import __version__


class InstalledCliTests(unittest.TestCase):
    def test_module_reports_version_without_requiring_subcommand(self):
        run = subprocess.run([sys.executable, '-m', 'batch_primer_design.advanced_cli', '--version'], capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertEqual(run.stdout.strip(), 'primer-workflow '+__version__)

    def test_installed_entry_exposes_current_local_and_online_parameters(self):
        executable = Path(sys.executable).parent / ('primer-workflow.exe' if sys.platform == 'win32' else 'primer-workflow')
        self.assertTrue(executable.is_file(), 'install the project with python -m pip install -e . before running tests')
        for args, required in [(['design','--help'], ['--gff', '--assembly']), (['reference','--help'], ['inspect'])]:
            run = subprocess.run([str(executable), *args], capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            for option in required:
                self.assertIn(option, run.stdout)
