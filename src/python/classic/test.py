r"""
Execute nosetests in all subdirectories, to run a series of quick
regression tests.

Sends output and result/errors to separate files to simplify checking
results and looking for errors.
"""

from pathlib import Path
import os
import sys
import subprocess
import importlib
import inspect
import shutil
import pytest

import numpy as np

import clawpack.clawutil.runclaw as runclaw
import clawpack.clawutil.claw_git_status as claw_git_status
import clawpack.pyclaw.solution as solution

# Clean library files whenever this module is used
if "CLAW" in os.environ:
    CLAW = Path(os.environ["CLAW"])
else:
    raise ValueError("Need to set CLAW environment variable.")

for lib_path in (CLAW / "classic" / "src" / "1d").glob("*.o"):
    lib_path.unlink()
for lib_path in (CLAW / "classic" / "src" / "2d").glob("*.o"):
    lib_path.unlink()
for lib_path in (CLAW / "classic" / "src" / "2d").glob("*.o"):
    lib_path.unlink()


class ClawpackClassicTestRunner:

    def __init__(self, path, caller_path):

        self.temp_path = path
        # :TODO: see if there's a way to get this automatically
        self.test_path = Path(caller_path).parent
        self.executable_name = 'xclaw'


    def set_data(self, setrun_module=None):

        sys.path.insert(0, self.test_path)
        if not setrun_module:
            setrun_module = 'setrun'
        if setrun_module in sys.modules:
            del(sys.modules[setrun_module])
        setrun = importlib.import_module(setrun_module)
        self.rundata = setrun.setrun()
        sys.path.pop(0)


    def write_data(self, path=None):

        if not path:
            path = self.temp_path
        self.rundata.write(out_dir=path)


    def build_executable(self, make_level='default', FFLAGS=None, LFLAGS=None):

        # Assumes GCC CLI
        if not FFLAGS:
            FFLAGS = os.environ.get('FFLAGS', "-O2 -fopenmp")
        if not LFLAGS:
            LFLAGS = os.environ.get('LFLAGS', FFLAGS)

        if make_level.lower() == "new":
            cmd = "".join((f"cd {self.test_path} ; make new ",
                           f"FFLAGS='{FFLAGS}' LFLAGS='{LFLAGS}'"))
        elif make_level.lower() == "default":
            # clean up *.o and *.mod files in test path only
            for path in self.test_path.glob("*.o"):
                path.unlink()
            for path in self.test_path.glob("*.mod"):
                path.unlink()
            cmd = "".join((f"cd {self.test_path} ; make .exe ",
                           f"FFLAGS='{FFLAGS}' LFLAGS='{LFLAGS}'"))

        elif make_level.lower() == "exe":
            cmd = "".join((f"cd {self.test_path} ; make .exe ",
                           f"FFLAGS='{FFLAGS}' LFLAGS='{LFLAGS}'"))
        else:
            raise ValueError(f"Invaled make_level={make_level} given.")

        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            self.clean_up()
            raise e

        shutil.move(self.test_path / self.executable_name, self.temp_path)


    def run_code(self):
        runclaw.runclaw(xclawcmd=self.temp_path / self.executable_name,
                        rundir=self.temp_path,
                        outdir=self.temp_path,
                        overwrite=True,
                        restart=False)


    def clean_up(self):
        pass


    def check_frame(self, frame, indices=(0), regression_path=None, save=False, **kwargs):

        if not(isinstance(indices, tuple) or isinstance(indices, list)):
            indices = tuple(indices)

        if not regression_path:
            regression_path = self.test_path / "regression_data"

        # Load test output data
        sol = solution.Solution(frame, path=self.temp_path)
        sol_sums = sol.q[indices, ...].sum(axis=1)

        # Load regression data
        regression_data = regression_path / f"frame{str(frame).zfill(4)}.txt"
        if save:
            np.savetxt(regression_data, sol_sums)
            claw_git_status.make_git_status_file(outdir=regression_path)
        regression_sum = np.loadtxt(regression_data)

        # Compare data
        kwargs.setdefault('rtol', 1e-14)
        kwargs.setdefault('atol', 1e-8)
        np.testing.assert_allclose(sol_sums, regression_sum, **kwargs)


    def check_gauge(self, gauge_id, indices=(0), regression_path=None, save=False, **kwargs):
        r"""Basic test to assert gauge equality

        :Input:
         - *save* (bool) - If *True* will save the output from this test to 
           the file *regresion_data.txt*.  Default is *False*.
         - *indices* (tuple) - Contains indices to compare in the gague 
           comparison.  Defaults to *(0)*.
         - *rtol* (float) - Relative tolerance used in the comparison, default 
           is *1e-14*.  Note that the old *tolerance* input is now synonymous 
           with this parameter.
         - *atol* (float) - Absolute tolerance used in the comparison, default
           is *1e-08*.
        """

        if not(isinstance(indices, tuple) or isinstance(indices, list)):
            indices = tuple(indices)

        if not regression_path:
            regression_path = self.test_path / "regression_data"

        # Load test output data
        gauge = gauges.GaugeSolution(gauge_id, path=self.temp_path)

        # Load regression data
        if save:
            shutil.copy(self.temp_path / f"gauge{str(gauge_id).zfill(5)}.txt",
                        regression_path)
            claw_git_status.make_git_status_file(outdir=regression_path)
        regression_gauge = gauges.GaugeSolution(gauge_id, path=regression_path)

        # Compare data
        kwargs.setdefault('rtol', 1e-14)
        kwargs.setdefault('atol', 1e-8)
        try:
            for n in indices:
                np.testing.assert_allclose(gauge.q[n, :],
                                              regression_gauge.q[n, :], 
                                              **kwargs)
        except AssertionError as e:
            err_msg = "\n".join((e.args[0], 
                                "Gauge Match Failed for gauge = %s" % gauge_id))
            err_msg = "\n".join((err_msg, "  failures in fields:"))
            failure_indices = []
            for n in indices:
                if ~np.allclose(gauge.q[n, :], regression_gauge.q[n, :], 
                                **kwargs):
                    failure_indices.append(str(n))
            index_str = ", ".join(failure_indices)
            raise AssertionError(" ".join((err_msg, index_str)))