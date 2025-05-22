"""
Regression tests for 3D heterogeneous acoustics problem.
"""

import sys
import unittest

import clawpack.classic.test as test


def test_acoustics_3d_heterogeneous(tmp_path, save):
    r"""Basic test for a 3D heterogeneous acoustics test."""

    ctr = test.ClawpackClassicTestRunner(tmp_path, __file__)

    ctr.set_data()
    ctr.rundata.clawdata.num_cells = [20, 20, 20]
    ctr.rundata.clawdata.num_output_times = 2
    ctr.rundata.clawdata.tfinal = 1.0
    ctr.write_data()

    ctr.executable_name = 'xclaw'
    ctr.build_executable()

    ctr.run_code()

    ctr.check_frame(1, indices=(0, 1, 2), save=save)
    ctr.check_frame(2, indices=(0, 1, 2), save=save)

if __name__=="__main__":
    pytest.main([__file__])