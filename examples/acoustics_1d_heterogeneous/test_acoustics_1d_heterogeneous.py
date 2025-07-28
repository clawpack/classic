"""
Regression tests for a 1D heterogeneous acoustics test
"""

from pathlib import Path
import pytest
import clawpack.classic.test as test

def test_acoustics_1d_heterogeneous(tmp_path: Path, save: bool):

    ctr = test.ClawpackClassicTestRunner(tmp_path)

    ctr.set_data()
    ctr.rundata.clawdata.num_output_times = 2
    ctr.rundata.clawdata.tfinal = 5.0
    ctr.rundata.clawdata.output_t0 = False
    ctr.write_data()

    ctr.executable_name = 'xclaw'
    ctr.build_executable()

    ctr.run_code()

    ctr.check_frame(1, indices=(0, 1), save=save)
    ctr.check_frame(2, indices=(0, 1), save=save)

if __name__=="__main__":
    pytest.main([__file__])
