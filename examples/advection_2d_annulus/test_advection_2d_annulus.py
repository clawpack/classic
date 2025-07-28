"""
Regression tests for 2D advection on an annulus.
"""

from pathlib import Path

import clawpack.classic.test as test

def test_advection_2d_annulus(tmp_path: Path, save: bool):

    ctr = test.ClawpackClassicTestRunner(tmp_path)

    ctr.set_data()
    ctr.rundata.clawdata.num_output_times = 2
    ctr.rundata.clawdata.tfinal = 0.5
    ctr.write_data()

    ctr.executable_name = 'xclaw'
    ctr.build_executable()

    ctr.run_code()

    ctr.check_frame(1, save=save)
    ctr.check_frame(2, save=save)

if __name__=="__main__":
    pytest.main([__file__])
