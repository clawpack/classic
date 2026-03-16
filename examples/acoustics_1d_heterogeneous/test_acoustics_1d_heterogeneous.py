#!/usr/bin/env python
"""
Regression tests for a 1D heterogeneous acoustics test
"""

from pathlib import Path
import pytest
import clawpack.classic.test as test

def test_acoustics_1d_heterogeneous(tmp_path: Path, save: bool):
    # Create an example-local runner.  tmp_path keeps all generated files and
    # solver output out of the source tree.
    runner = test.ClassicTestRunner(tmp_path, test_path=Path(__file__).parent)

    # Load the default example configuration from setrun.py and make the run
    # smaller/faster for regression testing.
    runner.set_data()
    runner.rundata.clawdata.num_output_times = 2
    runner.rundata.clawdata.tfinal = 5.0
    runner.rundata.clawdata.output_t0 = False
    
    # Write the input data files into the temporary run directory.
    runner.write_data()

    # Build using the example Makefile, then move the executable into tmp_path.
    runner.build_executable()

    # Run the code with both input files and output files in tmp_path.
    runner.run_code()

    # Compare selected output frames against saved regression baselines.
    runner.check_frame(1, indices=(0, 1), save=save)
    runner.check_frame(2, indices=(0, 1), save=save)

if __name__=="__main__":
    pytest.main([__file__])
