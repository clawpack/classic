#!/usr/bin/env python

from pathlib import Path
import pytest

import clawpack.classic.test as test


def test_acoustics_1d_example1(tmp_path: Path, save: bool):
    runner = test.ClassicTestRunner(tmp_path,
                                    test_path=Path(__file__).parent)

    runner.set_data()

    runner.rundata.clawdata.num_output_times = 2
    runner.rundata.clawdata.tfinal = 1.0
    runner.rundata.clawdata.output_t0 = False

    runner.write_data()

    runner.executable_name = "xclaw"
    runner.build_executable()
    runner.run_code()

    runner.check_frame(1, indices=(0, 1), save=save)
    runner.check_frame(2, indices=(0, 1), save=save)

if __name__=="__main__":
    raise SystemExit(pytest.main([__file__]))
