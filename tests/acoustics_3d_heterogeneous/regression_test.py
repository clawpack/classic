"""
Regression tests for 3D heterogeneous acoustics problem.
"""

from __future__ import absolute_import
import sys
import unittest

import clawpack.classic.test as test


class Acoustics3DHeterogeneousTest(test.ClassicRegressionTest):
    r"""Basic test for a 3D heterogeneous acoustics test."""


    def runTest(self, save=False):

        # Write out data files
        self.load_rundata()
        self.write_rundata_objects()

        # Run code
        self.run_code()

        # Perform tests
        self.check_frame(save=save, indices=[0, 1, 2], frame_num=1,
                         file_name='regression_data_test2.txt')
        self.check_frame(save=save, indices=[0, 1, 2], frame_num=2,
                         file_name='regression_data_test3.txt')

        self.success = True



if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = Acoustics3DHeterogeneousTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()