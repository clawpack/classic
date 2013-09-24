This directory is for regression tests.

To run all tests:
    python run_tests.py
or:
    make tests

To clean up afterwards, removing all executables, output, and test results:
    make clobber

Each test does a short run and computes a few numbers based on all solution
values to see if they agree with archived results.


Developers: To create new archived results for a test case, go into the
directory and type:
    python regression_tests.py True
Then 'git add' and issue a pull request if you believe the new results are 
more correct.

