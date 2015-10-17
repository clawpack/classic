This directory is for regression tests.

To run all tests, at the command line type:
    nosetests
or, for more verbose output:
    nosetests -sv

Each test does a short run and computes a few numbers based on all solution
values to see if they agree with archived results.


Developers: To create new archived results for a test case, go into the
directory and type:
    python regression_tests.py True
Then 'git add' and issue a pull request if you believe the new results are 
more correct.

You can also run `make .output` as usual in a subdirectory if you want to
examine all the output, which is removed automatically if nosetests runs a
test and it passes.  (If it fails, the output is preserved in a directory
with a name like `Acoustics1DHeterogeneousTest_output` if the failing test
was in drirectory `acoustics_1d_heterogeneous`.)

