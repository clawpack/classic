"""
Regression tests.  Execute via:
    python regression_tests.py
to test, or
    python regression_tests.py True
to create new regression data for archiving.
"""


from clawpack.visclaw import data
import os, sys
import numpy as np
import subprocess

# Create plotdata object for reading in frames in tests below:
plotdata = data.ClawPlotData()
plotdata.outdir = '_output'

def test1():
    """
    Compile and run the code
    """
    job = subprocess.Popen(['make', 'clean'])
    return_code = job.wait()
    assert return_code == 0, "Problem with 'make clean'"
    job = subprocess.Popen(['make', '.output'])
    return_code = job.wait()
    assert return_code == 0, "Problem with 'make .output'"


def test2(save_new_regression_data=False):
    """
    Check sum of all q values in frame 1
    """
    save_new_regression_data = (save_new_regression_data in [True,'True'])

    f = plotdata.getframe(1)
    psum = f.state.q[0,:].sum()
    usum = f.state.q[1,:].sum()
    
    fname_data = 'regression_data_test2.txt'
    
    if save_new_regression_data:
        np.savetxt(fname_data,np.array([psum,usum]))
        print "*** Created new regression_data file ", fname_data
        

    # Read in archived data for comparison:
    regression_data = np.loadtxt(fname_data)
    psum_expected = regression_data[0]
    usum_expected = regression_data[1]

    tol = 1e-14
    assert np.allclose(psum,psum_expected,tol), \
        "frame 1: psum = %s, expected: %s"  % (psum, psum_expected)
    assert np.allclose(usum,usum_expected,tol), \
        "frame 1: usum = %s, expected: %s"  % (usum, usum_expected)
    print "Frame 1 OK"
    

def test3(save_new_regression_data=False):
    """
    Check sum of all q values in frame 2
    """
    save_new_regression_data = (save_new_regression_data in [True,'True'])
    
    f = plotdata.getframe(2)
    psum = f.state.q[0,:].sum()
    usum = f.state.q[1,:].sum()
    
    fname_data = 'regression_data_test3.txt'

    if save_new_regression_data:
        np.savetxt(fname_data,np.array([psum,usum]))
        print "*** Created new regression_data file ", fname_data
        

    # Read in archived data for comparison:
    regression_data = np.loadtxt(fname_data)
    psum_expected = regression_data[0]
    usum_expected = regression_data[1]

    tol = 1e-14
    assert np.allclose(psum,psum_expected,tol), \
        "frame 2: psum = %s, expected: %s"  % (psum, psum_expected)
    assert np.allclose(usum,usum_expected,tol), \
        "frame 2: usum = %s, expected: %s"  % (usum, usum_expected)
    print "Frame 2 OK"

    
if __name__=="__main__":
    test1()
    print sys.argv
    test2(*sys.argv[1:])
    test3(*sys.argv[1:])
    


