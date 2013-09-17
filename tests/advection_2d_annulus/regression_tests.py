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
    my2 = int(f.state.q.shape[-1] / 2)
    qsum1 = f.state.q[0,:,:my2].sum()
    qsum2 = f.state.q[0,:,my2:].sum()
    
    fname_data = 'regression_data_test2.txt'
    
    if save_new_regression_data:
        np.savetxt(fname_data,np.array([qsum1,qsum2]))
        print "*** Created new regression_data file ", fname_data
        

    # Read in archived data for comparison:
    regression_data = np.loadtxt(fname_data)
    qsum1_expected = regression_data[0]
    qsum2_expected = regression_data[1]

    tol = 1e-14
    assert np.allclose(qsum1,qsum1_expected,tol), \
        "frame 1: qsum1 = %s, expected: %s"  % (qsum1, qsum1_expected)
    assert np.allclose(qsum2,qsum2_expected,tol), \
        "frame 1: qsum2 = %s, expected: %s"  % (qsum2, qsum2_expected)
    print "Frame 1 OK"
    

def test3(save_new_regression_data=False):
    """
    Check sum of all q values in frame 2
    """
    save_new_regression_data = (save_new_regression_data in [True,'True'])
    
    f = plotdata.getframe(2)
    my2 = int(f.state.q.shape[-1] / 2)
    qsum1 = f.state.q[0,:,:my2].sum()
    qsum2 = f.state.q[0,:,my2:].sum()
    
    fname_data = 'regression_data_test3.txt'

    if save_new_regression_data:
        np.savetxt(fname_data,np.array([qsum1,qsum2]))
        print "*** Created new regression_data file ", fname_data
        

    # Read in archived data for comparison:
    regression_data = np.loadtxt(fname_data)
    qsum1_expected = regression_data[0]
    qsum2_expected = regression_data[1]

    tol = 1e-14
    assert np.allclose(qsum1,qsum1_expected,tol), \
        "frame 2: qsum1 = %s, expected: %s"  % (qsum1, qsum1_expected)
    assert np.allclose(qsum2,qsum2_expected,tol), \
        "frame 2: qsum2 = %s, expected: %s"  % (qsum2, qsum2_expected)
    print "Frame 2 OK"

    
if __name__=="__main__":
    test1()
    print sys.argv
    test2(*sys.argv[1:])
    test3(*sys.argv[1:])
    


