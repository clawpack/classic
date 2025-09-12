r"""
Execute nosetests in all subdirectories, to run a series of quick
regression tests.

Sends output and result/errors to separate files to simplify checking
results and looking for errors.
"""

from pathlib import Path
import os

import clawpack.clawutil.test as test

# Clean library files whenever this module is used
if "CLAW" in os.environ:
    CLAW = Path(os.environ["CLAW"])
else:
    raise ValueError("Need to set CLAW environment variable.")

for lib_path in (CLAW / "classic" / "src" / "1d").glob("*.o"):
    lib_path.unlink()
for lib_path in (CLAW / "classic" / "src" / "2d").glob("*.o"):
    lib_path.unlink()
for lib_path in (CLAW / "classic" / "src" / "3d").glob("*.o"):
    lib_path.unlink()

class ClawpackClassicTestRunner(test.ClawpackTestRunner):

    def __init__(self, path: Path):
        super(ClawpackClassicTestRunner, self).__init__(path)
