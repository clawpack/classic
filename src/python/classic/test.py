r"""
Defines the Classic Clawpack Test Runner class for running PyTest based
regression tests in classic clawpack.

Refer to the documentation for PyTest to manage output and reporting.
"""

from pathlib import Path
from typing import Optional
import os

import clawpack.clawutil.test as test

# Clean library files whenever this module is used
if "CLAW" in os.environ:
    CLAW = Path(os.environ["CLAW"])
else:
    raise ValueError("Need to set CLAW environment variable.")

class ClassicTestRunner(test.ClawpackTestRunner):

    def __init__(self, path: Path, test_path: Optional[Path]=None):
        super(ClassicTestRunner, self).__init__(path, test_path=test_path)
