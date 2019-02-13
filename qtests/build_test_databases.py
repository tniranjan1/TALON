# This script builds the databases needed for the TALON testing suite to run

import os
import subprocess
import sys

try:
   subprocess.check_output(
       ["python", "../initialize_talon_database.py",
        "--f", "input_files/toy_transcript/toy_annot.gtf",
        "--a",  "toy_annot",
        "--g",  "toy_build", "--o", "scratch/toy"])
except Exception as e:
    print(e)
    sys.exit("Database initialization failed on toy artificial transcript")