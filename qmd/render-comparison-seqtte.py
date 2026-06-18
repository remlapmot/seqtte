#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "jupyterlab>=4.4.3",
#     "jupyterlab-stata-highlight2>=0.1.2",
#     "nbstata>=0.8.5",
# ]
# ///
import subprocess

cmd0 = "python -m nbstata.install --sys-prefix"
retval0 = subprocess.call(cmd0, shell=True)
print('returned value:', retval0)

cmd1 = "quarto render comparison-seqtte.qmd"
retval1 = subprocess.call(cmd1, shell=True)
print('returned value:', retval1)
