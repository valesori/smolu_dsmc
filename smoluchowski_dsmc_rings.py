"""Irreversible Smoluchowski equation with ring formation (Ref. 2). Kept for backward compatibility.

Equivalent to ``python -m smolu_dsmc rings [options]``; run with ``--help`` for the options.
"""

import sys

from smolu_dsmc.cli import main

if __name__ == "__main__":
    main(["rings", *sys.argv[1:]])
