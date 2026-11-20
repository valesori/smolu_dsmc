"""Smoluchowski equation with fragmentation and entanglement (Ref. 1).

Kept for backward compatibility. Equivalent to
``python -m smolu_dsmc frag --gamma 7e-6 [options]``; run with ``--help`` for the options.
"""

import sys

from smolu_dsmc.cli import main

if __name__ == "__main__":
    main(["frag", "--gamma", "7e-6", *sys.argv[1:]])
