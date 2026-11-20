"""Smoluchowski equation with fragmentation (Ref. 1). Kept for backward compatibility.

Equivalent to ``python -m smolu_dsmc frag [options]``; run with ``--help`` for the options.
"""

import sys

from smolu_dsmc.cli import main

if __name__ == "__main__":
    main(["frag", *sys.argv[1:]])
