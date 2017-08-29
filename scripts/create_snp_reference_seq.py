#!/usr/bin/env python

from __future__ import print_function

import subprocess
import sys

print("\nWarning: this script is deprecated and will be removed in a future SNP Pipeline release.", file=sys.stderr)
print("Use cfsan_snp_pipeline.\n", file=sys.stderr)

subprocess_args = ["cfsan_snp_pipeline", "snp_reference"] + sys.argv[1:]
ret = subprocess.call(subprocess_args)
exit(ret)
