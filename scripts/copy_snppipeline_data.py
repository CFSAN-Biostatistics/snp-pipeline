#!/usr/bin/env python

import subprocess
import sys

subprocess_args = ["cfsan_snp_pipeline", "data"] + sys.argv[1:]
ret = subprocess.call(subprocess_args)
exit(ret)
