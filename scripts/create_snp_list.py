#!/usr/bin/env python

import subprocess
import sys

subprocess_args = ["cfsan_snp_pipeline", "merge_sites"] + sys.argv[1:]
ret = subprocess.call(subprocess_args)
exit(ret)
