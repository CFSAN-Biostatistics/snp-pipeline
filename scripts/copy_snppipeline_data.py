#!/usr/bin/env python

import sys
from snppipeline import cfsan_snp_pipeline

argv = list(sys.argv)
sys.argv[0] = "cfsan_snp_pipeline"  # this text appears in the usage help
argv[0] = "data"
exit(cfsan_snp_pipeline.run_command(argv))
