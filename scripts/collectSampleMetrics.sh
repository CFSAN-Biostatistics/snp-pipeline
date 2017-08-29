#!/bin/bash
#

>&2 echo ""
>&2 echo "Warning: this script is deprecated and will be removed in a future SNP Pipeline release."
>&2 echo "Use cfsan_snp_pipeline."
>&2 echo ""

cfsan_snp_pipeline collect_metrics "$@"
