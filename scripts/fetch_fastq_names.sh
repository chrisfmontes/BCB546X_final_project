#!/bin/bash
set -e
set -u
set -o pipefail

esearch -db sra -query $1 | efetch --format runinfo |cut -d "," -f 1 > $2