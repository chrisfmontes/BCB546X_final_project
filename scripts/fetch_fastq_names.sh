#!/bin/bash
set -e #prevents the shell script from proceeding if one of its commands fails.
set -u #aborts a script if an unset variable is encountered
set -o pipefail #after "set -e" is provided, this indicates any error in a pipe should cause exit.

esearch -db sra -query $1 | efetch --format runinfo |cut -d "," -f 1 > $2
# fetch from sra (-db sra) all the information for the Bioproject accession indicated as $1, then get the run information
# and finally extract the column that has the SRR accessions. Output that information to the file $2.