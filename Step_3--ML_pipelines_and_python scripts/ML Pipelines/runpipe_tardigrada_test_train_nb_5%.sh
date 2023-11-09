#!/bin/bash

#
# Runs a simple analysis pipeline on the data in data/fulldata
#

### Safety standards
# -u : tell bash to fail if any variable is not set
# -e : tell bash to fail if any command fails (unless in an if)
# -o pipefail : tell bash to fail if part of a pipe fails (needs -e)
set -e
set -u
set -o pipefail

# save the directory that our script it in so that we can find
# the tools
SCRIPTDIR=`dirname $0`

# -x Turn tracing on to make clearer what is happening
set -x

# split the data
python3 ${SCRIPTDIR}/create_test_train_split.py \
        --tag "tardigrada-5%-nb" -l "Label" \
        tardigrada_data/Tardigrada_2MajOrder_oh_5per.csv

# calculate the projection for this split
python3 ${SCRIPTDIR}/calculate_data_projection_nb.py \
                "tardigrada-5%-nb-train-vs-test"

# run SVN on the datea
python3 ${SCRIPTDIR}/evaluate_nb.py --fig "tardigrada-5%-nb-train-vs-test"

