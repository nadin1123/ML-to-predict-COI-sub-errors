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
python3 ${SCRIPTDIR}/create_test_train_split_2datasets.py \
        --tag "annelida-2%-rf-2ds" -l "Label" \
        annelida_data/Annelida_3MinorOrders_oh_2per.csv \
        annelida_data/Annelida_2Orders_oh_2per.csv

# calculate the projection for this split
python3 ${SCRIPTDIR}/calculate_data_projection_rf.py \
                "annelida-2%-rf-2ds-train-vs-test"

# run RF on the data
python3 ${SCRIPTDIR}/evaluate_rf.py --fig "annelida-2%-rf-2ds-train-vs-test"