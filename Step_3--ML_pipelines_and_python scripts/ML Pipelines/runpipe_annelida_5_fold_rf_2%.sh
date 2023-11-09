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

# create the "folds"
python3 ${SCRIPTDIR}/create_folds.py --tag "annelida-2%-rf" -l Label 5 annelida_data/Annelida_2MajOrder_oh_2per.csv

for fold in 00 01 02 03 04
do
    echo " = Fold ${fold}"

    # create output file, override if already present
    output=annelida-2%-rf-folded-${fold}/output.txt

    # calculate the projection for each fold
    python3 ${SCRIPTDIR}/calculate_data_projection_rf.py "annelida-2%-rf-folded-${fold}"

    # run RF on the fold
    python3 ${SCRIPTDIR}/evaluate_rf.py --classes "1,2" --fig "annelida-2%-rf-folded-${fold}" | tee $output

    # leave two blank lines between folds
    echo "\n"
done

