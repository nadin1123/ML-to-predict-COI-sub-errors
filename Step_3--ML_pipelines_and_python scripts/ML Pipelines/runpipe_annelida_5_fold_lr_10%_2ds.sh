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
python3 ${SCRIPTDIR}/create_folds_2datasets.py --tag "annelida-10%-lr-2ds" -l Label 5 annelida_data/Annelida_3MinOrder_oh_10per.csv annelida_data/Annelida_2MajOrder_oh_10per.csv

for fold in 00 01 02 03 04
do
    echo " = Fold ${fold}"

    # create output file, override if already present
    output=annelida-10%-lr-2ds-folded-${fold}/output.txt

    # calculate the projection for each fold
    python3 ${SCRIPTDIR}/calculate_data_projection_lr.py "annelida-10%-lr-2ds-folded-${fold}"

    # run LR on the fold
    python3 ${SCRIPTDIR}/evaluate_lr.py --classes "1,2" --fig "annelida-10%-lr-2ds-folded-${fold}" | tee $output

    # leave two blank lines between folds
    echo "\n"
done

