#!/usr/bin/env python3

'''
Break the data up into training/testing
'''

import argparse
import os
import sys

import pandas as pd
import numpy as np

import random

from sklearn.model_selection import StratifiedKFold, StratifiedShuffleSplit


## get values from command line if given
argparser = argparse.ArgumentParser(
		description="Set up N fold cross validation files.",
		epilog="If a label column is not given, it will be assumed "
			"that the last column is the label column.")
argparser.add_argument("-t", "--tag", action="store",
		default=None, type=str,
		help="Give an additional tag used in output directory name "
			"creation.")
argparser.add_argument("-v", "--vsplit", action="store",
		default=0.5, type=float,
		help="Set the fraction of a 'training' block to be set "
			"aside for validation, where the rest is used for "
			"training.  There is always 1/N of the data used for "
			"testing in each fold -- this describes what happens "
			"to the (N-1)/N portion of data not used for testing.")
argparser.add_argument("-s", "--seed", action="store",
        default=None, type=int,
        help="Set the random seed -- if set to any value, every "
            "run of the program will use the same sequence of "
            "random values")
argparser.add_argument("-l", "--label", action="store",
        default=None, type=str,
        help="Indicate the label column -- if not given, the last "
			"column is assumed to be the label")
argparser.add_argument("N", default=10, type=int,
		help="This determines how many folds to divide the data into. "
			"An independent sample of 1/N of the data will be used for "
			"testing in each fold, and the remaining data will be used "
			"for training and validation.  See also '--vsplit'.")
argparser.add_argument("filename",
		help="The filename to use for input, which is assumed to contain "
			"all available data.")

args = argparser.parse_args(sys.argv[1:])



## set the seed if given for reproducible runs
if args.seed is not None:
    random.seed(argparse.seed)

if args.label is None:
	LABEL=None
else:
	LABEL=args.label

if args.tag is None:
	TAG = ""
else:
	TAG = "%s-" % args.tag

DATAFILENAME=args.filename
N_FOLDS=args.N
VALIDSPLIT=args.vsplit


###
### Start processing
###

print(f" . Loading '{DATAFILENAME}'")
labelled_data = pd.read_csv(DATAFILENAME)
print(" . Data loaded")
print(" . %d folds: %.2f split for validation" % (N_FOLDS, VALIDSPLIT))

if LABEL is None:
	LABEL = labelled_data.columns[-1]
	print(f" . Using '{LABEL}' as label")
elif LABEL not in labelled_data.columns:
	print(f"Error: Data from '{DATFILENAME}' has no column named '{LABEL}'",
			file=sys.stderr)
	sys.exit(1)


##
## Processing starts here
##



k_folds = StratifiedKFold(n_splits=N_FOLDS)

# Split off the label column into a separate vector
#
# At this point X is a matrix of measures (feature values),
# and y is the vector of labels describing the rows in the X matrix
X, y = labelled_data.drop(columns=[LABEL]), labelled_data[LABEL]


# Now create our data for each fold
fold_index = 0
for train_and_validate_indices, test_indices in k_folds.split(X,y):

	## At this point, we have a chunk of 1/N_FOLDS of the data
	## set aside for testing, however the remaining data is
	## required not only for testing, but if we are going to
	## do any parameter tuning or other validation, we need
	## to further divide this data up so that we have independent
	## estimates of the data distribution for training, validation
	## and testing.

	## We want to preserve our stratification when subdividing the
	## train_and_validation data into training and validation sets,
	## so we use StratifiedShuffleSplit to do this

	# generate a (set of one) split
	splitter = StratifiedShuffleSplit(n_splits=1, test_size=VALIDSPLIT)

	# Pull the data for the train and validate combo.
	# We can access the label values directly by index,
	# but for X we need to select the whole row by index
	# so therefore we use np.take() which does exactly this
	y_train_and_validate = y[train_and_validate_indices]
	X_train_and_validate = np.take(X, train_and_validate_indices, axis=0)

	# we only need the first one, so just use next() to pull the
	# value from the iterator and then discard the iterator
	validation_indices, training_indices = \
			next(splitter.split(
					X_train_and_validate,
					y_train_and_validate))

		
	# At this point we have the indices of data samples
	# for our testing, validation and training -- we just
	# need the data.  Again we pull directly from y, and
	# use np.take() to pull entire rows from X
	y_training = y[training_indices]
	y_validation = y[validation_indices]
	y_test = y[test_indices]

	X_training = np.take(X, training_indices, axis=0)
	X_validation = np.take(X, validation_indices, axis=0)
	X_test = np.take(X, test_indices, axis=0)


	# save the data into the directory identified by this fold
	data_dirname = "%sfolded-%02d" % (TAG, fold_index)
	if not os.path.exists(data_dirname):
		os.mkdir(data_dirname)

	print(f" . Storing fold {fold_index} in {data_dirname}")

	X_training.to_csv("%s/X_train.csv" % data_dirname, index=False)
	X_validation.to_csv("%s/X_validation.csv" % data_dirname, index=False)
	X_test.to_csv("%s/X_test.csv" % data_dirname, index=False)

	y_training.to_csv("%s/y_train.csv" % data_dirname, index=False)
	y_validation.to_csv("%s/y_validation.csv" % data_dirname, index=False)
	y_test.to_csv("%s/y_test.csv" % data_dirname, index=False)

	fold_index += 1

