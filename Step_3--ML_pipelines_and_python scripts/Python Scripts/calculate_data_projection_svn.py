#!/usr/bin/env python3

'''
Standardize the data based on training and apply PCA if requested
'''

import sys
import argparse
import pathlib

import pandas as pd
import numpy as np

import random

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

import matplotlib.pyplot as plt
import seaborn as sns



###
### Functions used below
###
def calculate_PCA_projection( \
            X_training_projected,
            X_validation,
            X_test,
            y_train,
            sufficient_variance_explained_threshold):

    print(" . Performing Dimensionality Reduction using PCA")

    ###
    ### Use PCA to decrease dimensionality as estimated also on the
    ### training data.
    ###


    pca_model = PCA() 
    X_training_projected = pca_model.fit_transform(X_training_projected)

    # In a similar fashion to the scaler, PCA can be "fit" to the data
    # and then used to transform other data sets with the same coordinates
    X_test = pca_model.transform(X_test)
    X_validation = pca_model.transform(X_validation)


    ## Now walk through the components, determining how much of the total
    ## variance is explained as we add each.  The PCA routine presents them
    ## in decreasing order, so the first will explain the most variance,
    ## the second the next most, etc.

    # Calculate the total variance as the sum of the variances in each
    # dimension
    total = sum(pca_model.explained_variance_)


    # Iterate summing each new explained variance until we have decided
    # that we have "enough"
    k=0
    variance_explained=0.0
    current_variance=0.0

    while variance_explained < sufficient_variance_explained_threshold:
        current_variance += pca_model.explained_variance_[k]
        variance_explained = current_variance / total
        k=k+1

    k_PCA_sufficient = k

    print(" . PCA completed: %.2f of variance explained at %d components of %d total"
            % (sufficient_variance_explained_threshold,
                    k_PCA_sufficient, X_n_cols))


    # Using only k_PCA_sufficent components, re-fit the data
    # (producing the same fit) but now transforming to the lower
    # k_PCA_sufficient dimensional space
    pca = PCA(n_components = k_PCA_sufficient)
    X_training_projected = pca.fit_transform(X_training_projected)
    X_test = pca.transform(X_test)
    X_validation = pca.transform(X_validation)


    # Extract the amount of variance explained as a cumulative sum
    # series over the number of components
    cumulative_variance_explained = pca.explained_variance_ratio_.cumsum()



    ###
    ### Create some (relatively arbitrary) headers so that the
    ### the resulting data can be stored as a dataframe
    ###
    new_headers = [ 'PC_%0d' % i for i in range(k_PCA_sufficient) ]

    if MAKE_PCA_FIG:
        print(" . Generating figures....")
        save_pca_explanation_figure(k_PCA_sufficient, variance_explained,
                cumulative_variance_explained, DATA_DIRNAME)
        save_pca_scatterplot_figure(X_training_projected, y_train,
                DATA_DIRNAME)


    return (new_headers, X_training_projected, X_validation, X_test)




def save_pca_explanation_figure(k_PCA_sufficient, \
        variance_explained, \
        cumulative_variance_explained, \
        data_dirname):
        

    '''
    Make a figure to explain PCA results
    '''

    # Make a bar plot -- limit it to only k_PCA_sufficient bars wide
    # rather than the total number of assay components so that it is
    # readable
    plt.bar(range(k_PCA_sufficient), cumulative_variance_explained)

    # Put on some nice labels and title
    plt.ylabel("Cumulative Explained Variance")
    plt.xlabel("Principal Components")
    plt.title("%.2f%% of variance (> 90%%) is explained by the first %d columns"
            % (variance_explained * 100.0, k_PCA_sufficient))

    # save it to a file
    fig_filename = "%s/PCA-variance-explained-svn.pdf" % data_dirname
    plt.savefig(fig_filename, bbox_inches="tight")

    # If you want to see it on the screen, uncomment this -- not that
    # in that case this program will stop here until you close the
    # window
    #plt.show()

def save_pca_scatterplot_figure(X, y, data_dirname):
        
    '''
    Make a figure showing the PCA based scatterplot
    '''


    # y is a 1 x N matrix, and the columns of X are 1 x N vectors,
    # so we have to combine them in two steps
    plotdata = pd.DataFrame({"X_0" : X[:, 0], "X_1" : X[:, 1]})
    plotdata['label'] = y

    # Make a scatter plot of the first two axes of the X data,
    # with the colour (hue) based on the y values
    sns.relplot(x="X_0", y="X_1", hue="label", data=plotdata)

    # Put on some nice labels and title
    plt.xlabel("X_0 Components")
    plt.ylabel("X_1")
    plt.title("Scatter by first to principal components")

    # save it to a file
    fig_filename = "%s/PCA-data-scatter-svn.pdf" % data_dirname
    plt.savefig(fig_filename, bbox_inches="tight")

    # If you want to see it on the screen, uncomment this -- not that
    # in that case this program will stop here until you close the
    # window
    plt.show()


###
### Mainline
###


## get values from command line if given
argparser = argparse.ArgumentParser(
        description="Standardize and project the data using PCA")
argparser.add_argument("-s", "--seed", action="store",
        default=None, type=int,
        help="Set the random seed -- if set to any value, every "
            "run of the program will use the same sequence of "
            "random values")
argparser.add_argument("-P", "--PCA", action="store_true",
        help="Perform PCA based dimensionality reduction. "
            "See also --threshold.")
argparser.add_argument("-t", "--threshold", action="store",
        default=0.90, type=float,
        help="Set the threshold for PCA based dimensionality reduction. "
            "Does nothing unless --PCA is also given.")
argparser.add_argument("-f", "--fig", action="store_true",
        help="Store the results of PCA dimensionality reduction "
            "in a PDF figure")
argparser.add_argument("dirname",
        type=pathlib.Path,
        help="Sets the data directory -- alterate to -I")

args = argparser.parse_args(sys.argv[1:])



DO_PCA = args.PCA
SUFFICIENT_VARIANCE_EXPLAINED = args.threshold
if SUFFICIENT_VARIANCE_EXPLAINED > 1.0 or SUFFICIENT_VARIANCE_EXPLAINED < 0:
    argparser.print_help()
    sys.exit(1)



MAKE_PCA_FIG = args.fig


## set the seed if given for reproducible runs
if args.seed is not None:
    random.seed(argparse.seed)

DATA_DIRNAME=args.dirname


###
### Procesing starts here
###

print(" . Loading....")
X_training_from_file = pd.read_csv("%s/X_train.csv" % DATA_DIRNAME)
X_validation = pd.read_csv("%s/X_validation.csv" % DATA_DIRNAME)
X_test = pd.read_csv("%s/X_test.csv" % DATA_DIRNAME)

y_train = pd.read_csv("%s/y_train.csv" % DATA_DIRNAME)

_, X_n_cols = X_training_from_file.shape


###
###  Perform scaling based on the training data to ensure
###  that some axes are not orders of magnitude bigger than
###  others
###
print(" . Performing Standardized Scaling")

# perform a standard scaling on all X
X_train_scaler = StandardScaler() #The result of standardization is that the features will be rescaled to ensure the mean and the standard deviation to be 0 and 1, respectively
X_training_projected = X_train_scaler.fit_transform(X_training_from_file) #The result of scaling (max-min normalization) is transforming your data so that it fits within a specific scale, like 0-100 or 0-1. You want to scale data when you're using methods based on measures of how far apart data points, like support vector machines, or SVM or k-nearest neighbors, or KNN.
#https://www.kdnuggets.com/2020/04/data-transformation-standardization-normalization.html

# Note that once you have created the scaler using fit_transform,
# you can scale a new data set to the same extent, so this will
# transform the test set into the same coordinate space
X_test = X_train_scaler.transform(X_test)
X_validation = X_train_scaler.transform(X_validation)
#we use fit_transform() on training data but transform() on the test data
#https://towardsdatascience.com/what-and-why-behind-fit-transform-vs-transform-in-scikit-learn-78f915cf96fe

if DO_PCA:
    ## Calculate the PCA transform (below) to reduce dimensions
    (new_headers, X_training_projected, X_validation, X_test) = \
            calculate_PCA_projection(
                    X_training_projected,
                    X_validation,
                    X_test,
                    y_train,
                    SUFFICIENT_VARIANCE_EXPLAINED)
else:

    # create headers that indicate we standardized
    new_headers = [ 'SCALED_%s' % header for header in X_training_from_file.columns ]


# create dataframes from the raw array data created above
df_X_training = pd.DataFrame(data=X_training_projected,
            columns=new_headers)
df_X_validation = pd.DataFrame(data=X_validation, columns=new_headers)
df_X_test = pd.DataFrame(data=X_test, columns=new_headers)


## Store the results
df_X_training.to_csv("%s/X_proj_train.csv" % DATA_DIRNAME, index=False)
df_X_validation.to_csv("%s/X_proj_validation.csv" % DATA_DIRNAME, index=False)
df_X_test.to_csv("%s/X_proj_test.csv" % DATA_DIRNAME, index=False)



