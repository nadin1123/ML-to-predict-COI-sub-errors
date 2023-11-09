
#!/usr/bin/env python3

'''
Format confusion matrices -- either print, or put in a heatmap
'''

import sys
import math

import pandas as pd
import numpy as np

import random

from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from sklearn import metrics

import matplotlib.pyplot as plt
import seaborn as sns


def fix_confusion_class_names_(confusion_size, class_names, y_training):

    # make sure that we have some kind of class names, even if just from y
    if class_names is None or len(class_names) != confusion_size:
        class_names = [ "%d" % i for i in sorted(set(y_training)) ]

    return class_names

def print_confusion_matrix(cm, class_names, y_training, filename=None):
    '''
    Print the confusion matrix out on the screen
    '''

    if filename is not None:
        try:
            ofile = open(filename, "w")
        except IOError as err:
            print("Error opening confusion matrix file '{}' : {}"
                    .format(filename, err), file=sys.stderr)
            raise err
    else:
        ofile = sys.stdout

    confusion_size, _ = cm.shape

    class_names = fix_confusion_class_names_(
            confusion_size, class_names, y_training)

    # find the longest name (for padding)
    name_len = max( [ len(s) for s in class_names ] )

    # update with longest number length
    name_len = max( [ name_len, len('%d' % np.max(cm)) ] )

    # print out the table - first the header of names
    header = "%*s .-> " % (name_len, "")
    for name in class_names:
        header = "%s %*s" % (header, name_len, name)
    print(header, file=ofile)

    # now the data - with names down the left
    for i in range(confusion_size):
        leader = "%*s  :  " % (name_len, class_names[i])
        for j in range(confusion_size):
            leader = "%s %*d" % (leader, name_len, cm[i,j])
        print(leader, file=ofile)


def confusion_matrix_heatmap(filename, cm, class_names, y_training):
    '''
    Store the confusion matrix in a heatmap
    '''

    confusion_size, _ = cm.shape

    class_names = fix_confusion_class_names_(confusion_size, class_names, y_training)


    fig, ax = plt.subplots()
    tick_marks = np.arange(len(class_names))
    plt.xticks(tick_marks, class_names)
    plt.yticks(tick_marks, class_names)

    df_cm = pd.DataFrame(cm,
                         index=['Error', 'Clean'],
                         columns=['Error', 'Clean'])
    sns.heatmap(
            df_cm,
            annot=True,
            cmap="rocket",
            fmt='g')

    ax.xaxis.set_label_position("top")
    plt.tight_layout()
    plt.title('Confusion matrix', y=1.1)
    plt.ylabel('Actual label')
    plt.xlabel('Predicted label')
    fig.savefig(filename, bbox_inches="tight")

    return plt
