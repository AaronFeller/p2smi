#!/usr/bin/env python
"""
Module containing helper functions which can be used by any other module
"""

import os
import operator
import pickle
import pickle
import tempfile
import math

try:
    import oasa
    import pybel
    import openbabel as ob
    import Image as PIL
    import tkinter as tk
    import ImageTk as piltk
except ImportError as e:
    print(e)
    print("Necessary mondule missing for image display function, draw_molecule")


import numpy as np


def get_saved_data(file, folder=""):
    """Opens a python object picked to a file, and returns that object

    Tries a few ways to open pickled files to avoid line ending/binary
    errors"""
    if folder != "":
        file = os.path.join(folder, file)
    try:
        try:
            handle = open(file, "rU")  # U: open in universal newlines
            data = pickle.load(handle)
        except (IndexError, ValueError, ImportError):
            try:
                handle = open(file, "rb")  # b: open as binary
                data = pickle.load(open(file, "rb"))
            except (IndexError, ValueError, ImportError):
                handle = open(file, "r")
                data = pickle.load(handle)
    except IOError:
        print(("Object file %s  not found" % file))
        raise

    handle.close()

    return data


def write_object(obj, file, folder=""):
    """Takes a python object and writes it to a file with cPickle"""
    if folder != "":
        file = os.path.join(folder, file)
    filehandle = open(file, "wb")
    pickle.dump(obj, filehandle, protocol=1)
    filehandle.close()
    return True


def normaliser(input_li, itermax=10):
    """Takes a list of triplets, corresponding to
    (rowname, colname, value),
    and transforms it into a 2D array.

    Normalises the array by rows and columns in turn, itermax times

    returns the list in the same format, but normalised.

    Carries out normalisation of docking scores (predicted affinities) as
    described in Casey et al (2009).

    Receptors (columns) and ligands (rows) are both normalised for. Each column
    (representing a single receptor) is normalised by calculating the mean and
    standard deviation of all docking scores in the column, subtracting the
    stdev from eachs core and dividing by the mean.
    This is repeated for all columns, then all rows, iteratively, a default of
    10 times. The default order is receptors, then ligands

    """
    input_list = list(input_li)  # Take a copy
    input_list.sort(key=operator.itemgetter(0, 1))

    # Build array - get sublist for each rowname. Each colname has already been
    # sorted
    start_array = []
    current_rowname = ""
    row = []
    ordered_rows = []
    ordered_cols = []
    seen_all_cols = 0

    print((input_list[:3], input_list[-3:]))

    for triplet in input_list:
        row_name, column_name, value = triplet

        if row_name not in ordered_rows:
            ordered_rows.append(row_name)
        if column_name not in ordered_cols:
            ordered_cols.append(column_name)

        if not row_name == current_rowname:
            current_rowname = row_name
            if len(row) > 0:  # If row has been filled, will be empty first time
                start_array.append(row)
                row = []
            row.append(value)
        else:
            row.append(value)
    start_array.append(row)

    try:
        array = np.array(start_array)
    except ValueError as e:
        print((len(start_array)))
        for i in start_array:
            print((len(i)))
        for i in start_array:
            for j in i:
                if not isinstance(j, float):
                    print("%s is not a float...")

        raise
    numrows, numcols = array.shape
    print((numrows, numcols))

    # Normalise
    for iteration in range(itermax):
        print(("Iteration:", iteration))
        print("Normalising rows")
        for i in range(numrows):
            mean = np.mean(array[i, :])
            stddev = np.std(array[i, :])
            array[i, :] = np.nan_to_num((array[i, :] - mean) / stddev)
        print("Normalising columns")
        for j in range(numcols):
            mean = np.mean(array[:, j])
            stddev = np.std(array[:, j])
            array[:, j] = np.nan_to_num((array[:, j] - mean) / stddev)

    # Back to list:
    new_out_list = []
    for i, col in enumerate(ordered_cols):
        for j, row in enumerate(ordered_rows):
            new_out_list.append((row, col, array[j, i]))

    return new_out_list
