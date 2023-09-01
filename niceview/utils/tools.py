"""Tools."""

import numpy as np


def txt_to_list(txt_file):
    """Read lines of a txt file to a list.

    Args:
        txt_file (str): txt file path

    Returns:
        lines (list of str): list of string of lines in the txt file
    """
    with open(txt_file, 'r') as txt:
        lines = txt.readlines()
        lines = [line.strip() for line in lines]
    return lines


def list_to_txt(lines, txt_file):
    """Write lines of a list to a txt file.

    Args:
        lines (list of str): list of string of lines to be written
        txt_file (str): txt file path
    
    Returns:
        txt_file (str): txt file path
    """
    with open(txt_file, 'w') as txt:
        for line in lines:
            txt.write(line)
            txt.write('\n')
    return txt_file


def select_col_from_name(matrix, name_list, name):
    """Select column from matrix by name.
    
    Args:
        matrix (np.ndarray): matrix of shape (row, col).
        name_list (list): list of names.
        name (str): name to select.
    
    Returns:
        np.ndarray: column of shape (row,).
    """
    idx = name_list.index(name)
    return matrix[:, idx]


def normalize_nonzero_to_range(arr, min_range=1, max_range=255, log=False):
    """Normalize log-transformed nonzero elements to a specified range.
    
    Args:
        arr (np.ndarray): array of shape (row, col).
        min_range (int): minimum value of the range.
        max_range (int): maximum value of the range.
        log (bool): whether to log transform the array.
    
    Returns:
        np.ndarray: normalized array of shape (row, col).
    """
    # Find the nonzero elements
    nonzero_elements = arr[arr != 0]
    
    # log transform
    if log:
        nonzero_elements = np.log(nonzero_elements)

    # Calculate the min and max of the nonzero elements
    min_val = np.min(nonzero_elements)
    max_val = np.max(nonzero_elements)

    # Normalize the nonzero elements to the specified range
    normalized_nonzero = (
        (nonzero_elements - min_val) / (max_val - min_val)
    ) * (max_range - min_range) + min_range

    # Create a new array with the same shape as the original array
    normalized_arr = np.zeros_like(arr, dtype=float)
    normalized_arr[arr != 0] = normalized_nonzero

    return normalized_arr
