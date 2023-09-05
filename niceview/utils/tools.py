"""Tools."""


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
    return matrix.tocsr()[:, idx].todense()
