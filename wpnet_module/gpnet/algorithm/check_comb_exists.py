import numpy as np


def build_column_masks(pop_matrix):
    # n, m = matrix.shape
    matrix_T = pop_matrix.astype(bool).T
    col_masks = []
    col_popcnt = []
    for col in matrix_T:
        bytes_arr = np.packbits(col, bitorder="little").tobytes()
        mask = int.from_bytes(bytes_arr, byteorder="little")
        col_masks.append(mask)
        col_popcnt.append(bin(mask).count("1"))

    return col_masks, col_popcnt


def exists_covering_row(query_cols, col_masks, col_popcnt):
    if not query_cols:
        return True
    if not col_masks:
        return True
    if not col_popcnt:
        return True

    sorted_cols = sorted(query_cols, key=lambda x: col_popcnt[x])
    res = col_masks[sorted_cols[0]]
    if not res:
        return False

    for i in sorted_cols[1:]:
        res &= col_masks[i]
        if not res:
            return False

    return True
