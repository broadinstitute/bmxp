"""
Correlation calculations in python; a fallback in case of a missing binary file or
unsupported operating system
"""

import numpy as np
from scipy import stats


def backfill(a, b):
    """Backfill NaNs in each array based on the non-NaN values in the other. wraps
    _backfill for more convenient use

    :param a: array-like, first data array
    :param b: array-like, second data array
    """
    a, b = np.nan_to_num((a, b))
    a_filled = _backfill(a, b)
    b_filled = _backfill(b, a)  # pylint: disable=arguments-out-of-order

    a_filled[a_filled == 0] = np.nan
    b_filled[b_filled == 0] = np.nan

    return a_filled, b_filled


def _backfill(a, b):
    """Backfill NaNs in a based on non-NaN values in b

    :param a: array-like, first data array, array to backfill
    :param b: array-like, second data array
    """
    zeroes = len(a) - np.count_nonzero(a)
    fill_limit = min(np.count_nonzero(a), np.count_nonzero(b))
    if zeroes == 0 or fill_limit == 0:  # nothing to fill
        return a

    sorted_index = np.lexsort((b, a))
    fill_value = -1

    for i in range(zeroes, max(zeroes - fill_limit, 0), -1):
        index = sorted_index[i - 1]
        if b[index] == 0:  # out of non-NaN values
            break
        a[index] = fill_value
        fill_value -= 1

    return a


def free_p(_):
    """Empty function to match function name of C implementation"""


def pearson(x, y, _):
    """Pearson correlation with NaNs dropped

    :param x: array-like, first data array
    :param y: array-like, second data array
    """
    not_nan = np.logical_and(~np.isnan(x), ~np.isnan(y))
    x = x[not_nan]
    y = y[not_nan]
    return np.corrcoef(x, y)[0][1]


def spearman(x, y, _, drop_nan, legacy_mode):
    """Spearman correlation calculation

    :param x: array-like, first data array
    :param y: array-like, second data array
    :param drop_nan: boolean, when False, missing values are backfilled
    :param legacy_mode: boolean, when True, NaNs are filled with 0s when correlating
    features with very few overlapping non-NaN values
    """
    if legacy_mode and np.logical_and(~np.isnan(x), ~np.isnan(y)).sum() < 10:
        # ensure we don't modify the original array
        x = np.copy(x)
        y = np.copy(y)
        x[np.isnan(x)] = 0
        y[np.isnan(y)] = 0
    if not drop_nan:
        x, y = backfill(x, y)
    not_nan = np.logical_and(~np.isnan(x), ~np.isnan(y))
    x = x[not_nan]
    y = y[not_nan]

    # this is about 2x faster than scipy.stats.spearmanr
    rank_x = stats.rankdata(x)
    rank_y = stats.rankdata(y)
    return np.corrcoef(rank_x, rank_y)[0][1]


def pearson_array(arr, r1_start, r2_start, size):
    """Pearson method that accepts a 2D array and relevant indexes

    :param arr: 2D array-like, full set of data
    :param r1_start: int, index of first data array
    :param r2_start: int, index of second data array
    :param size: unused
    """
    return pearson(arr[r1_start], arr[r2_start], size)


def spearman_array(arr, r1_start, r2_start, size, drop_nan, legacy_mode):
    """Spearman method that accepts a 2D array and relevant indexes

    :param arr: 2D array-like, full set of data
    :param r1_start: int, index of first data array
    :param r2_start: int, index of second data array
    :param size: unused
    :param drop_nan: boolean, when False, missing values are backfilled
    :param legacy_mode: boolean, when True, NaNs are filled with 0s when correlating
    features with very few overlapping non-NaN values
    """
    return spearman(arr[r1_start], arr[r2_start], size, drop_nan, legacy_mode)
