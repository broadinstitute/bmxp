#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct Ranking {
  double value;
  int32_t index;
} Ranking;

double sum_m(double* x, int32_t size, int* mask) {
  double total = 0;
  for (int32_t i = 0; i < size; i++) {
    if (mask[i]) total += x[i];
  }
  return total;
}

double covariance_m(double* x, double* y, double meanX, double meanY, int32_t size, int* mask) {
  double result = 0;

  for (int32_t i = 0; i < size; i++) {
    if (mask[i]) result += (x[i] - meanX) * (y[i] - meanY);
  }
  return result;
}

double stdev_m(double* x, double meanX, int32_t size, int* mask, int32_t mask_size) {
  double stdevSquared = 0;
  for (int32_t i = 0; i < size; i++) {
    if (mask[i]) stdevSquared += pow(x[i] - meanX, 2);
  }
  stdevSquared = stdevSquared / mask_size;

  return sqrt(stdevSquared);
}

int comparator(const void* a, const void* b) {
  double valA = ((struct Ranking*)a)->value;
  double valB = ((struct Ranking*)b)->value;
  if (valA > valB) {
    return 1;
  }
  if (valA < valB) {
    return -1;
  }
}

int32_t* nanMask(double* x, double* y, int32_t size) {
  int32_t* mask = malloc(sizeof(int32_t) * size);
  for (int32_t i = 0; i < size; i++) {
    mask[i] = (x[i] == x[i]) && (y[i] == y[i]); // Both must be not nan
  }
  return mask;
}

double pearson(double* x, double* y, int32_t size) {
  int32_t* mask = nanMask(x, y, size);
  int32_t mask_size = 0;
  for (int32_t i = 0; i < size; i++) {
    mask_size += mask[i];
  }
  double meanX = sum_m(x, size, mask) / mask_size;
  double meanY = sum_m(y, size, mask) / mask_size;
  double cov = covariance_m(x, y, meanX, meanY, size, mask) / mask_size;
  double stdevX = stdev_m(x, meanX, size, mask, mask_size);
  double stdevY = stdev_m(y, meanY, size, mask, mask_size);
  free(mask);
  return cov / (stdevX * stdevY);
}

void free_p(void* ptr) { free(ptr); }

double spearman(double* x, double* y, int32_t size, int32_t dropNan) {
  Ranking* sortedX = malloc(sizeof(Ranking) * size);
  Ranking* sortedY = malloc(sizeof(Ranking) * size);
  int32_t* eitherNan = malloc(sizeof(int32_t) * size);
  double* xRanks = malloc(sizeof(double) * size);
  double* yRanks = malloc(sizeof(double) * size);
  double* xCopy = NULL;
  double* yCopy = NULL;
  int32_t i = 0;
  int32_t j = 0;

  if (dropNan) {
    xCopy = malloc(sizeof(double) * size);
    yCopy = malloc(sizeof(double) * size);

    // copy x and y so we don't modify original values, then point "x" and "y" pointers to the copies for ease of use
    for (int32_t i = 0; i < size; i++) {
      memcpy(&xCopy[i], &x[i], 8);
      memcpy(&yCopy[i], &y[i], 8);
    }
    x = xCopy;
    y = yCopy;

    // match behavior of old clustering algorithm
    int32_t totalNansToDrop = 0;
    for (i = 0; i < size; i++) {
      totalNansToDrop += (x[i] != x[i]) || (y[i] != y[i]);
    }
    if (size - totalNansToDrop < 10) {
      for (i = 0; i < size; i++) {
        if (y[i] != y[i]) y[i] = 0;
        if (x[i] != x[i]) x[i] = 0;
      }
    }
  }

  for (int32_t i = 0; i < size; i++) {
    memcpy(&xRanks[i], &x[i], 8);
    memcpy(&yRanks[i], &y[i], 8);
  }

  int32_t x_pos = 0;
  int32_t y_pos = 0;
  int32_t y_nan = -1;
  int32_t x_nan = -1;
  int32_t x_is_nan = 0;
  int32_t y_is_nan = 0;
  for (i = 0; i < size; i++) {
    // move to back
    if (x[i] != x[i]) {
      sortedX[size + x_nan].value = NAN;
      sortedX[size + x_nan].index = i;
      x_nan--;
    } else {
      memcpy(&sortedX[x_pos].value, &x[i], 8);
      sortedX[x_pos].index = i;
      x_pos++;
    }
    if (y[i] != y[i]) {
      sortedY[size + y_nan].value = NAN;
      sortedY[size + y_nan].index = i;
      y_nan--;
    } else {
      memcpy(&sortedY[y_pos].value, &y[i], 8);
      sortedY[y_pos].index = i;
      y_pos++;
    }
  }

  qsort(sortedX, x_pos, sizeof(Ranking), comparator); // index and values of sorted x
  qsort(sortedY, y_pos, sizeof(Ranking), comparator); // index and values of sorted y

  // Iterate through X starting at highest non-nan. Record the index of Ys to fill in, highest number to lowest
  // Track how many are filled, update nan-mask
  // x, y, and mask are in original order
  int32_t fill_limit = 0;
  if (x_pos < y_pos) {
    fill_limit = x_pos;
  } else {
    fill_limit = y_pos;
  }
  int32_t num_ys_to_fill = 0;
  int32_t num_xs_to_fill = 0;
  int32_t* y_idx_to_fill = malloc(sizeof(int32_t) * fill_limit);
  int32_t* x_idx_to_fill = malloc(sizeof(int32_t) * fill_limit);

  int32_t index;
  double value;
  int32_t next_x = 1;
  int32_t next_y = 1;
  // iterate backwards through sortedx to find Ys to fill
  if (!dropNan) {
    for (i = x_pos; (i > 0) && (num_ys_to_fill < fill_limit); i--) {
      index = sortedX[i - 1].index;
      value = sortedX[i - 1].value;
      if (y[index] != y[index]) { // if Y is nan, and if we have remaining Nans to fill
        // the order of y indices to fill. y_idx_to_fill = [68, 48, 49]; y[68] = 3, y[48] = 2, y[49]=1
        y_idx_to_fill[num_ys_to_fill] = index;
        // mask[index] = 1;
        num_ys_to_fill++;
      }
    }
    for (i = y_pos; (i > 0) && (num_xs_to_fill < fill_limit); i--) {
      index = sortedY[i - 1].index;
      value = sortedY[i - 1].value;
      if (x[index] != x[index]) { // if Y is nan, and if we have remaining Nans to fill
        // the order of y indices to fill. y_idx_to_fill = [68, 48, 49]; y[68] = 3, y[48] = 2, y[49]=1
        x_idx_to_fill[num_xs_to_fill] = index;
        // mask[index] = 1;
        num_xs_to_fill++;
      }
    }

    // iterate backwards through x_idx_to_fill, the last entry is 1
    for (i = num_xs_to_fill; i > 0; i--) {
      xRanks[x_idx_to_fill[i - 1]] = next_x;
      next_x++;
    }
    // iterate backwards through x_idx_to_fill
    for (i = num_ys_to_fill; i > 0; i--) {
      yRanks[y_idx_to_fill[i - 1]] = next_y;
      next_y++;
    }
  }

  // find where we will fill in remaining
  for (i = 0; i < size; i++) {
    eitherNan[i] = (yRanks[i] != yRanks[i]) || (xRanks[i] != xRanks[i]);
  }

  // find positions where neither x or y are NAN
  double currentValue;
  double nextValue;
  double tieRank;
  int32_t numToFill;
  int32_t fillIndex;
  // fill in remaining X's, where there is a real value in y
  for (i = 0; i < x_pos; i++) {
    index = sortedX[i].index;
    if (eitherNan[index]) {
      continue;
    }
    numToFill = 0;
    nextValue = NAN;
    currentValue = x[sortedX[i].index];
    // get the num to fill
    for (j = 0; j < x_pos - i; j++) {
      if (eitherNan[sortedX[i + j].index]) {
        continue;
      }
      nextValue = x[sortedX[i + j].index];
      if (nextValue != currentValue) {
        break;
      } else {
        numToFill++;
      }
    }
    tieRank = next_x + (numToFill - 1) / 2.0;
    // fill the values
    int filled = 0;
    for (j = 0; filled < numToFill; j++) {
      fillIndex = sortedX[i + j].index;
      if (eitherNan[fillIndex]) {
        continue;
      }
      filled += 1;
      xRanks[fillIndex] = tieRank;
    }
    i += j - 1;
    next_x += numToFill;
  }

  // fill in remaining X's, where there is a real value in y
  for (i = 0; i < y_pos; i++) {
    index = sortedY[i].index;
    if (eitherNan[index]) {
      continue;
    }
    numToFill = 0;
    nextValue = NAN;
    currentValue = y[sortedY[i].index];
    for (j = 0; j < y_pos - i; j++) {
      if (eitherNan[sortedY[i + j].index]) {
        continue;
      }
      nextValue = y[sortedY[i + j].index];
      if (nextValue != currentValue) {
        break;
      } else {
        numToFill++;
      }
    }
    tieRank = next_y + (numToFill - 1) / 2.0;
    int filled = 0;
    for (j = 0; filled < numToFill; j++) {
      fillIndex = sortedY[i + j].index;
      if (eitherNan[fillIndex]) {
        continue;
      }
      filled += 1;
      yRanks[fillIndex] = tieRank;
    }
    i += j - 1;
    next_y += numToFill;
  }

  double result = pearson(xRanks, yRanks, size);
  free(sortedX);
  free(sortedY);
  free(xRanks);
  free(yRanks);
  free(eitherNan);
  free(y_idx_to_fill);
  free(x_idx_to_fill);
  free(xCopy);
  free(yCopy);
  return result;
}

double spearman_array(double* arr, int32_t r1_start, int32_t r2_start, int32_t size, int32_t dropNan) {
  return spearman(&arr[r1_start * size], &arr[r2_start * size], size, dropNan);
}

double pearson_array(double* arr, int32_t r1_start, int32_t r2_start, int32_t size) {
  return pearson(&arr[r1_start * size], &arr[r2_start * size], size);
}
