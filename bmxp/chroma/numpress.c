
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct CharResult {
  int64_t len;
  unsigned char* data;

} CharResult;

typedef struct DoubleResult {
  int64_t len;
  double* data;

} DoubleResult;
const int IS_LITTLE_ENDIAN = 1;
const int THROW_ON_OVERFLOW = 1;

double optimalSlofFixedPoint(const double* data, size_t dataSize) {
  if (dataSize == 0) return 0;

  double maxDouble = 1;
  double x;
  double fp;

  for (size_t i = 0; i < dataSize; i++) {
    x = log(data[i] + 1);
    maxDouble = fmax(maxDouble, x);
  }

  fp = floor(0xFFFF / maxDouble);

  return fp;
}

void encodeFixedPoint(double fixedPoint, unsigned char* result) {
  int i;
  unsigned char* fp = (unsigned char*)&fixedPoint;
  for (i = 0; i < 8; i++) {
    result[i] = fp[IS_LITTLE_ENDIAN ? (7 - i) : i];
  }
}

double decodeFixedPoint(const unsigned char* data) {
  int i;
  double fixedPoint;
  unsigned char* fp = (unsigned char*)&fixedPoint;

  for (i = 0; i < 8; i++) {
    fp[i] = data[IS_LITTLE_ENDIAN ? (7 - i) : i];
  }
  return fixedPoint;
}

static void encodeInt(const int64_t x, unsigned char* res, size_t* res_length) {
  // get the bit pattern of a signed int x_inp
  unsigned int m;
  unsigned char i, l; // numbers between 0 and 9

  unsigned int mask = 0xf0000000;
  int64_t init = x & mask;

  if (init == 0) {
    l = 8;
    for (i = 0; i < 8; i++) {
      m = mask >> (4 * i);
      if ((x & m) != 0) {
        l = i;
        break;
      }
    }
    res[0] = l;
    for (i = l; i < 8; i++) {
      res[1 + i - l] = (unsigned char)(x >> (4 * (i - l)));
    }
    *res_length += 1 + 8 - l;

  } else if (init == mask) {
    l = 7;
    for (i = 0; i < 8; i++) {
      m = mask >> (4 * i);
      if ((x & m) != m) {
        l = i;
        break;
      }
    }
    res[0] = l + 8;
    for (i = l; i < 8; i++) {
      res[1 + i - l] = (unsigned char)(x >> (4 * (i - l)));
    }
    *res_length += 1 + 8 - l;

  } else {
    res[0] = 0;
    for (i = 0; i < 8; i++) {
      res[1 + i] = (unsigned char)(x >> (4 * i));
    }
    *res_length += 9;
  }
}

static void decodeInt(const unsigned char* data, size_t* di, size_t max_di, size_t* half, unsigned int* res) {
  size_t n, i;
  unsigned int mask, m;
  unsigned char head;
  unsigned char hb;

  // Extract the first half byte, specifying the number of leading zero half
  // bytes of the final integer.
  // If half is zero, we look at the first half byte, otherwise we look at
  // the second (lower) half byte and advance the counter to the next char.
  if (*half == 0) {
    head = data[*di] >> 4;
  } else {
    head = data[*di] & 0xf;
    (*di)++;
  }

  *half = 1 - (*half); // switch to other half byte
  *res = 0;

  if (head <= 8) {
    n = head;
  } else { // we have n leading ones, fill n half bytes in res with 0xf
    n = head - 8;
    mask = 0xf0000000;
    for (i = 0; i < n; i++) {
      m = mask >> (4 * i);
      *res = *res | m;
    }
  }

  if (n == 8) {
    return;
  }

  if (*di + ((8 - n) - (1 - *half)) / 2 >= max_di) {
    printf("[MSNumpress::decodeInt] Corrupt input data!\n");
    *res = -1;
    return;
  }

  for (i = n; i < 8; i++) {
    if (*half == 0) {
      hb = data[*di] >> 4;
    } else {
      hb = data[*di] & 0xf;
      (*di)++;
    }
    *res = *res | (unsigned int)(hb) << ((i - n) * 4);
    *half = 1 - (*half);
  }
}

double optimalLinearFixedPoint(const double* data, size_t dataSize) {
  /*
   * safer impl - apparently not needed though
   *
  if (dataSize == 0) return 0;

  double maxDouble = 0;
  double x;

  for (size_t i=0; i<dataSize; i++) {
          x = data[i];
          maxDouble = max(maxDouble, x);
  }

  return floor(0xFFFFFFFF / maxDouble);
  */
  if (dataSize == 0) return 0;
  if (dataSize == 1) return floor(0x7FFFFFFFl / data[0]);
  double maxDouble = fmax(data[0], data[1]);
  double extrapol;
  double diff;

  for (size_t i = 2; i < dataSize; i++) {
    extrapol = data[i - 1] + (data[i - 1] - data[i - 2]);
    diff = data[i] - extrapol;
    maxDouble = fmax(maxDouble, ceil(abs(diff) + 1));
  }

  return floor(0x7FFFFFFFl / maxDouble);
}

double optimalLinearFixedPointMass(const double* data, size_t dataSize, double mass_acc) {
  if (dataSize < 3) return 0; // we just encode the first two points as floats

  // We calculate the maximal fixedPoint we need to achieve a specific mass
  // accuracy. Note that the maximal error we will make by encoding as int is
  // 0.5 due to rounding errors.
  double maxFP = 0.5 / mass_acc;

  // There is a maximal value for the FP given by the int length (32bit)
  // which means we cannot choose a value higher than that. In case we cannot
  // achieve the desired accuracy, return failure (-1).
  double maxFP_overflow = optimalLinearFixedPoint(data, dataSize);
  if (maxFP > maxFP_overflow) return -1;

  return maxFP;
}

CharResult encodeLinear(const double* data, int64_t dataSize, double fixedPoint) {
  // given data, length, and fixedPoint,
  CharResult result;
  int64_t ints[3];
  size_t i, ri;
  unsigned char halfBytes[10];
  size_t halfByteCount;
  size_t hbi;
  int64_t extrapol;
  int64_t diff;
  result.len = dataSize * 5 + 8;
  result.data = malloc(result.len * 8);

  // allocate maximum of the result

  // printf("Encoding %d doubles with fixed point %f\n", (int)dataSize, fixedPoint);
  encodeFixedPoint(fixedPoint, result.data);

  if (dataSize == 0) {
    result.len = 8;
    return result;
  }

  ints[1] = (int64_t)(data[0] * fixedPoint + 0.5);
  for (i = 0; i < 4; i++) {
    result.data[8 + i] = (ints[1] >> (i * 8)) & 0xff;
  }

  if (dataSize == 1) {
    result.len = 12;
    return result;
  }

  ints[2] = (int64_t)(data[1] * fixedPoint + 0.5);
  for (i = 0; i < 4; i++) {
    result.data[12 + i] = (ints[2] >> (i * 8)) & 0xff;
  }

  halfByteCount = 0;
  ri = 16;

  for (i = 2; i < dataSize; i++) {
    ints[0] = ints[1];
    ints[1] = ints[2];
    if (THROW_ON_OVERFLOW && data[i] * fixedPoint + 0.5 > LLONG_MAX) {
      printf("[MSNumpress::encodeLinear] Next number overflows LLONG_MAX.");
      result.len = -1;
      return result;
    }

    ints[2] = (int64_t)(data[i] * fixedPoint + 0.5);
    extrapol = ints[1] + (ints[1] - ints[0]);

    if (THROW_ON_OVERFLOW && (ints[2] - extrapol > INT_MAX || ints[2] - extrapol < INT_MIN)) {
      printf("[MSNumpress::encodeLinear] Cannot encode a number that exceeds the bounds of [-INT_MAX, INT_MAX].");
      result.len = -1;
      return result;
    }

    diff = (ints[2] - extrapol);
    // printf("%lu %lu %lu,   extrapol: %ld    diff: %d \n", ints[0], ints[1], ints[2], extrapol, diff);
    encodeInt((diff), &halfBytes[halfByteCount], &halfByteCount);
    /*
    printf("%d (%d):  ", diff, (int)halfByteCount);
    for (size_t j=0; j<halfByteCount; j++) {
            printf("%x ", halfBytes[j] & 0xf);
    }
    printf("\n");
    */

    for (hbi = 1; hbi < halfByteCount; hbi += 2) {
      result.data[ri] = (unsigned char)((halfBytes[hbi - 1] << 4) | (halfBytes[hbi] & 0xf));
      // printf("%x \n", result[ri]);
      ri++;
    }
    if (halfByteCount % 2 != 0) {
      halfBytes[0] = halfBytes[halfByteCount - 1];
      halfByteCount = 1;
    } else {
      halfByteCount = 0;
    }
  }
  if (halfByteCount == 1) {
    result.data[ri] = (unsigned char)(halfBytes[0] << 4);
    ri++;
  }
  result.len = ri;
  return result;
}

DoubleResult decodeLinear(const unsigned char* data, const size_t dataSize) {
  size_t i;
  size_t ri = 0;
  int64_t init, buff;
  int diff;
  int64_t ints[3];
  // double d;
  size_t di;
  size_t half;
  int64_t extrapol;
  int64_t y;
  double fixedPoint;
  DoubleResult result;
  result.data = malloc(dataSize * 2 * 8);

  // printf("Decoding %d bytes with fixed point %f\n", (int)dataSize, fixedPoint);

  if (dataSize == 8) {
    result.len = 0;
    return result;
  };

  if (dataSize < 8) {
    printf("[MSNumpress::decodeLinear] Corrupt input data: not enough bytes to read fixed point!\n");
    result.len = -1;
    return result;
  }

  fixedPoint = decodeFixedPoint(data);

  if (dataSize < 12) {
    printf("[MSNumpress::decodeLinear] Corrupt input data: not enough bytes to read first value!\n");
    result.len = -1;
    return result;
  }

  ints[1] = 0;
  for (i = 0; i < 4; i++) {
    ints[1] = ints[1] | ((0xff & (init = data[8 + i])) << (i * 8));
  }
  result.data[0] = ints[1] / fixedPoint;

  if (dataSize == 12) {
    result.len = 1;
    return result;
  }
  if (dataSize < 16) {
    printf("[MSNumpress::decodeLinear] Corrupt input data: not enough bytes to read second value!\n");
    result.len = -1;
    return result;
  }

  ints[2] = 0;
  for (i = 0; i < 4; i++) {
    ints[2] = ints[2] | ((0xff & (init = data[12 + i])) << (i * 8));
  }
  result.data[1] = ints[2] / fixedPoint;

  half = 0;
  ri = 2;
  di = 16;

  // printf("   di     ri      half    int[0]    int[1]    extrapol   diff\n");

  while (di < dataSize) {
    if (di == (dataSize - 1) && half == 1) {
      if ((data[di] & 0xf) == 0x0) {
        break;
      }
    }
    // printf("%7d %7d %7d %lu %lu %ld", di, ri, half, ints[0], ints[1], extrapol);

    ints[0] = ints[1];
    ints[1] = ints[2];
    decodeInt(data, &di, dataSize, &half, (unsigned int*)&buff);
    diff = (int)(buff);

    extrapol = ints[1] + (ints[1] - ints[0]);
    y = extrapol + diff;
    // printf(" %d \n", diff);
    result.data[ri++] = y / fixedPoint;
    ints[2] = y;
  }
  result.len = ri;
  return result;
}

CharResult encodePic(const double* data, size_t dataSize) {
  size_t i, ri;
  int64_t x;
  unsigned char halfBytes[10];
  size_t halfByteCount;
  size_t hbi;

  CharResult result;
  result.data = malloc(dataSize * 5 * 8);

  // printf("Encoding %d doubles\n", (int)dataSize);

  halfByteCount = 0;
  ri = 0;

  for (i = 0; i < dataSize; i++) {

    if (THROW_ON_OVERFLOW && (data[i] + 0.5 > INT_MAX || data[i] < -0.5)) {
      printf("[MSNumpress::encodePic] Cannot use Pic to encode a number larger than INT_MAX or smaller than 0.");
      result.len = -1;
    }
    x = (int64_t)(data[i] + 0.5);
    // printf("%d %d %d,   extrapol: %d    diff: %d \n", ints[0], ints[1], ints[2], extrapol, diff);
    encodeInt(x, &halfBytes[halfByteCount], &halfByteCount);

    for (hbi = 1; hbi < halfByteCount; hbi += 2) {
      result.data[ri] = (unsigned char)((halfBytes[hbi - 1] << 4) | (halfBytes[hbi] & 0xf));
      // printf("%x \n", result[ri]);
      ri++;
    }
    if (halfByteCount % 2 != 0) {
      halfBytes[0] = halfBytes[halfByteCount - 1];
      halfByteCount = 1;
    } else {
      halfByteCount = 0;
    }
  }
  if (halfByteCount == 1) {
    result.data[ri] = (unsigned char)(halfBytes[0] << 4);
    ri++;
  }
  result.len = ri;
  return result;
}

DoubleResult decodePic(const unsigned char* data, const size_t dataSize) {
  size_t ri;
  unsigned int x;
  size_t di;
  size_t half;

  DoubleResult result;
  result.data = malloc(dataSize * 2 * 8);
  half = 0;
  ri = 0;
  di = 0;

  while (di < dataSize) {
    if (di == (dataSize - 1) && half == 1) {
      if ((data[di] & 0xf) == 0x0) {
        break;
      }
    }

    decodeInt(data, &di, dataSize, &half, &x);
    if (x == -1) {
      result.len = -1;
      return result;
    }

    // printf("%7d %7d %7d %7d %d\n", ri, di, half, dataSize, x);

    // printf("count: %d \n", count);
    result.data[ri++] = (double)(x);
  }
  result.len = ri;
  return result;
}

CharResult encodeSlof(const double* data, size_t dataSize, double fixedPoint) {
  CharResult result;
  size_t i, ri;
  double temp;
  unsigned short x;
  result.data = malloc(dataSize * 5 + 8);
  encodeFixedPoint(fixedPoint, result.data);

  ri = 8;
  for (i = 0; i < dataSize; i++) {
    temp = log(data[i] + 1) * fixedPoint;

    if (THROW_ON_OVERFLOW && temp > USHRT_MAX) {
      result.len = -1;
      printf("[MSNumpress::encodeSlof] Cannot encode a number that overflows USHRT_MAX.");
      return result;
    }

    x = (unsigned short)(temp + 0.5);
    result.data[ri++] = x & 0xff;
    result.data[ri++] = (x >> 8) & 0xff;
  }
  result.len = ri;
  return result;
}

DoubleResult decodeSlof(const unsigned char* data, const size_t dataSize) {
  DoubleResult result;
  size_t i, ri;
  unsigned short x;
  double fixedPoint;
  result.data = malloc(dataSize * 2 * 8);
  decodeFixedPoint(data);

  if (dataSize < 8) {
    printf("[MSNumpress::decodeSlof] Corrupt input data: not enough bytes to read fixed point!\n");
    result.len = -1;
    return result;
  }

  ri = 0;
  fixedPoint = decodeFixedPoint(data);

  for (i = 8; i < dataSize; i += 2) {
    x = (unsigned short)(data[i] | (data[i + 1] << 8));
    result.data[ri++] = exp(x / fixedPoint) - 1;
  }
  result.len = ri;
  return result;
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
// Tests

void tOptimalLinearFixedPoint() {
  srand(123459);

  size_t n = 1000;
  double mzs[1000];
  mzs[0] = 300 + (rand() % 1000) / 1000.0;
  for (size_t i = 1; i < n; i++)
    mzs[i] = mzs[i - 1] + (rand() % 1000) / 1000.0;
  printf("+                    max val: %f\n", mzs[n - 1]);
  printf("+          optimal linear fp: %f\n", optimalLinearFixedPoint(mzs, n));
  printf("+ pass    optimalLinearFixedPoint \n\n");
}

void tOptimalLinearFixedPointMass() {
  srand(123459);

  size_t n = 1000;
  double mzs[1000];
  mzs[0] = 300 + (rand() % 1000) / 1000.0;
  for (size_t i = 1; i < n; i++)
    mzs[i] = mzs[i - 1] + (rand() % 1000) / 1000.0;

  double fixedP;

  fixedP = optimalLinearFixedPointMass(mzs, n, 0.1);
  assert(abs(5.0 - fixedP) < 0.000005);

  fixedP = optimalLinearFixedPointMass(mzs, n, 0.01);
  assert(abs(50.0 - fixedP) < 0.000005);

  fixedP = optimalLinearFixedPointMass(mzs, n, 0.001);
  assert(abs(500.0 - fixedP) < 0.000005);

  printf("+          optimal linear fp: %f\n", fixedP);
  printf("+ pass    optimalLinearFixedPointMass \n\n");
}

void tEncodeLinear1() {
  size_t nMzs = 1;
  double mzs[1];
  mzs[0] = 100.0;
  CharResult encodedBytes = encodeLinear(mzs, nMzs, 100000.0);

  assert(12 == encodedBytes.len);
  assert(0x80 == encodedBytes.data[8]);
  assert(0x96 == encodedBytes.data[9]);
  assert(0x98 == encodedBytes.data[10]);
  assert(0x00 == encodedBytes.data[11]);

  printf("+ pass    encodeLinear1\n\n");
  free(encodedBytes.data);
}

void tEncodeLinear4() {
  size_t nMzs = 1;
  double mzs[4];

  mzs[0] = 100.0;
  mzs[1] = 200.0;
  mzs[2] = 300.00005;
  mzs[3] = 400.00010;

  nMzs = 4;
  unsigned char encoded[20];
  CharResult encodedBytes = encodeLinear(mzs, nMzs, 100000.0);

  assert(18 == encodedBytes.len);
  assert(0x80 == encodedBytes.data[8]);
  assert(0x96 == encodedBytes.data[9]);
  assert(0x98 == encodedBytes.data[10]);
  assert(0x00 == encodedBytes.data[11]);
  assert(0x75 == encodedBytes.data[16]);
  assert(0x80 == encodedBytes.data[17]);
  printf("+ pass    encodeLinear\n\n");
  free(encodedBytes.data);
}

void tDecodeLinearNice() {

  double mzs[4];

  mzs[0] = 100.0;
  mzs[1] = 200.0;
  mzs[2] = 300.00005;
  mzs[3] = 400.00010;

  size_t nMzs = 4;
  double fixedPoint = optimalLinearFixedPoint(mzs, nMzs);
  CharResult encodedBytes = encodeLinear(mzs, nMzs, fixedPoint);

  DoubleResult decodedBytes = decodeLinear(encodedBytes.data, encodedBytes.len);
  assert(4 == decodedBytes.len);
  assert(abs(100.0 - decodedBytes.data[0]) < 0.000005);
  assert(abs(200.0 - decodedBytes.data[1]) < 0.000005);
  assert(abs(300.00005 - decodedBytes.data[2]) < 0.000005);
  assert(abs(400.00010 - decodedBytes.data[3]) < 0.000005);
  printf("+ pass    decodeLinearNice \n");

  free(encodedBytes.data);
  free(decodedBytes.data);
}

void tDecodeLinearNiceLowFP() {

  double mzs[7];

  mzs[0] = 100.0;
  mzs[1] = 200.0;
  mzs[2] = 300.00005;
  mzs[3] = 400.00010;
  mzs[4] = 450.00010;
  mzs[5] = 455.00010;
  mzs[6] = 700.00010;

  size_t nMzs = 7;
  // check for fixed points
  {
    double fixedPoint;
    fixedPoint = optimalLinearFixedPointMass(mzs, nMzs, 0.1);
    assert(abs(5 - fixedPoint) < 0.000005);

    fixedPoint = optimalLinearFixedPointMass(mzs, nMzs, 1e-3);
    assert(abs(500 - fixedPoint) < 0.000005);

    fixedPoint = optimalLinearFixedPointMass(mzs, nMzs, 1e-5);
    assert(abs(50000 - fixedPoint) < 0.000005);

    fixedPoint = optimalLinearFixedPointMass(mzs, nMzs, 1e-7);
    assert(abs(5000000 - fixedPoint) < 0.000005);

    // cannot fulfill accuracy of 1e-8
    fixedPoint = optimalLinearFixedPointMass(mzs, nMzs, 1e-8);
    assert(abs(-1 - fixedPoint) < 0.000005);
  }

  {
    double fixedPoint = optimalLinearFixedPointMass(mzs, nMzs, 0.001);
    CharResult encodedBytes = encodeLinear(mzs, nMzs, fixedPoint);

    DoubleResult decodedBytes = decodeLinear(encodedBytes.data, encodedBytes.len);
    assert(25 == encodedBytes.len);
    assert(7 == decodedBytes.len);

    assert(abs(100.0 - decodedBytes.data[0]) < 0.001);
    assert(abs(200.0 - decodedBytes.data[1]) < 0.001);
    assert(abs(300.00005 - decodedBytes.data[2]) < 0.001);
    assert(abs(400.00010 - decodedBytes.data[3]) < 0.001);
    free(encodedBytes.data);
    free(decodedBytes.data);
  }

  double mz_err[5];
  double encodedLength[5];

  // for higher accuracy, we get longer encoded lengths
  mz_err[0] = 0.1;
  encodedLength[0] = 22;
  mz_err[1] = 1e-3;
  encodedLength[1] = 25;
  mz_err[2] = 1e-5;
  encodedLength[2] = 29;
  mz_err[3] = 1e-6;
  encodedLength[3] = 30;
  mz_err[4] = 1e-7;
  encodedLength[4] = 31;

  for (int k = 0; k < 4; k++) {
    double fixedPoint = optimalLinearFixedPointMass(mzs, nMzs, mz_err[k]);
    CharResult encodedBytes = encodeLinear(mzs, nMzs, fixedPoint);
    DoubleResult decodedBytes = decodeLinear(encodedBytes.data, encodedBytes.len);
    assert(encodedLength[k] == encodedBytes.len);
    assert(7 == decodedBytes.len);
    assert(abs(100.0 - decodedBytes.data[0]) < mz_err[k]);
    assert(abs(200.0 - decodedBytes.data[1]) < mz_err[k]);
    assert(abs(300.00005 - decodedBytes.data[2]) < mz_err[k]);
    assert(abs(400.00010 - decodedBytes.data[3]) < mz_err[k]);
    assert(abs(450.00010 - decodedBytes.data[4]) < mz_err[k]);
    assert(abs(455.00010 - decodedBytes.data[5]) < mz_err[k]);
    assert(abs(700.00010 - decodedBytes.data[6]) < mz_err[k]);
    free(encodedBytes.data);
    free(decodedBytes.data);
  }

  printf("+ pass    decodeLinearNiceFP \n\n");
}

void tDecodeLinearWierd() {
  double mzs[4];

  mzs[0] = 100.0;
  mzs[1] = 200.0;
  mzs[2] = 300.00005;
  mzs[3] = 0.00010;

  size_t nMzs = 4;
  double fixedPoint = optimalLinearFixedPoint(mzs, nMzs);
  CharResult encodedBytes = encodeLinear(mzs, nMzs, fixedPoint);
  DoubleResult decodedBytes = decodeLinear(encodedBytes.data, encodedBytes.len);
  assert(4 == decodedBytes.len);
  assert(abs(100.0 - decodedBytes.data[0]) < 0.000005);
  assert(abs(200.0 - decodedBytes.data[1]) < 0.000005);
  assert(abs(300.00005 - decodedBytes.data[2]) < 0.000005);
  assert(abs(0.00010 - decodedBytes.data[3]) < 0.000005);

  printf("+ pass    decodeLinearWierd \n\n");
  free(encodedBytes.data);
  free(decodedBytes.data);
}

void tDecodeLinearCorrupt1() {
  unsigned char encoded[20];
  DoubleResult decodedBytes = decodeLinear(encoded, 11);
  if (decodedBytes.len >= 0) {
    printf("- fail    decodeLinearCorrupt1: didn't throw exception for corrupt input \n\n");
    assert(0 == 1);
  }
  free(decodedBytes.data);
  decodedBytes = decodeLinear(encoded, 14);
  if (decodedBytes.len >= 0) {
    printf("- fail    decodeLinearCorrupt1: didn't throw exception for corrupt input \n\n");
    assert(0 == 1);
  }

  printf("+ pass    decodeLinearCorrupt 1 \n\n");
  free(decodedBytes.data);
}

void tDecodeLinearCorrupt2() {

  double mzs[4];

  mzs[0] = 100.0;
  mzs[1] = 200.0;
  mzs[2] = 300.00005;
  mzs[3] = 0.00010;

  size_t nMzs = 4;
  double fixedPoint = optimalLinearFixedPoint(mzs, nMzs);
  CharResult encodedBytes = encodeLinear(mzs, nMzs, fixedPoint);
  DoubleResult decodedBytes = decodeLinear(encodedBytes.data, 11);
  if (decodedBytes.len >= 0) {
    printf("- fail    decodeLinearCorrupt2: didn't throw exception for corrupt input \n\n");
    assert(0 == 1);
  }

  printf("+ pass    decodeLinearCorrupt 2 \n\n");
  free(encodedBytes.data);
  free(decodedBytes.data);
}

void tEncodeDecodeLinearStraight() {
  size_t n = 15;
  double mzs[15];
  for (size_t i = 0; i < n; i++)
    mzs[i] = i;

  double fixedPoint = optimalLinearFixedPoint(mzs, n);
  CharResult encodedBytes = encodeLinear(mzs, n, fixedPoint);
  DoubleResult decodedBytes = decodeLinear(encodedBytes.data, encodedBytes.len);

  assert(n == decodedBytes.len);

  double m = 0;
  double error;
  double mLim = 0.000005;

  for (size_t i = 0; i < n; i++) {
    error = abs(mzs[i] - decodedBytes.data[i]);
    m = fmax(m, error);
    if (error >= mLim) {
      printf("error   %f above limit %f\n", error, mLim);
      assert(error < mLim);
    }
  }
  printf("+     size compressed: %ld%%\n", encodedBytes.len / (n * 8) * 100);
  printf("+           max error: %f   limit: %f\n", m, mLim);
  printf("+ pass    encodeDecodeLinearStraight \n\n");
  free(encodedBytes.data);
  free(decodedBytes.data);
}

void tEncodeDecodeLinear() {
  srand(123459);

  size_t n = 1000;
  double mzs[1000];
  mzs[0] = 300 + rand() / RAND_MAX;
  for (size_t i = 1; i < n; i++)
    mzs[i] = mzs[i - 1] + rand() / RAND_MAX;

  double fixedPoint = optimalLinearFixedPoint(mzs, n);
  CharResult encodedBytes = encodeLinear(mzs, n, fixedPoint);
  DoubleResult decodedBytes = decodeLinear(encodedBytes.data, encodedBytes.len);

  assert(n == decodedBytes.len);

  double m = 0;
  double error;
  double mLim = 0.000005;

  for (size_t i = 0; i < n; i++) {
    error = abs(mzs[i] - decodedBytes.data[i]);
    m = fmax(m, error);
    if (error >= mLim) {
      printf("error   %f above limit %f\n", error, mLim);
      assert(error < mLim);
    }
  }
  printf("+     size compressed: %ld%%\n", encodedBytes.len / (n * 8) * 100);
  printf("+           max error: %f   limit: %f\n", m, mLim);
  printf("+ pass    encodeDecodeLinear \n\n");
  free(encodedBytes.data);
  free(decodedBytes.data);
}

void tEncodeDecodePic() {
  srand(123459);

  size_t n = 1000;
  double ics[1000];
  for (size_t i = 0; i < n; i++)
    ics[i] = rand() % 1000000;

  CharResult encodedBytes = encodePic(ics, n);
  DoubleResult decodedBytes = decodePic(encodedBytes.data, encodedBytes.len);

  assert(n == decodedBytes.len);

  for (size_t i = 0; i < n; i++)
    assert(abs(ics[i] - decodedBytes.data[i]) < 0.5);

  printf("+ pass    encodeDecodePic \n");
  free(encodedBytes.data);
  free(decodedBytes.data);
}

void tEncodeDecodePic5() {
  srand(123459);

  size_t n = 1000;
  double ics[1000];
  for (size_t i = 0; i < n; i++)
    ics[i] = rand() % 1000000;

  ics[1] = 0.0;
  ics[2] = 0.0001;
  ics[3] = 0.0002;
  DoubleResult firstDecoded;
  firstDecoded.len = n;
  firstDecoded.data = malloc(firstDecoded.len * 8);
  CharResult encodedBytes = encodePic(ics, n);
  DoubleResult decodedBytes = decodePic(encodedBytes.data, encodedBytes.len);
  for (size_t i = 0; i < n; i++)
    firstDecoded.data[i] = decodedBytes.data[i];
  for (size_t i = 0; i < 5; i++) {
    free(encodedBytes.data);
    encodedBytes = encodePic(decodedBytes.data, decodedBytes.len);
    free(decodedBytes.data);
    decodedBytes = decodePic(encodedBytes.data, encodedBytes.len);
  }

  assert(n == decodedBytes.len);

  for (size_t i = 0; i < n; i++)
    if (firstDecoded.data[i] != decodedBytes.data[i]) {
      printf("%f %f\n", firstDecoded.data[i], decodedBytes.data[i]);
      assert(firstDecoded.data[i] == decodedBytes.data[i]);
    }

  printf("+ pass    encodeDecodePic5 \n");
  free(encodedBytes.data);
  free(decodedBytes.data);
  free(firstDecoded.data);
}

void tEncodeDecodeLinear5() {
  srand(123662);

  size_t n = 1000;
  double mzs[1000];
  mzs[0] = 100 + (rand() % 1000) / 1000.0;
  for (size_t i = 1; i < n; i++)
    mzs[i] = mzs[i - 1] + (rand() % 1000) / 1000.0;

  DoubleResult firstDecoded;
  firstDecoded.len = n;
  firstDecoded.data = malloc(firstDecoded.len * 8);

  double fixedPoint = optimalLinearFixedPoint(mzs, n);

  CharResult encodedBytes = encodeLinear(mzs, n, fixedPoint);
  DoubleResult decodedBytes = decodeLinear(encodedBytes.data, encodedBytes.len);

  for (size_t i = 0; i < n; i++)
    firstDecoded.data[i] = decodedBytes.data[i];

  for (size_t i = 0; i < 5; i++) {
    free(encodedBytes.data);
    encodedBytes = encodeLinear(decodedBytes.data, decodedBytes.len, fixedPoint);
    free(decodedBytes.data);
    decodedBytes = decodeLinear(encodedBytes.data, encodedBytes.len);
  }

  double afterFixedPoint = optimalLinearFixedPoint(decodedBytes.data, n);

  assert(fixedPoint == afterFixedPoint);
  assert(n == decodedBytes.len);

  for (size_t i = 0; i < n; i++)
    if (firstDecoded.data[i] != decodedBytes.data[i]) {
      printf("%f %f\n", firstDecoded.data[i], decodedBytes.data[i]);
      assert(firstDecoded.data[i] == decodedBytes.data[i]);
    }

  printf("+ pass    encodeDecodeLinear5 \n");
  free(encodedBytes.data);
  free(decodedBytes.data);
  free(firstDecoded.data);
}

void tOptimalSlofFixedPoint() {

  srand(123459);

  size_t n = 1000;
  double ics[1000];
  for (size_t i = 0; i < n; i++)
    ics[i] = rand() % 1000000;

  printf("+           optimal slof fp: %f\n", optimalSlofFixedPoint(ics, n));
  printf("+ pass    optimalSlofFixedPoint \n");
}

void tEncodeDecodeSlof() {
  srand(123459);

  size_t n = 1000;
  double ics[1000];
  for (size_t i = 0; i < n; i++)
    ics[i] = rand() % 1000000;

  ics[1] = 0.0;
  ics[2] = 0.0001;
  ics[3] = 0.0002;
  double fixedPoint = optimalSlofFixedPoint(ics, n);
  CharResult encodedBytes = encodeSlof(ics, n, fixedPoint);
  DoubleResult decodedBytes = decodeSlof(encodedBytes.data, encodedBytes.len);
  assert(n == decodedBytes.len);

  double m = 0;
  double rm = 0;
  double mLim = 0.0005;
  double rmLim = 0.0005;
  double error;

  for (size_t i = 0; i < n; i++)
    if (ics[i] < 1.0) {
      error = abs(ics[i] - decodedBytes.data[i]);
      m = fmax(m, error);
      if (error >= mLim) {

        printf("\n%f %f\n", ics[i], decodedBytes.data[i]);
        assert(error < mLim);
      }
    } else {
      error = abs((ics[i] - decodedBytes.data[i]) / ((ics[i] + decodedBytes.data[i]) / 2));
      rm = fmax(rm, error);
      if (error >= rmLim) {
        printf("\n%f %f\n", ics[i], decodedBytes.data[i]);
        assert(error < rmLim);
      }
    }
  printf("+               max error: %f  limit: %f\n", m, mLim);
  printf("+           max rel error: %f  limit: %f\n", rm, rmLim);
  printf("+ pass    encodeDecodeSlof \n\n");
  free(encodedBytes.data);
  free(decodedBytes.data);
}

void tEncodeDecodeSlof5() {
  srand(123459);

  size_t n = 1000;
  double ics[1000];
  for (size_t i = 0; i < n; i++)
    ics[i] = rand() % 1000000;

  ics[1] = 0.0;
  ics[2] = 0.0001;
  ics[3] = 0.0002;

  double fixedPoint = optimalSlofFixedPoint(ics, n);

  CharResult encodedBytes = encodeSlof(ics, n, fixedPoint);
  DoubleResult decodedBytes = decodeSlof(encodedBytes.data, encodedBytes.len);
  DoubleResult firstDecoded;
  firstDecoded.len = n;
  firstDecoded.data = malloc(firstDecoded.len * 8);

  for (size_t i = 0; i < n; i++)
    firstDecoded.data[i] = decodedBytes.data[i];

  for (size_t i = 0; i < 5; i++) {
    free(encodedBytes.data);
    encodedBytes = encodeSlof(decodedBytes.data, n, fixedPoint);
    free(decodedBytes.data);
    decodedBytes = decodeSlof(encodedBytes.data, encodedBytes.len);
  }

  double afterFixedPoint = optimalSlofFixedPoint(decodedBytes.data, n);

  assert(n == decodedBytes.len);
  assert(fixedPoint == afterFixedPoint);

  for (size_t i = 0; i < n; i++)
    if (firstDecoded.data[i] != decodedBytes.data[i]) {
      printf("\n%f %f", firstDecoded.data[i], decodedBytes.data[i]);
      assert(firstDecoded.data[i] == decodedBytes.data[i]);
    }

  printf("+ pass    encodeDecodeSlof5 \n");
  free(firstDecoded.data);
  free(encodedBytes.data);
  free(decodedBytes.data);
}

void testErroneousDecodePic() {
  // set data to [  100, 102, 140, 92, 33, 80, 145  ]; // Base64 is "ZGaMXCFQkQ=="
  int64_t n = 32;
  unsigned char* data = calloc(n, 1);
  data[0] = 100;
  data[1] = 102;
  data[2] = 140;
  data[3] = 92;
  data[4] = 33;
  data[5] = 80;
  data[6] = 145;

  DoubleResult result = decodePic(data, n);
  if (result.len >= 0) {
    printf("- fail    testErroneousDecodePic: didn't throw exception for corrupt input %ld \n\n", result.len);
    assert(0 == 1);
  }

  printf("+ pass    testErroneousDecodePic \n\n");
  free(data);
  free(result.data);
}

// void main() {
//   tOptimalLinearFixedPoint();
//   tOptimalLinearFixedPointMass();
//   tEncodeLinear1();
//   tEncodeLinear4();
//   tDecodeLinearNice();
//   tDecodeLinearNiceLowFP();
//   tDecodeLinearWierd();
//   tDecodeLinearCorrupt1();
//   tDecodeLinearCorrupt2();
//   tEncodeDecodeLinearStraight();
//   tEncodeDecodeLinear();
//   tEncodeDecodePic();
//   //   encodeDecodeSafeStraight();
//   //   encodeDecodeSafe();
//   tOptimalSlofFixedPoint();
//   tEncodeDecodeSlof();
//   tEncodeDecodeLinear5();
//   tEncodeDecodePic5();
//   tEncodeDecodeSlof5();
//   testErroneousDecodePic();
// }

// /////////////////////////////////////////////////////////////
// Not sure how to implement the "safe" stuff in C

// int64_t encodeSafe(const double* data, const size_t dataSize, unsigned char* result) {
//   size_t i, j, ri = 0;
//   double latest[3];
//   double extrapol, diff;
//   const unsigned char* fp;
//   CharResult result;

//   // printf("d0 d1 d2 extrapol diff\n");

//   if (dataSize == 0) {

//     return ri;
//   }

//   latest[1] = data[0];
//   fp = (unsigned char*)data;
//   for (i = 0; i < 8; i++) {
//     result[ri++] = fp[IS_LITTLE_ENDIAN ? (7 - i) : i];
//   }

//   if (dataSize == 1) return ri;

//   latest[2] = data[1];
//   fp = (unsigned char*)&(data[1]);
//   for (i = 0; i < 8; i++) {
//     result[ri++] = fp[IS_LITTLE_ENDIAN ? (7 - i) : i];
//   }

//   fp = (unsigned char*)&diff;
//   for (i = 2; i < dataSize; i++) {
//     latest[0] = latest[1];
//     latest[1] = latest[2];
//     latest[2] = data[i];
//     extrapol = latest[1] + (latest[1] - latest[0]);
//     diff = latest[2] - extrapol;
//     // printf("%f %f %f %f %f\n", latest[0], latest[1], latest[2], extrapol, diff);
//     for (j = 0; j < 8; j++) {
//       result[ri++] = fp[IS_LITTLE_ENDIAN ? (7 - j) : j];
//     }
//   }

//   return ri;
// }

// int64_t decodeSafe(const unsigned char* data, const size_t dataSize, double* result) {
//   size_t i, di, ri;
//   double extrapol, diff;
//   double latest[3];
//   unsigned char* fp;

//   if (dataSize % 8 != 0)
//     printf("[MSNumpress::decodeSafe] Corrupt input data: number of bytes needs to be multiple of 8! ");
//   return -1;

//   // printf("d0 d1 extrapol diff\td2\n");

//   fp = (unsigned char*)&(latest[1]);
//   for (i = 0; i < 8; i++) {
//     fp[i] = data[IS_LITTLE_ENDIAN ? (7 - i) : i];
//   }
//   result[0] = latest[1];

//   if (dataSize == 8) return 1;

//   fp = (unsigned char*)&(latest[2]);
//   for (i = 0; i < 8; i++) {
//     fp[i] = data[8 + (IS_LITTLE_ENDIAN ? (7 - i) : i)];
//   }
//   result[1] = latest[2];

//   ri = 2;

//   fp = (unsigned char*)&diff;
//   for (di = 16; di < dataSize; di += 8) {
//     latest[0] = latest[1];
//     latest[1] = latest[2];

//     for (i = 0; i < 8; i++) {
//       fp[i] = data[di + (IS_LITTLE_ENDIAN ? (7 - i) : i)];
//     }

//     extrapol = latest[1] + (latest[1] - latest[0]);
//     latest[2] = extrapol + diff;

//     // printf("%f %f %f %f\t%f \n", latest[0], latest[1], extrapol, diff, latest[2]);

//     result[ri++] = latest[2];
//   }
//   return ri;
// }
