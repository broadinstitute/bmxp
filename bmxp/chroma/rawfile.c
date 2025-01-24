#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct Xic {
  float* rt;
  float* intensity;
  int32_t length;
  float* mzs;
} Xic;

typedef struct Scan {
  float* prIntensities;
  float* prMzs;
  float* centIntensities;
  float* centMzs;
  uint32_t prTotal;   // granular total of all mz/ints
  uint32_t centTotal; // delete// number of centroided peaks.
  double time;
  uint32_t index;
  uint32_t filter;

} Scan;

typedef struct Chromatogram {
  char* id;
  float* intensities;
  float* rts;
  uint32_t length; // granular total of all mz/ints
  float precursorMz;
  float productMz;
} Chromatogram;

typedef struct ScanFilter {
  uint8_t polarity; // 0 neg, 1 pos
  uint8_t analyzer; // which analyzer
  uint8_t msLevel;  // level of ms, i.e. how many fragmentations per scan
  uint32_t nPrecursors;
  double* precursorMz;
  double* energy;
  double lowMass;
  double highMass;
  uint8_t scanMode; // 0 - no profile, 1 - has profile
  uint8_t scanType; // 0:Full, 1:Zoom, 2:SIM, 3:SRM, 4:CRM, 5:undefined, 6:Q1, 7:Q3

  // future:
  // uint32_t numCollisions;    // how many collisions per each mslevel (mslevel-1)
  // float** collisionEnergies; // list of collision energies
  // uint32_t numIsolations;    // len MSLevel
  // float*** isolations;       // don't fear the triple pointer
} ScanFilter;

typedef struct RawFile {
  char* instrumentModel;
  char* fileName;
  uint32_t nControllers;
  Scan* scans;
  char* data;
  uint32_t numScans;
  uint32_t FirstScanNumber;
  uint32_t LastScanNumber;
  double startTime;
  double endTime;
  ScanFilter* scanFilters;
  uint32_t nFilters;
  uint64_t timestamp;
  Chromatogram* chromatograms;
  uint32_t numChroms;
} RawFile;

uint8_t compareFilters(ScanFilter* a, ScanFilter* b) {
  if (a->polarity != b->polarity || a->scanMode != b->scanMode || a->msLevel != b->msLevel ||
      a->scanType != b->scanType || a->analyzer != b->analyzer || a->nPrecursors != b->nPrecursors ||
      a->lowMass != b->lowMass || a->highMass != b->highMass) {
    return 0;
  }
  if (a->msLevel > 1) {
    for (uint8_t i = 0; i < a->nPrecursors; i++) {
      if (a->precursorMz[i] != b->precursorMz[i] || a->energy[i] != b->energy[i]) {
        return 0;
      }
    }
  }
  return 1;
}
uint64_t file_size(char* filename) {
  // return size;
  uint64_t filesize = 0; /* or unsigned long for C89 compatability*/
  unsigned char buffer[1048576];
  FILE* f = fopen(filename, "rb");
  while (1) {
    size_t numRead = fread(buffer, 1, 1048576, f);
    if (numRead >= 1)
      filesize += numRead;
    else
      break;
  }
  fclose(f);
  return filesize;
}

char* copyData(char* filename, uint64_t size) {
  if (size == 0) size = file_size(filename);
  FILE* f = fopen(filename, "rb");
  fseek(f, 0, SEEK_END);
  char* fileArray = (char*)malloc(sizeof(char) * size); // this becomes RawFile->data; freed 482
  rewind(f);
  fread(fileArray, 1, size, f);
  fclose(f);
  return fileArray;
}

int findFirstGE(Scan* arr, uint32_t n, float target) // replace with binary
{
  if (arr[0].time >= target) {
    return 0;
  }
  if (n == 1) {
    return 1;
  }
  uint32_t n2 = n / 2;
  if (arr[n2].time >= target) {
    return findFirstGE(arr, n2, target);
  }
  return n2 + findFirstGE(&arr[n2], n - n2, target);
}

int findLastLE(Scan* arr, uint32_t n, float target) // replace with binary
{
  if (arr[0].time > target) {
    return -1;
  }
  if (n == 1) {
    return 0;
  }
  uint32_t n2 = n / 2;
  if (arr[n2 - 1].time <= target) {
    return n2 + findLastLE(&arr[n2], n - n2, target);
  }
  return findLastLE(arr, n2, target);
}

int findFirstGEFloat(float* arr, uint32_t n, float target) // replace with binary
{
  if (n == 0 || arr[0] >= target) {
    return 0;
  }
  if (n == 1) {
    return 1;
  }
  uint32_t n2 = n / 2;
  if (arr[n2] >= target) {
    return findFirstGEFloat(arr, n2, target);
  }
  return n2 + findFirstGEFloat(&arr[n2], n - n2, target);
}

int findLastLEFloat(float* arr, uint32_t n, float target) // replace with binary
{
  if (n == 0 || arr[0] > target) {
    return -1;
  }
  if (n == 1) {
    return 0;
  }
  uint32_t n2 = n / 2;
  if (arr[n2 - 1] <= target) {
    return n2 + findLastLEFloat(&arr[n2], n - n2, target);
  }
  return findLastLEFloat(arr, n2, target);
}

float float_min(float* arr, uint32_t len) {

  if (len <= 0 || arr == NULL) {
    return 0;
  }
  float result = arr[0];
  for (uint32_t i = 1; i < len; i++)
    if (arr[i] < result) result = arr[i];

  return result;
}

float float_max(float* arr, uint32_t len) {

  if (len <= 0 || arr == NULL) {
    return 0;
  }
  float result = arr[0];
  for (uint32_t i = 1; i < len; i++)
    if (arr[i] > result) result = arr[i];

  return result;
}

Xic Pull_xic(RawFile* rawFile, float rt1, float rt2, float mz, float ppm, int32_t filter, int32_t centroid,
             int pullMzs) {
  // filter = index of target scan filter, or -1 for no scan filter
  Scan* scans = rawFile->scans;

  int32_t start = findFirstGE(scans, rawFile->numScans, rt1);
  int32_t stop = findLastLE(scans, rawFile->numScans, rt2);
  int32_t length = stop - start + 1;
  float* intensities = (float*)calloc(length, sizeof(float));
  float* rts = (float*)malloc(length * sizeof(float));
  float* mzs = NULL;
  if (pullMzs) {
    mzs = (float*)malloc(length * (sizeof(float)));
  }
  int32_t peakStart = 0;
  int32_t peakStop = 0;

  int32_t pos = 0;
  for (int32_t i = 0; i < length; i++) {
    if (filter < 0 || (int)scans[start + i].filter == filter) {
      rts[pos] = (float)scans[start + i].time;
      // iterate through mzs and add them up
      float mzAvgs = 0;
      float newIntensity;
      float* mzsPtrStart;
      float* intensityPtrStart;
      int32_t length;
      if (centroid) {
        mzsPtrStart = scans[start + i].centMzs;
        intensityPtrStart = scans[start + i].centIntensities;
        length = scans[start + i].centTotal;
      } else {
        mzsPtrStart = scans[start + i].prMzs;
        intensityPtrStart = scans[start + i].prIntensities;
        length = scans[start + i].prTotal;
      }
      // get total intensity and mz weighted avg
      peakStart = findFirstGEFloat(mzsPtrStart, length, mz - (mz / 1000000) * ppm);
      peakStop = findLastLEFloat(mzsPtrStart, length, mz + (mz / 1000000) * ppm);
      for (int32_t j = peakStart; j < peakStop + 1; j++) {
        intensities[pos] += intensityPtrStart[j];
        if (pullMzs) mzAvgs += intensityPtrStart[j] * mzsPtrStart[j];
      }
      if (pullMzs)
        if (intensities[pos] > 0)
          mzs[pos] = mzAvgs / intensities[pos];
        else
          mzs[pos] = 0;
      pos++;
    }
  }
  return (Xic){rts, intensities, pos, mzs};
}

Xic Pull_chrom_xic(RawFile* rawFile, float rt1, float rt2, float precursor, float product, int pullMzs) {
  Chromatogram* chrom;
  // find the right chromatogram
  for (uint32_t i = 0; i < rawFile->numChroms; i++) {
    chrom = &(rawFile->chromatograms[i]);
    if ((fabs(chrom->precursorMz - precursor) <= .5) && (fabs(chrom->productMz - product) <= .5) &&
        (rt2 >= float_min(chrom->rts, chrom->length)) && (rt1 <= float_max(chrom->rts, chrom->length))) {
      break;
    } else
      chrom = NULL;
  }

  if (chrom == NULL) {
    return (Xic){NULL, NULL, 0, NULL};
  }

  int32_t start = findFirstGEFloat(chrom->rts, chrom->length, rt1);
  int32_t stop = findLastLEFloat(chrom->rts, chrom->length, rt2);
  int32_t length = stop - start + 1;
  float* intensities = (float*)calloc(length, sizeof(float));
  float* rts = (float*)calloc(length, sizeof(float));
  float* mzs = NULL;
  if (pullMzs) {
    mzs = (float*)calloc(length, (sizeof(float)));
    for (uint32_t i = 0; i < length; i++)
      mzs[i] = product;
  }
  memcpy(rts, &(chrom->rts[start]), length * 4);
  memcpy(intensities, &(chrom->intensities[start]), length * 4);

  return (Xic){rts, intensities, length, mzs};
}

int freeFilter(ScanFilter* filter) {
  free(filter->precursorMz);
  free(filter->energy);
  return 0;
}

int freeScan(Scan* scan) {
  free(scan->prIntensities);
  free(scan->prMzs);
  free(scan->centIntensities);
  free(scan->centMzs);
  return 0;
}

int freeChrom(Chromatogram* chrom) {
  free(chrom->id);
  free(chrom->rts);
  free(chrom->intensities);
  return 0;
}

int16_t FreeRawFile(RawFile* rf) {
  free(rf->fileName);
  free(rf->instrumentModel);
  for (uint32_t i = 0; i < rf->nFilters; i++) {
    freeFilter(&rf->scanFilters[i]);
  }
  free(rf->scanFilters);
  for (uint32_t i = 0; i < rf->numScans; i++) {
    freeScan(&rf->scans[i]);
  }
  free(rf->scans);

  for (uint32_t i = 0; i < rf->numChroms; i++) {
    freeChrom(&rf->chromatograms[i]);
  }
  free(rf->chromatograms);

  // rf->data  already freed
  free(rf);
}

int16_t FreeP(void* p) { free(p); }
