#include "rawfile.c"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parsefinnigan.c"
#include "parsemzml.c"

// Shared entrypoints

RawFile* FromBytes(char* data, int profile, int centroid, char* format) {
  int16_t success;
  RawFile* rawFile = malloc(sizeof(RawFile));

  rawFile->numScans = 0;
  rawFile->scans = NULL;
  rawFile->numChroms = 0;
  rawFile->chromatograms = NULL;
  rawFile->nFilters = 0;
  rawFile->scanFilters = NULL;

  rawFile->fileName = NULL;
  rawFile->instrumentModel = NULL;

  rawFile->timestamp = 0;
  rawFile->data = data;
  if (strcmp(format, "rawfile") == 0) success = initializeRawfile(rawFile, profile, centroid);
  if (strcmp(format, "mzml") == 0) success = initializeMzml(rawFile, profile, centroid);
  if (success < 0) {
    free(rawFile);
    return NULL;
  }
  return rawFile;
}

RawFile* Open(char* filename, int profile, int centroid, char* format, uint64_t size) {
  // printf("\nOpen called***********************.\n");
  char* data = copyData(filename, size);
  // printf("opening...");
  RawFile* rawFile = FromBytes(data, profile, centroid, format);
  free(data);
  return rawFile;
}

int16_t main() {
  char* filename = "D:\\Work\\TFE Tests\\CN Small\\Agilent-6495-QQQ.mzML";
  RawFile* rawFile;
  rawFile = Open(filename, 0, 1, "mzml", 0);
  Xic results;
  // for (int i = 0; i < rawFile->nFilters; i++) {
  //   if (i == 10) break;
  //   ScanFilter f = rawFile->scanFilters[i];
  //   printf("%d %d %f %f\n", i, f.scanMode, f.highMass, f.lowMass);
  // }
  printf("opened2\n");
  results = Pull_chrom_xic(rawFile, 6.8, 8.8, 90.06, 44.296, 1);
  printf("opened3 %d\n", results.length);
  for (int i = 0; i < results.length; i++) {
    printf("%f ", results.intensity[i]);
  }
}
