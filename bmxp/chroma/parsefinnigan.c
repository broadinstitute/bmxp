#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

char* errorString;

char* GetErrorString() { return errorString; }

typedef struct ProfileChunk {
  uint32_t Firstbin;
  uint32_t Nbins; //   uint32
  float Fudge;    //
  float* Signal;  //   []float32
} ProfileChunk;

typedef struct FinniganInfo {
  double lowMz;            // finnigan
  double highMz;           // finnigan
  uint64_t ScanindexAddr;  // finnigan
  uint64_t ScanparamsAddr; // finnigan
  uint64_t dataAddr;       // finnigan
  uint64_t scanTrailerAddr;
} FinniganInfo;

typedef struct FinniganScan {
  ProfileChunk* chunks;  // delete
  double FirstValue;     // delete
  double Step;           // delete
  uint32_t PeakCount;    // delete// number of chunks (peaks?)
  uint64_t offset;       // delete
  double A;              // delete
  double B;              // delete
  double C;              // delete
  uint32_t PeaklistSize; // delete// uint32
  uint32_t Layout;       // delete// uint32
  uint32_t Nparam;
  uint32_t Nbins;       // uint32
  uint32_t ProfileSize; // uint32
  unsigned char* flags;

} FinniganScan;

char* fromPascal(char* data, uint64_t pos) {
  uint64_t length = 0;
  memcpy(&length, &data[pos], 4);
  pos += 4;
  char* returnString = (char*)malloc(length + 1);
  for (uint64_t i = 0; i < length; i++) {
    memcpy(&returnString[i], &data[pos], 1);
    pos += 2;
  }
  returnString[length] = '\0';
  return returnString;
}

uint64_t movePascal(char* data, uint64_t pos, uint64_t max) {
  uint64_t length = 0;
  for (uint64_t i = 0; i < max; i++) {
    memcpy(&length, &data[pos], 4);
    pos += 4 + length * 2;
  }
  return pos;
}

int16_t fillScanEvents(RawFile* rawFile, FinniganInfo finniganInfo, FinniganScan* finniganScans) {
  uint64_t pos = finniganInfo.scanTrailerAddr;
  uint32_t k = 0;

  ScanFilter* filters = (ScanFilter*)malloc(sizeof(ScanFilter) * rawFile->numScans); // freed 467

  rawFile->scanFilters = filters;
  Scan* scans = rawFile->scans;
  char* data = rawFile->data;
  pos += 4;
  for (uint32_t i = 0; i < rawFile->numScans; i++) {

    // printf("Parsing %d\n", rawFile->numScans);
    // if (i % 1000 == 0) {
    //   printf("Scan: %d ", i);
    // }
    // printf("filter pos: %lld ", pos + 4);
    memcpy(&filters[k].polarity, &data[pos + 4], 1);
    // printf("polarity: %d ", filters[k].polarity);
    memcpy(&filters[k].scanMode, &data[pos + 5], 1);
    // printf("scanMode: %d ", filters[k].scanMode);
    memcpy(&filters[k].msLevel, &data[pos + 6], 1);
    // printf("msLevel: %d ", filters[k].msLevel);
    memcpy(&filters[k].scanType, &data[pos + 7], 1);
    uint8_t inSource = data[pos + 8];
    // printf("scanType: %d ", filters[k].scanType);
    // memcpy(&filters[k].dependentScans, &data[pos + 10], 1);
    memcpy(&filters[k].analyzer, &data[pos + 40], 1);
    // printf("filters[k].analyzer: %d ", filters[k].analyzer);
    pos += 136;
    memcpy(&filters[k].nPrecursors, &data[pos], 4);
    // printf("nprecursors: %ld ", filters[k].nPrecursors);
    pos += 4;
    filters[k].precursorMz = (double*)calloc(filters[k].nPrecursors, sizeof(double)); // freed 448
    filters[k].energy = (double*)calloc(filters[k].nPrecursors, sizeof(double));      // freed 449
    for (uint8_t n = 0; n < filters[k].nPrecursors; n++) {
      if (filters[k].msLevel > 1) {
        memcpy(&filters[k].precursorMz[n], &data[pos], 8);
        memcpy(&filters[k].energy[n], &data[pos + 16], 8);
      }
      pos += 56;
    }

    memcpy(&filters[k].lowMass, &data[pos + 4], 8);
    memcpy(&filters[k].highMass, &data[pos + 12], 8);
    // printf("low high: %f %f ", filters[k].lowMass, filters[k].highMass);
    int32_t match = -1;
    for (uint32_t fn = 0; fn < k; fn++) {
      if (compareFilters(&filters[k], &filters[fn])) {
        match = fn;
        break;
      }
    }
    if (match > -1) {
      scans[i].filter = match;
      // free precursor mz/energy lists
      // all other fields will be overwritten on the next loop
      free(filters[k].precursorMz);
      free(filters[k].energy);
    } else {
      scans[i].filter = k;
      k++;
    }
    memcpy(&finniganScans[i].Nparam, &data[pos + 20], 4);
    // printf("Nparam: %d ", finniganScans[i].Nparam);
    pos += 24;
    if (finniganScans[i].Nparam > 0) {
      memcpy(&finniganScans[i].A, &data[pos + 16], 8);
      memcpy(&finniganScans[i].B, &data[pos + 24], 8);
      memcpy(&finniganScans[i].C, &data[pos + 32], 8);
    }
    // printf("Nparam: %s ", rawFile->instrumentModel);
    // imaging files fail here. They need 52 bytes for some reason
    // pos += 52;
    if (strcmp(rawFile->instrumentModel, "Orbitrap Exploris 240") == 0 && inSource == 1) {
      pos += 48;
      pos = movePascal(rawFile->data, pos, 1);
    } else if (strcmp(rawFile->instrumentModel, "") == 0) {
      pos += 52;
    } else if (finniganScans[i].Nparam == 7) {
      pos += 68;
    } else if (finniganScans[i].Nparam == 5) {
      pos += 60;
    } else if (finniganScans[i].Nparam == 0) {
      pos += 12;
    }
  }

  rawFile->nFilters = k;
}

int16_t fillScanIndices(RawFile* rawFile, FinniganInfo finniganInfo, FinniganScan* finniganScans) {
  uint64_t pos = finniganInfo.ScanindexAddr;
  char* data = rawFile->data;
  Scan* scans = rawFile->scans;
  for (uint32_t i = 0; i < rawFile->numScans; i++) {
    memcpy(&scans[i].index, &data[pos + 4], 4);
    memcpy(&scans[i].time, &data[pos + 24], 8);
    memcpy(&finniganScans[i].offset, &data[pos + 72], 8);
    pos += 88;
  }
}

int16_t readScanDataPacket(RawFile* rawFile, uint64_t index, Scan* scan, int profile, FinniganInfo finniganInfo,
                           FinniganScan* finniganScan) {
  char* data = rawFile->data;
  finniganScan->flags = NULL;

  memcpy(&finniganScan->ProfileSize, &data[index + 4], 4);
  memcpy(&finniganScan->PeaklistSize, &data[index + 8], 4);
  memcpy(&finniganScan->Layout, &data[index + 12], 4);
  index += 40;

  // skip intensive parts if profile is false. still need to skip ahead though
  if (finniganScan->ProfileSize > 0) {
    memcpy(&finniganScan->FirstValue, &data[index], 8);
    memcpy(&finniganScan->Step, &data[index + 8], 8);
    memcpy(&finniganScan->PeakCount, &data[index + 16], 4);
    memcpy(&finniganScan->Nbins, &data[index + 20], 4);
    index += 24;
    finniganScan->chunks = (ProfileChunk*)malloc(sizeof(ProfileChunk) * finniganScan->PeakCount); // freed 299
    if (!profile) {
      for (int32_t i = 0; i < finniganScan->PeakCount; i++) {
        memcpy(&finniganScan->chunks[i].Nbins, &data[index + 4], 4);
        index += 8;
        if (finniganScan->Layout > 0) {
          index += 4;
        }
        index += 4 * finniganScan->chunks[i].Nbins;
      }
    } else {
      for (int32_t i = 0; i < finniganScan->PeakCount; i++) {
        memcpy(&finniganScan->chunks[i].Firstbin, &data[index], 4);
        memcpy(&finniganScan->chunks[i].Nbins, &data[index + 4], 4);
        index += 8;
        if (finniganScan->Layout > 0) {
          memcpy(&finniganScan->chunks[i].Fudge, &data[index], 4);
          index += 4;
        }
        // just point to the file for signals -- they'll need to be copied later anyways before we delete the data
        finniganScan->chunks[i].Signal = (float*)&data[index];
        index += 4 * finniganScan->chunks[i].Nbins;
      }
    }
  }
  // copying centroid is fast
  if (finniganScan->PeaklistSize > 0) {
    memcpy(&scan->centTotal, &data[index], 4);
    scan->centMzs = (float*)malloc(sizeof(float) * scan->centTotal);         // freed 456
    scan->centIntensities = (float*)malloc(sizeof(float) * scan->centTotal); // freed 456
    index += 4;
    for (int32_t i = 0; i < scan->centTotal; i++) {
      memcpy(&scan->centMzs[i], &data[index], 4);
      memcpy(&scan->centIntensities[i], &data[index + 4], 4);
      index += 8;
    }

    // get the flags
    finniganScan->flags = (unsigned char*)malloc(sizeof(unsigned char) * scan->centTotal);
    for (int32_t i = 0; i < scan->centTotal; i++) {
      index += 2;
      memcpy(&finniganScan->flags[i], &data[index], 1);
      index += 2;
    }
  }
}

int16_t fillScans(RawFile* rawFile, int profile, int centroid, FinniganInfo finniganInfo, FinniganScan* finniganScans) {
  Scan* scan;
  FinniganScan* finniganScan;
  uint64_t pos;

  double v;
  double tmpmz;
  int32_t k = 0; // tracks the individual bin in a
  for (uint32_t sn = 0; sn < rawFile->numScans; sn++) {
    uint32_t peakTotal = 0;
    scan = &rawFile->scans[sn];
    finniganScan = &finniganScans[sn];
    k = 0;
    readScanDataPacket(rawFile, finniganScan->offset, scan, profile, finniganInfo, finniganScan);

    // remove exception peaks
    for (uint32_t i = 0; i < finniganScan->PeakCount; i++) {
      if (finniganScan->flags[i] & (1 << 4)) continue;
      peakTotal += 1;
    }

    if (peakTotal != finniganScan->PeakCount) {
      float* fCentMzs = (float*)malloc(sizeof(float) * peakTotal);     // freed 456
      float* fIntensities = (float*)malloc(sizeof(float) * peakTotal); // freed 456
      uint32_t newIndex = 0;
      for (uint32_t i = 0; i < finniganScan->PeakCount; i++) {
        if (finniganScan->flags[i] & (1 << 4)) continue;

        fCentMzs[newIndex] = scan->centMzs[i];
        fIntensities[newIndex] = scan->centIntensities[i];
        newIndex += 1;
      }
      free(scan->centMzs);
      free(scan->centIntensities);
      scan->centMzs = fCentMzs;
      scan->centIntensities = fIntensities;
      scan->centTotal = peakTotal;
    }

    if (!profile) {
      free(finniganScan->chunks);
      free(finniganScan->flags);
      scan->prIntensities = NULL;
      scan->prMzs = NULL;
      continue;
    }

    // count number of bins in a scan (), and make it a flat []peak
    peakTotal = 0;
    for (int32_t pN = 0; pN < finniganScan->PeakCount; pN++) {
      if ((finniganScan->flags[pN] & (1 << 4))) continue;
      peakTotal += finniganScan->chunks[pN].Nbins;
    }
    if (peakTotal == 0) {
      peakTotal = scan->centTotal;
    }

    scan->prTotal = peakTotal;
    scan->prIntensities = (float*)malloc(sizeof(float) * peakTotal); // freed 454
    scan->prMzs = (float*)malloc(sizeof(float) * peakTotal);         // freed 455

    // convert the readings from Hz into mz/intensities
    for (int32_t pN = 0; pN < finniganScan->PeakCount; pN++) {
      if (finniganScan->flags[pN] & (1 << 4)) continue;
      for (int32_t binN = 0; binN < finniganScan->chunks[pN].Nbins; binN++) {
        // individual bins
        v = finniganScan->FirstValue + (finniganScan->chunks[pN].Firstbin + binN) * finniganScan->Step;
        switch (finniganScan->Nparam) {
        case 4:
          scan->prMzs[k] =
            finniganScan->A + finniganScan->B / v + finniganScan->C / v / v + finniganScan->chunks[pN].Fudge;
          break;
        case 5:
        case 7:
          scan->prMzs[k] = finniganScan->A + finniganScan->B / v / v + finniganScan->C / v / v / v / v +
            finniganScan->chunks[pN].Fudge;
          break;
        default:
          scan->prMzs[k] = v + finniganScan->chunks[pN].Fudge;
          break;
        }
        scan->prIntensities[k] = finniganScan->chunks[pN].Signal[binN];
        ++k;
      }
    }

    // put centroided data in mz/intensity lists if there's no profile data
    // if (scan->PeakCount == 0) {
    //   for (uint32_t centroidN = 0; centroidN < scan->centTotal; centroidN++) {
    //     scan->prMzs[k] = scan->centMzs[centroidN];
    //     scan->prIntensities[k] = scan->centIntensities[centroidN];
    //     ++k;
    //   }
    // }
    free(finniganScan->chunks);
    free(finniganScan->flags);
  }
}

int16_t initializeRawfile(RawFile* rawFile, int profile, int centroid) {
  // Skip to the version number in FileHeader
  // printf("In!\n");
  uint64_t pos = 36;
  uint32_t version;
  FinniganInfo finniganInfo;
  memcpy(&version, &rawFile->data[pos], 4);
  memcpy(&rawFile->timestamp, &rawFile->data[pos] + 4, 8);
  if (version < 66) {
    errorString = "Versions under 66 are not supported.";
    return -1;
  }
  pos = 1420;                                         // start right after InjectionData in SequencerRow
  pos = movePascal(rawFile->data, pos, 11);           // skip labels, comments, method files
  rawFile->fileName = fromPascal(rawFile->data, pos); // freed 471
  // printf("%s ", rawFile->fileName);
  pos = movePascal(rawFile->data, pos, 5); // skip filename, path, vial name
  pos += 4;
  pos = movePascal(rawFile->data, pos, 15);
  pos += 24; // AutoSamplerInfo
  pos = movePascal(rawFile->data, pos, 1);
  pos += 28; // skip the methodfile and timestamp
  memcpy(&rawFile->nControllers, &rawFile->data[pos], 4);
  pos += 780; // Skip other stuff in the preamble and the padding
  memcpy(&finniganInfo.dataAddr, &rawFile->data[pos], 8);
  pos += 16;                                                                          // dataaddr, unknown6
  int64_t* runHeaderAddr = (int64_t*)malloc(sizeof(int64_t) * rawFile->nControllers); // freed 340
  for (int16_t i = 0; i < rawFile->nControllers; i++) {
    memcpy(&runHeaderAddr[i], &rawFile->data[pos], 8);
    pos += 16; // skip runHeaderAddr and the unknown 7
  }
  finniganInfo.scanTrailerAddr = 0;
  for (int16_t i = 0; i < rawFile->nControllers && finniganInfo.scanTrailerAddr == 0; i++) {
    pos = runHeaderAddr[i];
    memcpy(&rawFile->FirstScanNumber, &rawFile->data[pos + 8], 4);
    memcpy(&rawFile->LastScanNumber, &rawFile->data[pos + 12], 4);
    // memcpy(&rawFile->lowMz, &rawFile->data[pos + 56], 8);
    // memcpy(&rawFile->highMz, &rawFile->data[pos + 64], 8);
    memcpy(&rawFile->startTime, &rawFile->data[pos + 72], 8);
    memcpy(&rawFile->endTime, &rawFile->data[pos + 80], 8);
    memcpy(&finniganInfo.ScanindexAddr, &rawFile->data[pos + 7408], 8);
    memcpy(&finniganInfo.scanTrailerAddr, &rawFile->data[pos + 7448], 8);
    memcpy(&finniganInfo.ScanparamsAddr, &rawFile->data[pos + 7456], 8);
    pos += 7588;
    pos = movePascal(rawFile->data, pos, 1);
    rawFile->instrumentModel = fromPascal(rawFile->data, pos); // freed 463
  }
  free(runHeaderAddr);
  // pull the scan Events
  rawFile->numScans = rawFile->LastScanNumber - rawFile->FirstScanNumber + 1;

  // detect incomplete files and exit instead of crashing later
  if (rawFile->nControllers == 0 || rawFile->numScans == 0) {
    // clean up after ourselves, since FreeRawFile won't work on a partially initialized file
    free(rawFile->fileName);
    if (rawFile->nControllers > 0) {
      free(rawFile->instrumentModel);
    }
    errorString = "Invalid raw file.";
    return -1;
  }
  Scan* scans = (Scan*)calloc(rawFile->numScans, sizeof(Scan)); // freed 471
  FinniganScan* finniganScans = (FinniganScan*)calloc(rawFile->numScans, sizeof(FinniganScan));
  rawFile->scans = scans;
  fillScanEvents(rawFile, finniganInfo, finniganScans);
  fillScanIndices(rawFile, finniganInfo, finniganScans);
  for (uint32_t i = 0; i < rawFile->numScans; i++) {
    finniganScans[i].offset += finniganInfo.dataAddr;
  }
  fillScans(rawFile, profile, centroid, finniganInfo, finniganScans);
  free(finniganScans);

  return 0;
}
