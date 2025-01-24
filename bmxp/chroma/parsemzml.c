#include "base64.h"
#include "numpress.c"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

typedef struct XmlElement {
  int hasClosing;
  int isClosing;
  char* name;
  char* attributes;
} XmlElement;

void getNextElement(XmlElement* element, char* fileArray, int64_t* pos) {
  // Gets the element and advances the pointer to the closing tag
  free(element->attributes);
  free(element->name);
  element->attributes = NULL;
  element->name = NULL;
  size_t result;
  char currentChar;
  char* rtnString;
  element->hasClosing = 0;
  element->isClosing = 0;
  int64_t startPos;
  int inSQuote = 0;
  int inDQuote = 0;
  // find the opening tag
  while (1) {
    currentChar = fileArray[*pos];
    if (currentChar == '<') {
      break;
    }
    *pos += 1;
  }

  ///////
  // find the name
  ///////

  // skip whitespace between opening symbol and name
  while (1) {
    *pos += 1;
    currentChar = fileArray[*pos];
    if (currentChar == ' ' || currentChar == '.' || currentChar == '\t' || currentChar == '\n') {
      continue;
    }
    startPos = *pos;
    break;
  }

  // if we see a slash before the name, it's a closing tag. label and strip whitespace
  if (currentChar == '/') {
    element->isClosing = 1;
    while (1) {
      *pos += 1;
      currentChar = fileArray[*pos];
      if (currentChar == ' ' || currentChar == '.' || currentChar == '\t' || currentChar == '\n') {
        continue;
      }
      startPos = *pos;
      break;
    }
  }

  // start the name
  while (1) {
    *pos += 1;
    currentChar = fileArray[*pos];

    // next whitepsace, slash or closing is the end of the name
    if (currentChar == ' ' || currentChar == '.' || currentChar == '\t' || currentChar == '\n' || currentChar == '>' ||
        currentChar == '/') {
      rtnString = malloc(*pos - startPos + 1); // has enough for the null terminator
      memcpy(rtnString, &fileArray[startPos], *pos - startPos);

      //   result = fread(rtnString, 1, finalPos - startPos, f);
      rtnString[*pos - startPos] = '\0';
      element->name = rtnString;
      break;
    }
  }

  // check the last character for whitespace
  while (1) {
    if (currentChar == '>') {
      // pos += 1;
      return;
    }
    if (currentChar == '/') element->hasClosing = 1;
    if (currentChar == ' ' || currentChar == '.' || currentChar == '\t' || currentChar == '\n' || currentChar == '/') {
      *pos += 1;
      currentChar = fileArray[*pos];
    }

    else
      break;
  }
  startPos = *pos;
  // pos and current are pointing to first character
  // check the next one

  while (1) {
    // result = fread(&currentChar, 1, 1, f);
    *pos += 1;
    currentChar = fileArray[*pos];

    // check if its quoted
    if (inSQuote && currentChar != '\'') continue;
    if (inDQuote && currentChar != '"') continue;
    if (currentChar == '\'') {
      inSQuote = !inSQuote;
      continue;
    }
    if (currentChar == '"') {
      inDQuote = !inDQuote;
      continue;
    }
    if (currentChar == '/') element->hasClosing = 1;

    // continue until we find the closing tag
    if (currentChar == '>') {
      rtnString = malloc(*pos - startPos + 1);
      memcpy(rtnString, &fileArray[startPos], *pos - startPos);
      rtnString[*pos - startPos] = '\0';
      element->attributes = rtnString;
      break;
    }
  }

  return;
}

char* getAttribute(char* attrString, char* attribute) {
  char* value = NULL;
  uint64_t pos = 0;
  int inSQuote = 0;
  int inDQuote = 0;

  char currentChar;
  uint32_t len = strlen(attribute);
  uint32_t attrStrLen = strlen(attrString);

  // find the attribute location, pos will be the starting position
  while (1) {
    // exit if there isn't enough space for =""
    if (len + pos + 3 >= attrStrLen) return NULL;
    int found = 1;
    currentChar = attrString[pos];

    // check if its quoted
    if (inSQuote && currentChar != '\'') {
      pos++;
      continue;
    };
    if (inDQuote && currentChar != '"') {
      pos++;
      continue;
    };
    if (currentChar == '\'') {
      pos++;
      inSQuote = !inSQuote;
      continue;
    }
    if (currentChar == '"') {
      pos++;
      inDQuote = !inDQuote;
      continue;
    }

    // set to not found if there is a missing equality
    for (uint32_t i = 0; i < len; i++) {
      if (attrString[pos + i] != attribute[i]) {
        found = 0;
        break;
      }
    }
    if (!found) {
      pos++;
      continue;
    }
    // if found, check that the terminator matches up with = or  ' '. otherwise it's a false hit
    if (attrString[pos + len] == ' ' || attrString[pos + len] == '.' || attrString[pos + len] == '\t' ||
        attrString[pos + len] == '\n' || attrString[pos + len] == '=') {
    } else {
      found = 0;
    }
    // iterate and continue

    if (!found) {
      pos++;
      continue;
    }

    pos = pos + len;

    // find first double or single quote
    currentChar = attrString[pos];
    while (currentChar != '"' && currentChar != '\'') {
      pos++;
      if (pos >= attrStrLen) return NULL;
      currentChar = attrString[pos];
    }
    char quoteChar = currentChar;
    uint32_t start = pos;
    pos++;
    if (pos >= attrStrLen) return NULL;
    currentChar = attrString[pos];
    while (currentChar != quoteChar) {
      pos++;
      if (pos >= attrStrLen) return NULL;
      currentChar = attrString[pos];
    }
    char* rtnString = malloc(pos - start);
    memcpy(rtnString, &attrString[start + 1], pos - start);
    rtnString[pos - start - 1] = '\0';
    return rtnString;
  }
}

void mzml_handleFileDescription(RawFile* rawFile, XmlElement* element, uint64_t* pos) {
  while (1) {
    getNextElement(element, rawFile->data, pos);
    // maybe do something here eventually
    if (strcmp(element->name, "fileDescription") == 0 && element->isClosing) {
      return;
    }
  }
}

void mzml_handleReferenceableGroup(RawFile* rawFile, XmlElement* element, uint64_t* pos) {
  rawFile->instrumentModel = "Unknown";

  while (1) {
    getNextElement(element, rawFile->data, pos);
    if (strcmp(element->name, "cvParam") == 0) {
      char* accession = getAttribute(element->attributes, "accession");
      free(accession);
    }
    // maybe do something here eventually
    if (strcmp(element->name, "referenceableParamGroupList") == 0 && element->isClosing) {
      return;
    }
  }
}

void mzml_fillSoftware(RawFile* rawFile, XmlElement* element, uint64_t* pos) {
  while (1) {
    getNextElement(element, rawFile->data, pos);
    // maybe do something here eventually
    if (strcmp(element->name, "softwareList") == 0 && element->isClosing) {
      return;
    }
  }
}

void mzml_handleInstrumentConfig(RawFile* rawFile, XmlElement* element, uint64_t* pos) {
  while (1) {
    getNextElement(element, rawFile->data, pos);
    // maybe do something here eventually
    if (strcmp(element->name, "instrumentConfigurationList") == 0 && element->isClosing) {
      return;
    }
  }
}

void mzml_handleDataProcessing(RawFile* rawFile, XmlElement* element, uint64_t* pos) {
  while (1) {
    getNextElement(element, rawFile->data, pos);
    // maybe do something here eventually
    if (strcmp(element->name, "dataProcessingList") == 0 && element->isClosing) {
      return;
    }
  }
}

char* getData(char* fileArray, int64_t* pos) {
  // gets data and advances to the last character
  char* data;
  *pos += 1; // go to the actual start
  uint64_t start = *pos;
  char currentChar = fileArray[*pos];
  while (currentChar != '<') {
    *pos += 1;
    currentChar = fileArray[*pos];
  }
  data = malloc(*pos - start + 1);
  memcpy(data, &fileArray[start], *pos - start);
  data[*pos - start] = '\0';
  return data;
}

float* decodeBinary(char* b64data, char* compression, char* units, uint32_t arrLength) {
  float* values; // final result
  uLongf decLen;
  size_t b64DecLen;
  unsigned char* b64Dec = NULL;
  unsigned char* decompressed = NULL;

  // everything needs to be base64 decoded
  b64DecLen = b64_decoded_size(b64data);
  b64Dec = malloc(b64DecLen);
  b64_decode(b64data, b64Dec, b64DecLen);

  // uncompressed, no zlib
  if (strcmp(compression, "MS:1000576") == 0 || strcmp(compression, "MS:1002312") == 0 ||
      strcmp(compression, "MS:1002313") == 0 || strcmp(compression, "MS:1002314") == 0) {
    decompressed = b64Dec;
    decLen = b64DecLen;
    b64Dec = NULL;
    // printf("Decoded size: %d %d ", b64DecLen, arrLength);

  }
  // zlib
  else {
    decLen = 8 * arrLength; // maximum theoretical
    decompressed = (unsigned char*)malloc(decLen);
    int res = uncompress(decompressed, &decLen, b64Dec, b64DecLen);
  }

  // if no numpress, cast to correct size
  if (strcmp(compression, "MS:1000574") == 0 || strcmp(compression, "MS:1000576") == 0) {
    if (strcmp(units, "MS:1000521") == 0) {
      values = malloc(arrLength * 4);
      memcpy(values, decompressed, arrLength * 4);
    } else if (strcmp(units, "MS:1000523") == 0) {
      double* doubleValues = (double*)decompressed;
      values = malloc(arrLength * 4);
      for (uint32_t i = 0; i < arrLength; i++) {
        values[i] = doubleValues[i];
      }
    }
  }
  // if it's numpress, uncompress and cast as 32bit
  else {
    DoubleResult npResult;
    if (strcmp(compression, "MS:1002312") == 0 || strcmp(compression, "MS:1002746") == 0) // numpress linear
      npResult = decodeLinear(decompressed, decLen);
    else if (strcmp(compression, "MS:1002313"), 0 || strcmp(compression, "MS:1002747") == 0) // needs positive int
      npResult = decodePic(decompressed, decLen);
    else if (strcmp(compression, "MS:1002314"), 0 || strcmp(compression, "MS:1002748") == 0) // short logged float
      npResult = decodeSlof(decompressed, decLen);
    values = malloc(npResult.len * 4);
    for (uint32_t i = 0; i < npResult.len; i++) {
      values[i] = (float)npResult.data[i];
    }
    free(npResult.data);
  }
  free(b64Dec);
  free(decompressed);
  return values;
}

void mzml_handleScanCvParams(RawFile* rawFile, uint64_t* pos, XmlElement* element, Scan* scan, ScanFilter* scanFilter,
                             int* numCEs) {

  char* accession = getAttribute(element->attributes, "accession");
  // found in spectrumlist
  if (strcmp(accession, "MS:1000511") == 0) { //  ms level
    char* tempStr = getAttribute(element->attributes, "value");
    scanFilter->msLevel = atoi(tempStr);
    free(tempStr);
  } else if (strcmp(accession, "MS:1000127") == 0) { // centroid
    scanFilter->scanMode = 0;
  } else if (strcmp(accession, "MS:1000128") == 0) { // profile
    scanFilter->scanMode = 1;
  } else if (strcmp(accession, "MS:1000129") == 0) { // negative polarity
    scanFilter->polarity = 0;
  } else if (strcmp(accession, "MS:1000130") == 0) { // positive polarity
    scanFilter->polarity = 1;
  }
  // found in scanList
  else if (strcmp(accession, "MS:1000501") == 0) { //  low mz
    char* tempStr = getAttribute(element->attributes, "value");
    scanFilter->lowMass = atof(tempStr);
    free(tempStr);
  } else if (strcmp(accession, "MS:1000500") == 0) { //  high mz
    char* tempStr = getAttribute(element->attributes, "value");
    scanFilter->highMass = atof(tempStr);
    free(tempStr);
  } else if (strcmp(accession, "MS:1000016") == 0) { // retention time
    char* tempStr = getAttribute(element->attributes, "value");
    scan->time = atof(tempStr);
    free(tempStr);
  } else if (strcmp(accession, "MS:1000827") == 0) { // isolation window
    scanFilter->nPrecursors = scanFilter->nPrecursors + 1;
    scanFilter->precursorMz = realloc(scanFilter->precursorMz, scanFilter->nPrecursors * 8);
    char* tempStr = getAttribute(element->attributes, "value");
    scanFilter->precursorMz[scanFilter->nPrecursors - 1] = atof(tempStr);
    free(tempStr);
  } else if (strcmp(accession, "MS:1000045") == 0) { // collision energy
    *numCEs = *numCEs + 1;
    scanFilter->energy = realloc(scanFilter->energy, *numCEs * 8);
    char* tempStr = getAttribute(element->attributes, "value");
    scan->time = atof(tempStr);
    scanFilter->energy[*numCEs - 1] = atof(tempStr);
    free(tempStr);
  }
  free(accession);
}

void mzml_handleBinaryData(RawFile* rawFile, uint64_t* pos, XmlElement* element, Scan* scan, ScanFilter* scanFilter,
                           uint32_t arrLength) {
  // Binary Data for Scans
  char* compression = NULL;
  char* units = NULL;
  int numpress = 0;
  char* arrayType = NULL;
  unsigned char* b64Data;
  while (!(strcmp(element->name, "binaryDataArray") == 0 && element->isClosing)) {
    getNextElement(element, rawFile->data, pos);
    if (element->isClosing) {
      continue;
    }

    if (strcmp(element->name, "cvParam") == 0) {
      char* accession = getAttribute(element->attributes, "accession");

      // get the compression info and array type
      if (strcmp(accession, "MS:1000574") == 0 || strcmp(accession, "MS:1000576") == 0 ||
          strcmp(accession, "MS:1002312") == 0 || strcmp(accession, "MS:1002313") == 0 ||
          strcmp(accession, "MS:1002314") == 0 || strcmp(accession, "MS:1002746") == 0 ||
          strcmp(accession, "MS:1002747") == 0 || strcmp(accession, "MS:1002748") == 0) {
        free(compression);
        compression = accession;
      }
      // units, i.e. float, double, int54
      else if (strcmp(accession, "MS:1000521") == 0 || strcmp(accession, "MS:1000523") == 0 ||
               strcmp(accession, "MS:1000519") == 0 || strcmp(accession, "MS:1000522") == 0 ||
               strcmp(accession, "MS:1001479") == 0) {
        free(units);
        units = accession;
      }
      // type of data, i.e. 514 m/z, 515 intensity
      else if (strcmp(accession, "MS:1000515") == 0 || strcmp(accession, "MS:1000514") == 0) {
        free(arrayType);
        arrayType = accession;
      }
      // otherwise ignore
      else {
        free(accession);
      }
    }
    if (strcmp(element->name, "binary") == 0 && !element->isClosing) b64Data = getData(rawFile->data, pos);
  }

  // with everything retrieved, fill in the Raw data
  // printf(" %s %s %s %d %s\n", units, arrayType, compression, arrLength, b64Data);
  float* decoded = decodeBinary(b64Data, compression, units, arrLength);
  if (scanFilter->scanMode == 0) {
    if (strcmp(arrayType, "MS:1000515") == 0) { // intensity
      scan->centIntensities = decoded;
    } else if (strcmp(arrayType, "MS:1000514") == 0) { // mz
      scan->centMzs = decoded;
    }
    scan->centTotal = arrLength;
  } else if (scanFilter->scanMode == 1) {
    if (strcmp(arrayType, "MS:1000515") == 0) // intensity
    {
      scan->prIntensities = decoded;
    } else if (strcmp(arrayType, "MS:1000514") == 0) // mz
      scan->prMzs = decoded;

    scan->prTotal = arrLength;
  }
  free(b64Data);
  free(arrayType);
  free(units);
  free(compression);
}

void mzml_handleChromBinaryData(RawFile* rawFile, uint64_t* pos, XmlElement* element, Chromatogram* chrom,
                                uint32_t arrLength) {
  // Binary Data for Chromatograms
  char* compression = NULL;
  char* units = NULL;
  int numpress = 0;
  char* arrayType = NULL;
  unsigned char* b64Data = NULL;

  // as long as it's not the closing tag for 'binaryDataArry'
  while (!(strcmp(element->name, "binaryDataArray") == 0 && element->isClosing)) {
    getNextElement(element, rawFile->data, pos);
    if (element->isClosing) {
      continue;
    }

    // get the compression info and array type
    if (strcmp(element->name, "cvParam") == 0) {
      char* accession = getAttribute(element->attributes, "accession");
      if (strcmp(accession, "MS:1000574") == 0 || strcmp(accession, "MS:1000576") == 0 ||
          strcmp(accession, "MS:1002312") == 0 || strcmp(accession, "MS:1002313") == 0 ||
          strcmp(accession, "MS:1002314") == 0 || strcmp(accession, "MS:1002746") == 0 ||
          strcmp(accession, "MS:1002747") == 0 || strcmp(accession, "MS:1002748") == 0) {
        free(compression);
        compression = accession;
      }
      // unit type -- float32, etc...
      else if (strcmp(accession, "MS:1000521") == 0 || strcmp(accession, "MS:1000523") == 0 ||
               strcmp(accession, "MS:1000519") == 0 || strcmp(accession, "MS:1000522") == 0 ||
               strcmp(accession, "MS:1001479") == 0) {
        free(units);
        units = accession;
      }
      // data type -- rts, intensities, etc..
      else if (strcmp(accession, "MS:1000515") == 0 || strcmp(accession, "MS:1000514") == 0 ||
               strcmp(accession, "MS:1000595") == 0) {
        free(arrayType);
        arrayType = accession;
      } else {
        free(accession);
      }
    }
    if (strcmp(element->name, "binary") == 0 && !element->isClosing) b64Data = getData(rawFile->data, pos);
  }
  // skip if not intensity or rt
  if (arrayType == NULL) {
    free(b64Data);
    free(arrayType);
    free(units);
    free(compression);
    return;
  }
  // printf("decoding... %d %s %s %s %s", arrLength, arrayType, compression, units, b64Data);
  // with everything retrieved, fill in the Raw data
  float* decoded = decodeBinary(b64Data, compression, units, arrLength);

  if (strcmp(arrayType, "MS:1000515") == 0) { // intensity
    chrom->intensities = decoded;
  } else if (strcmp(arrayType, "MS:1000595") == 0) { // rts
    chrom->rts = decoded;
  }
  free(b64Data);
  free(arrayType);
  free(units);
  free(compression);
}

float mzml_getIsolation(RawFile* rawFile, uint64_t* pos, XmlElement* element) {
  float results = 0;
  char* openingTag = malloc(strlen(element->name) + 1);
  strcpy(openingTag, element->name);
  while (!(strcmp(element->name, openingTag) == 0 && element->isClosing)) {
    getNextElement(element, rawFile->data, pos);
    if (element->isClosing) {
      continue;
    }
    if (strcmp(element->name, "isolationWindow") == 0) {
      while (!(strcmp(element->name, "isolationWindow") == 0 && element->isClosing)) {
        getNextElement(element, rawFile->data, pos);
        if (strcmp(element->name, "cvParam") == 0) {
          char* accession = getAttribute(element->attributes, "accession");
          if (strcmp(accession, "MS:1000827") == 0) {
            char* floatStr = getAttribute(element->attributes, "value");
            results = atof(floatStr);
            free(floatStr);
          }
          free(accession);
        }
      }
    }
  }
  free(openingTag);
  return results;
}

void mzml_createSpectrum(RawFile* rawFile, XmlElement* element, uint64_t* pos) {
  // iterate through the spectrumlist elements
  uint32_t scanFilterIndex = 0;

  while (!(strcmp(element->name, "spectrumList") == 0 && element->isClosing)) {
    getNextElement(element, rawFile->data, pos);
    if (strcmp(element->name, "spectrum") == 0 && !element->hasClosing) {
      // get scan index
      char* scanIndexStr = getAttribute(element->attributes, "index");
      uint32_t scanIndex = atoi(scanIndexStr);
      free(scanIndexStr);
      Scan* scan = &rawFile->scans[scanIndex];
      ScanFilter scanFilter;
      scanFilter.precursorMz = NULL;
      scanFilter.energy = NULL;

      // length of rts/mzs
      char* arrayLengthStr = getAttribute(element->attributes, "defaultArrayLength");
      uint32_t arrLength = atoi(arrayLengthStr);
      free(arrayLengthStr);

      // clunkily find the scan number
      char* idStr = getAttribute(element->attributes, "id");
      char* start = strstr(idStr, "scan=");
      uint32_t charPos = 5;
      for (charPos; charPos < strlen(start); charPos++) {
        if (start[charPos] == ' ' || start[charPos] == '"' || start[charPos] == '\'') break;
      }
      start[charPos] = '\0';
      uint32_t scanNumber = atoi(&start[5]);
      free(idStr);

      // defaults
      int* numCEs = malloc(sizeof(int));
      *numCEs = 0;

      scanFilter.msLevel = 1;             // no ms2
      scanFilter.nPrecursors = 0;         // no precursors
      scanFilter.precursorMz = malloc(0); // no ms2
      scanFilter.energy = malloc(0);      // no ms2

      // defaults
      scanFilter.polarity = 2;
      scanFilter.scanType = 0; // full
      scanFilter.scanMode = 0; // centroid
      scanFilter.analyzer = 6;
      // scanFilter.dependentScans = 0;
      scanFilter.lowMass = 0;
      scanFilter.highMass = 0;

      // iterate through the spectrum stuff
      while (!(strcmp(element->name, "spectrum") == 0 && element->isClosing)) {
        getNextElement(element, rawFile->data, pos);
        if (strcmp(element->name, "cvParam") == 0)
          mzml_handleScanCvParams(rawFile, pos, element, scan, &scanFilter, numCEs);
        else if (strcmp(element->name, "binaryDataArray") == 0) {
          mzml_handleBinaryData(rawFile, pos, element, scan, &scanFilter, arrLength);
        }
      }

      // check and add scan filter
      int found = 0;
      for (uint32_t k = 0; k < scanFilterIndex; k++) {
        if (compareFilters(&scanFilter, &rawFile->scanFilters[k])) {
          // already exists, clear and
          scan->filter = k;
          free(scanFilter.precursorMz);
          free(scanFilter.energy);
          found = 1;
          break;
        }
      }
      if (!found) {
        // add new filter, don't clear the arrays
        // printf("Adding\n");
        scan->filter = scanFilterIndex;
        rawFile->scanFilters[scanFilterIndex] = scanFilter;
        scanFilterIndex++;
      }
      free(numCEs);
    }
  }

  rawFile->nFilters = scanFilterIndex;
}

void mzml_createChromatograms(RawFile* rawFile, XmlElement* element, uint64_t* pos) {
  // iterate through the chromatogramList elements
  while (!(strcmp(element->name, "chromatogramList") == 0 && element->isClosing)) {
    getNextElement(element, rawFile->data, pos);
    if (strcmp(element->name, "chromatogram") == 0 && !element->hasClosing) {
      // get scan index
      char* chromIndexStr = getAttribute(element->attributes, "index");
      uint32_t chromIndex = atoi(chromIndexStr);
      free(chromIndexStr);
      Chromatogram* chrom = &rawFile->chromatograms[chromIndex];

      // get scan id
      char* chromId = getAttribute(element->attributes, "id");
      chrom->id = chromId;

      // length of rts/mzs
      char* arrayLengthStr = getAttribute(element->attributes, "defaultArrayLength");
      uint32_t arrLength = atoi(arrayLengthStr);
      free(arrayLengthStr);
      chrom->length = arrLength;

      // find relevant info in the chromatogram section
      while (!(strcmp(element->name, "chromatogram") == 0 && element->isClosing)) {
        getNextElement(element, rawFile->data, pos);
        if (strcmp(element->name, "precursor") == 0) {
          chrom->precursorMz = mzml_getIsolation(rawFile, pos, element);
        }
        if (strcmp(element->name, "product") == 0) {
          chrom->productMz = mzml_getIsolation(rawFile, pos, element);
        }
        if (strcmp(element->name, "binaryDataArray") == 0)
          mzml_handleChromBinaryData(rawFile, pos, element, chrom, arrLength);
      }
    }
  }
}

void mzml_fillScans(RawFile* rawFile, XmlElement* element, uint64_t* pos) {
  while (1) {
    getNextElement(element, rawFile->data, pos);

    // create scans if it has a spectrumList
    if (strcmp(element->name, "spectrumList") == 0 && !element->isClosing) {
      char* value = getAttribute(element->attributes, "count");
      uint32_t numScans = atoi(value);
      free(value);

      rawFile->numScans = numScans;
      rawFile->nFilters = numScans; // just allocate scans, no comparison for right now
      rawFile->scans = calloc(numScans, sizeof(Scan));
      rawFile->scanFilters = malloc(sizeof(ScanFilter) * numScans);

      for (uint32_t i = 0; i < numScans; i++) {
        rawFile->scanFilters[i].precursorMz = NULL; // no ms2
        rawFile->scanFilters[i].energy = NULL;      // no ms2
        rawFile->scans[i].prIntensities = NULL;
        rawFile->scans[i].prMzs = NULL;
        rawFile->scans[i].centIntensities = NULL;
        rawFile->scans[i].centMzs = NULL;
        rawFile->scans[i].prTotal = 0;
        rawFile->scans[i].centTotal = 0;
      }

      mzml_createSpectrum(rawFile, element, pos);
    }

    // create chromatograms if it has a chromatogramList
    if (strcmp(element->name, "chromatogramList") == 0 && !element->isClosing) {
      char* value = getAttribute(element->attributes, "count");
      uint32_t numChroms = atoi(value);
      free(value);
      rawFile->numChroms = numChroms;
      rawFile->chromatograms = calloc(numChroms, sizeof(Chromatogram));
      for (uint32_t i = 0; i < numChroms; i++) {
        rawFile->chromatograms[i].rts = NULL;
        rawFile->chromatograms[i].intensities = NULL;
        rawFile->chromatograms[i].id = NULL;
        rawFile->chromatograms[i].length = 0;
        rawFile->chromatograms[i].precursorMz = 0;
        rawFile->chromatograms[i].productMz = 0;
      }
      mzml_createChromatograms(rawFile, element, pos);
    }

    // escape when we see the run tag
    if (strcmp(element->name, "run") == 0 && element->isClosing) {
      return;
    }
  }
}

int16_t initializeMzml(RawFile* rawFile, int profile, int centroid) {
  // printf("Initializing mzml...");
  XmlElement* element = malloc(sizeof(XmlElement));
  element->name = NULL;
  element->attributes = NULL;

  int64_t pos = 0;
  while (1) {
    getNextElement(element, rawFile->data, &pos);
    if (strcmp(element->name, "fileDescription") == 0 && (!element->hasClosing))
      mzml_handleFileDescription(rawFile, element, &pos);
    else if (strcmp(element->name, "referenceableParamGroupList") == 0 && (!element->hasClosing))
      mzml_handleReferenceableGroup(rawFile, element, &pos);
    else if (strcmp(element->name, "softwareList") == 0 && (!element->hasClosing))
      mzml_fillSoftware(rawFile, element, &pos);
    else if (strcmp(element->name, "instrumentConfigurationList") == 0 && (!element->hasClosing))
      mzml_handleInstrumentConfig(rawFile, element, &pos);
    else if (strcmp(element->name, "dataProcessingList") == 0 && (!element->hasClosing))
      mzml_handleDataProcessing(rawFile, element, &pos);
    else if (strcmp(element->name, "run") == 0 && (!element->hasClosing))
      mzml_fillScans(rawFile, element, &pos);

    if (strcmp(element->name, "mzML") == 0 && element->isClosing) break;
  }
  // fill these in...
  rawFile->fileName = malloc(9);
  strcpy(rawFile->fileName, "Filename");
  rawFile->instrumentModel = malloc(6);
  strcpy(rawFile->instrumentModel, "Model");

  free(element->name);
  free(element->attributes);
  free(element);
  return 1;
}
//
//