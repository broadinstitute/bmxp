# Chroma

## Quick use
1. Install with `pip install bmxp` (requires Python 3)
2. Get a rawfile or mzml file
3. Run this code:

```python
from bmxp.chroma import RawFile
rf = RawFile("test.raw")
xic = rf.xic(90.0550, 7.5, 8.5)
print(xic.rt)
print(xic.intensity)
print(xic.mz)
```

This will print the xic for m/z = 90.0550 from 7.5 to 8.5 minutes. The default ppm threshold is 20 ppm.

Chroma is a rapid parser for rawfiles and mzml files. It is designed for speed to rapidly parse data in an acquisition file into a usable format, where XICs and scan information can be extracted. For this, it is not the best tool for deep exploration of metadata.

Chroma is experimental, but currently it:
* Reads Thermo Rawfiles
* Reads .mzml, supporting numpress and zlib
* Reads spectrum data (orbi, tof)
* Reads chromatogram data (QQQ)
* Parses xics
* Parses scans and scan filters
* Supports MS2
* Runs on Windows and Linux (we can help compile for other platforms)

There are many improvements yet to make.
