# pylint: disable=too-many-instance-attributes,missing-class-docstring
"""
Package to open and analyze raw ms files
"""
import ctypes
from ctypes import (
    CDLL,
    c_void_p,
    c_char_p,
    c_float,
    c_double,
    c_uint8,
    c_int16,
    c_int32,
    c_bool,
    c_uint32,
    c_uint64,
    Structure,
    POINTER,
    cast,
)
import os
import platform
from enum import Enum
import numpy as np

__version__ = "0.0.3"


dir_path = os.path.dirname(os.path.realpath(__file__))
if platform.system() == "Windows":
    rawfilereader = CDLL(os.path.join(dir_path, "rawfilereader.dll"))
else:
    rawfilereader = CDLL(os.path.join(dir_path, "rawfilereader.so"))


class Empty:
    pass


class Analyzer(Enum):
    ITMS = 0
    TQMS = 1
    SQMS = 2
    TQFMS = 3
    FTMS = 4
    SECTOR = 5
    UNDEFINED = 6


class ScanMode(Enum):
    CENTROID = 0
    PROFILE = 1
    UNDEFINED = 2


class ScanType(Enum):
    FULL = 0
    ZOOM = 1
    SIM = 2
    SRM = 3
    CRM = 4
    UNDEFINED = 5
    Q1 = 6
    Q3 = 7


class Polarity(Enum):
    NEGATIVE = 0
    POSITIVE = 1
    UNDEFINED = 2


class CXic(Structure):  # pylint: disable=too-few-public-methods
    """
    C Struct for the XIC
    """

    _fields_ = [
        ("rt", c_void_p),
        ("intensity", c_void_p),
        ("length", c_int32),
        ("mzs", c_void_p),
    ]


class CScanFilter(Structure):  # pylint: disable=too-few-public-methods
    """
    C Struct for the ScanFilter
    """

    _fields_ = [
        ("polarity", c_uint8),
        ("analyzer", c_uint8),
        ("ms_level", c_uint8),
        ("n_precursors", c_uint32),
        ("precursor_mz", c_void_p),
        ("energy", c_void_p),
        ("low_mass", c_double),
        ("high_mass", c_double),
        ("scan_mode", c_uint8),
        ("scan_type", c_uint8),
        # ("dependent_scans", c_uint8),
    ]


class CRawFile(Structure):  # pylint: disable=too-few-public-methods
    """
    C Struct for the RawFile
    """

    _fields_ = [
        ("instrument_model", c_char_p),
        ("filename", c_char_p),
        ("n_controllers", c_uint32),
        ("scans", c_void_p),
        ("data", c_void_p),  # pointer to full data, instead of raw bytes
        ("num_scans", c_uint32),
        ("first_scan_number", c_uint32),
        ("last_scan_number", c_uint32),
        ("start_time", c_double),
        ("end_time", c_double),
        ("scan_filters", c_void_p),
        ("n_filters", c_uint32),
        ("timestamp", c_uint64),
        ("chromatograms", c_void_p),
        ("num_chroms", c_uint32),
        # ("scan_trailer_addr", c_uint64),
        # ("low_mz", c_double),
        # ("high_mz", c_double),
        # ("scan_index_addr", c_uint64),
        # ("scan_params_addr", c_uint64),
        # ("version", c_uint32),
        # ("data_addr", c_uint64),
    ]


class CScan(Structure):  # pylint: disable=too-few-public-methods
    """
    C Struct for Scan
    """

    _fields_ = [
        ("prIntensities", c_void_p),
        ("prMzs", c_void_p),
        ("centIntensities", c_void_p),
        ("centMzs", c_void_p),
        ("prTotal", c_uint32),
        ("centTotal", c_uint32),
        ("time", c_double),
        ("index", c_uint32),
        ("filter", c_uint32),
    ]


class CChrom(Structure):  # pylint: disable=too-few-public-methods
    """
    C Struct for Chromatogram
    """

    _fields_ = [
        ("id", c_void_p),
        ("intensities", c_void_p),
        ("rts", c_void_p),
        ("length", c_uint32),
        ("precursorMz", c_float),
        ("productMz", c_float),
    ]


# Define the arg and return types of the function
rawfilereader.Open.argtypes = [c_char_p, c_int32, c_int32, c_char_p, c_uint64]
rawfilereader.Open.restype = c_void_p  # RawFile, just store as a pointer

rawfilereader.FromBytes.argtypes = [c_char_p, c_int32, c_int32, c_char_p]
rawfilereader.FromBytes.restype = c_void_p  # RawFile, just store as a pointer

rawfilereader.FreeRawFile.argtypes = [c_void_p]
rawfilereader.FreeRawFile.restype = c_int16

rawfilereader.Pull_xic.argtypes = [
    c_void_p,
    c_float,
    c_float,
    c_float,
    c_float,
    c_int32,
    c_bool,
    c_bool,
]
rawfilereader.Pull_xic.restype = CXic

rawfilereader.FreeP.argtypes = [c_void_p]

rawfilereader.GetErrorString.restype = c_char_p


# Define our usable structures
class Xic:  # pylint: disable=too-few-public-methods
    """
    An Extracted Ion Chromatogram. Mainly, it stores numpy arrays for
    retention times and intensities
    """

    def __init__(self, c_xic=None):
        """Initialize from a C_Xic struct"""
        if c_xic is None:
            self.length = 0
            self.rt = np.array([])
            self.intensity = np.array([])
            self.mz = np.array([])
            return
        if isinstance(c_xic, dict):
            self.rt = c_xic["rt"]
            self.intensity = c_xic["intensity"]
            self.length = c_xic["length"]
            self.mz = np.array([c_xic["mz"]] * c_xic["length"])
            return

        rt_p = c_xic.rt
        intensity_p = c_xic.intensity
        mz_p = c_xic.mzs
        self.length = c_xic.length
        self.rt = np.copy(
            np.frombuffer((c_float * c_xic.length).from_address(rt_p), np.float32)
        )
        self.intensity = np.copy(
            np.frombuffer(
                (c_float * c_xic.length).from_address(intensity_p), np.float32
            )
        )
        if mz_p:
            self.mz = np.copy(
                np.frombuffer((c_float * c_xic.length).from_address(mz_p), np.float32)
            )
        else:
            self.mz = np.zeros(self.length, np.float32)

        rawfilereader.FreeP(rt_p)
        rawfilereader.FreeP(intensity_p)
        rawfilereader.FreeP(mz_p)


class ScanFilter:  # pylint: disable=too-few-public-methods
    """
    A scan filter, an assortment of properties that a scan can have
    """

    def __init__(self, scanf):
        self.polarity = Polarity(scanf["polarity"])
        self.scan_mode = ScanMode(scanf["scan_mode"])
        # self.dependent_scans = scanf["dependent_scans"]
        self.ms_level = scanf["ms_level"]
        self.scan_type = ScanType(scanf["scan_type"])
        self.analyzer = Analyzer(scanf["analyzer"])
        self.n_precursors = scanf["n_precursors"]
        self.low_mass = scanf["low_mass"]
        self.high_mass = scanf["high_mass"]

        # with such short lists, this seems to be faster than np.frombuffer
        self.precursor_mz = list(
            (c_double * self.n_precursors).from_address(int(scanf["precursor_mz"]))
        )
        self.energy = list(
            (c_double * self.n_precursors).from_address(int(scanf["energy"]))
        )

    def __str__(self):
        if self.ms_level == 1:
            return f"MS1; Range:{self.low_mass}-{self.high_mass}"
        return (
            f"MS{self.ms_level}; Precursors:{self.precursor_mz}; Energy:{self.energy};"
            f" Range:{self.low_mass}-{self.high_mass}"
        )

    def __repr__(self):
        return self.__str__()


class Scan:  # pylint: disable=too-few-public-methods
    """
    An individual scan
    """

    def __init__(self, c_scan, profile=False, centroid=True):
        self.profile = profile
        self.centroid = centroid
        self.time = c_scan["time"]
        self.index = c_scan["index"]
        self.cent_total = int(c_scan["centTotal"])
        self.pr_total = int(c_scan["prTotal"])
        self.scan_filter = c_scan["filter"]

        if self.centroid:
            self.cent_mzs = np.frombuffer(
                (c_float * c_scan["centTotal"]).from_address(int(c_scan["centMzs"])),
                c_float,
            )
            self.cent_intensities = np.frombuffer(
                (c_float * c_scan["centTotal"]).from_address(
                    int(c_scan["centIntensities"])
                ),
                c_float,
            )
        if self.profile:
            self.profile_mzs = np.frombuffer(
                (c_float * c_scan["prTotal"]).from_address(int(c_scan["prMzs"])),
                c_float,
            )
            self.profile_intensities = np.frombuffer(
                (c_float * c_scan["prTotal"]).from_address(
                    int(c_scan["prIntensities"])
                ),
                c_float,
            )


class Chromatogram:  # pylint: disable=too-few-public-methods
    """
    An individual scan
    """

    def __init__(self, c_chrom):
        self.id = ctypes.string_at(int(c_chrom["id"]))
        self.length = int(c_chrom["length"])
        self.precursor_mz = c_chrom["precursorMz"]
        self.product_mz = c_chrom["productMz"]

        self.rts = np.frombuffer(
            (c_float * c_chrom["length"]).from_address(int(c_chrom["rts"])),
            c_float,
        )
        self.intensities = np.frombuffer(
            (c_float * c_chrom["length"]).from_address(int(c_chrom["intensities"])),
            c_float,
        )


class RawFile:
    """
    The Open rawfile. Constructor accepts a filename string or file contents as bytes.
    """

    def __init__(self, file, profile=False, centroid=True, file_format=None):
        if isinstance(file, str):
            file_size = os.path.getsize(file)
            if file_format is None:
                if file.lower().endswith(".raw"):
                    file_format = "rawfile"
                else:
                    file_format = "mzml"
            format_char = file_format.encode("utf-8")
            filename_char = file.encode("utf-8")
            ptr = rawfilereader.Open(
                filename_char, profile, centroid, format_char, file_size
            )
        elif isinstance(file, bytes):
            if file_format is None:
                file_format = "rawfile"
            format_char = file_format.encode("utf-8")
            ptr = rawfilereader.FromBytes(file, profile, centroid, format_char)
        else:
            raise TypeError(
                "RawFile can be constructed from a string (filename) or bytes "
                f"(file contents). You provided {type(file)}."
            )
        if not ptr:  # null pointer if reading the raw file failed
            raise RuntimeError(rawfilereader.GetErrorString().decode("UTF-8"))
        rf = cast(ptr, POINTER(CRawFile)).contents  # dereference pointer

        self._data_p = ptr
        try:
            self.instrument_model = rf.instrument_model.decode("UTF-8")
        except:
            self.instrument_model = "unknown"
        try:
            self.filename = rf.filename.decode("UTF-8")
        except UnicodeDecodeError:
            self.filename = ""
        self.profile = profile
        self.centroid = centroid
        self.data = rf.data  # pointer
        self.num_scans = rf.num_scans
        self.num_chroms = rf.num_chroms

        # self.scan_trailer_addr = rf.scan_trailer_addr
        self.first_scan_number = rf.first_scan_number
        self.last_scan_number = rf.last_scan_number
        # self.low_mz = rf.low_mz
        # self.high_mz = rf.high_mz
        self.start_time = rf.start_time
        self.end_time = rf.end_time
        # self.scan_index_addr = rf.scan_index_addr
        # self.scan_params_addr = rf.scan_params_addr
        # self.version = rf.version
        self.n_controllers = rf.n_controllers
        # self.data_addr = rf.data_addr
        self.n_filters = rf.n_filters
        self.timestamp = rf.timestamp
        self.scan_filters = self._get_scan_filters(
            rf.scan_filters
        )  # number of filters,
        if self.num_scans > 0:
            scans = np.frombuffer(
                (CScan * self.num_scans).from_address(rf.scans), CScan
            )
            self.scans = [Scan(s, profile, centroid) for s in scans]
        else:
            self.scans = []

        if self.num_chroms > 0:
            chromatograms = np.frombuffer(
                (CChrom * self.num_chroms).from_address(rf.chromatograms),
                CChrom,
            )
            self.chromatograms = [Chromatogram(c) for c in chromatograms]

        else:
            self.chromatograms = []

    def __del__(self):
        try:
            rawfilereader.FreeRawFile(self._data_p)
        except AttributeError:
            # raw file failed to initialize, just move on
            pass

    def _get_scan_filters(self, ptr):
        if not ptr:
            return []
        filters = np.copy(
            np.frombuffer((CScanFilter * self.n_filters).from_address(ptr), CScanFilter)
        )
        return [ScanFilter(f) for f in filters]

    def chrom_xic(self, mz, precursor, rt1, rt2):
        """
        Given a product mz and precursor mz and rt start and stop,
        returns a chromtogram. Look for +/- 0.5 mz for ranges
        """
        if self.num_chroms == 0:
            return Xic(None)
        for chrom in self.chromatograms:
            if (
                abs(chrom.precursor_mz - precursor) < 0.5
                and abs(chrom.product_mz - mz) < 0.5
            ):
                start = np.searchsorted(chrom.rts, rt1, side="left")
                stop = np.searchsorted(chrom.rts, rt2, side="right")
                if start == len(chrom.rts) or stop == 0:
                    return Xic(None)
                return Xic(
                    {
                        "rt": chrom.rts[start:stop],
                        "intensity": chrom.intensities[start:stop],
                        "length": stop - start,
                        "mz": chrom.product_mz,
                    }
                )
        return Xic(None)

    def xic(
        self,
        mz,
        rt1,
        rt2,
        ppm=20,
        scan_filter=None,
        centroid=True,
        pull_mzs=False,
        precursor=None,
    ):
        """
        Pulls an XIC
        """

        if self.num_scans == 0:
            return self.chrom_xic(mz, precursor, rt1, rt2)
        if scan_filter is None:
            scan_filter = -1

        if (centroid and not self.centroid) or (not centroid and not self.profile):
            empty = Empty()
            empty.intensity = np.array([])
            empty.rt = np.array([])
            empty.length = 0
            return empty
        return Xic(
            rawfilereader.Pull_xic(
                self._data_p, rt1, rt2, mz, ppm, scan_filter, centroid, pull_mzs
            )
        )

    def ms2(self, mz, rt1, rt2, ppm=5, centroid=True):
        """
        Pulls MS2 fragments for a given precursor MZ and RT range
        Returns fragment data and info about the scans and scan filters used
        """
        if not centroid:
            return None

        mz1 = mz - (mz / 1000000) * ppm
        mz2 = mz + (mz / 1000000) * ppm

        filters = []
        for i, scan_filter in enumerate(self.scan_filters):
            if (
                len(scan_filter.precursor_mz) > 0
                and mz1 <= scan_filter.precursor_mz[0] <= mz2
            ):
                filters.append(i)

        scan_list = []
        polarity = set()
        for i, scan in enumerate(self.scans):
            if scan.scan_filter in filters and rt1 <= scan.time <= rt2:
                scan_list.append(scan)
                polarity.add(self.scan_filters[scan.scan_filter].polarity)

        scan_info = {
            # filename without full path, if possible
            "rawfile": [self.filename.split("\\")[-1]],
            "scanList": [scan.index for scan in scan_list],
            "scanRts": [scan.time for scan in scan_list],
            "massRange": [0, 0],
            "precursorRange": [mz1, mz2],
            "resolution": [],
            "polarity": [p.name.lower() for p in polarity],
        }

        raw = []
        for scan in scan_list:
            for ms2_mz, intensity in zip(scan.cent_mzs, scan.cent_intensities):
                scan_filter = self.scan_filters[scan.scan_filter]
                scan_info["massRange"] = [scan_filter.low_mass, scan_filter.high_mass]
                raw.append(
                    {
                        "scanNumber": scan.index,
                        "rt": scan.time,
                        "mz": ms2_mz,
                        "intensity": intensity,
                        "chemicalFormula": "",
                        "formulaPpm": "",
                        "NCE": scan_filter.energy[0],
                        "aggregateNumber": "",
                    }
                )

        return raw, scan_info

    def export_ms2(
        self,
        mz,
        rt1,
        rt2,
        ppm=5,
        filename=None,
        rt_prec=4,
        mz_prec=5,
        int_prec=0,
        min_scans=1,
        include_raw=True,
    ):
        """
        Extracts, aggregates, and exports MS2 data in a standard MXP reporting format.
        Writes to csv if filename is provided, including extension (e.g. "results.csv").
        """
        raw, header = self.ms2(mz, rt1, rt2, ppm)
        raw, aggregated = self.aggregate(raw, ppm)

        output = [["HEADER"], ["version", 1]]
        for key, value in header.items():
            if isinstance(value, list):
                output.append([key, *value])
            else:
                output.append([key, value])

        output.append(["RAW"])
        if include_raw and len(raw) > 0:
            output.append(list(raw[0].keys()))
            for scan in raw:
                temp_output = list(scan.values())
                temp_output[1] = f"{temp_output[1]:.{rt_prec}f}"
                temp_output[2] = f"{temp_output[2]:.{mz_prec}f}"
                temp_output[3] = f"{temp_output[3]:.{int_prec}f}"
                output.append(temp_output)

        output.append(["AGGREGATED", "sum_v2"])
        if len(aggregated) > 0:
            output.append(list(aggregated[0].keys()))
            for group in aggregated:
                temp_output = list(group.values())
                if temp_output[6] < min_scans:
                    continue
                temp_output[1] = f"{temp_output[1]:.{mz_prec}f}"
                temp_output[2] = f"{temp_output[2]:.{int_prec}f}"
                output.append(temp_output)

        if filename:
            with open(filename, "w", encoding="ISO-8859-1") as file:
                for line in output:
                    file.write(",".join([str(item) for item in line]) + "\n")

        return output

    def aggregate(self, raw, ppm):
        """
        Aggregate MSn fragments with similar masses and equal NCEs.
        :param raw: list of dicts; MSn fragment data, including keys: "mz", "intensity",
                    "chemicalFormula", "formulaPpm", "NCE" at minimum
        :param ppm: float, permitted delta ppm for grouping MZs
        :return:
            raw: list of dicts, raw input updated with aggregateNumber
            aggregated: list of dicts, aggregated fragments
        """
        aggregated = []
        if not raw:
            return raw, aggregated
        raw_copy = raw.copy()  # Clearer than raw[:]
        raw_copy.sort(reverse=True, key=lambda e: e["intensity"])
        active = np.ones(len(raw_copy), dtype=bool)
        group_idx = 0
        data = np.array([(e["mz"], e["intensity"], e["NCE"]) for e in raw_copy])

        for i, raw_feature in enumerate(raw_copy):
            if not active[i]:
                continue

            center_mz = data[i, 0]
            ppm_factor = (center_mz / 1_000_000) * ppm
            min_mz, max_mz = center_mz - ppm_factor, center_mz + ppm_factor
            energy = data[i, 2]
            matched = (
                active
                & (data[:, 2] == energy)
                & (data[:, 0] >= min_mz)
                & (data[:, 0] <= max_mz)
            )
            if not matched.any():
                continue

            average_mz = np.average(data[matched, 0], weights=data[matched, 1])
            aggregated.append(
                {
                    "aggregateNumber": group_idx,
                    "mz": average_mz,
                    "intensity": data[matched, 1].sum(),
                    "chemicalFormula": raw_feature["chemicalFormula"],
                    "formulaPpm": raw_feature["formulaPpm"],
                    "NCE": energy,
                    "scanCount": matched.sum(),
                }
            )
            matched_indices = np.where(matched)[0]
            for j in matched_indices:
                raw_copy[j]["aggregateNumber"] = group_idx
            group_idx += 1
            active[matched_indices] = False
        return raw, aggregated


# testing functions
if __name__ == "__main__":
    from datetime import datetime, timedelta

    rf = RawFile(r"C:\Users\danie\Downloads\prepare.raw", profile=True)
    print(rf.timestamp)
    print(datetime(1601, 1, 1) + timedelta(microseconds=rf.timestamp / 10))
    print(datetime(1601, 1, 1) + timedelta(microseconds=rf.timestamp // 10))
