"""Definitions of dataclasses used by the python rawfilereader"""

# pylint: disable=too-many-instance-attributes,missing-class-docstring
from dataclasses import dataclass, field
from typing import List
from rawfilereader import Analyzer, ScanMode, ScanType, Polarity


@dataclass
class Xic:
    rt: float
    intensity: float
    length: int
    mz: float


@dataclass
class ProfileChunk:
    first_bin: int = 0
    n_bins: int = 0
    fudge: float = 0
    signal: List[float] = field(default_factory=list)


@dataclass
class Scan:
    intensities: List[float] = field(default_factory=list)
    mzs: List[float] = field(default_factory=list)
    peak_total: int = 0  # granular total of all mz/ints

    cent_intensities: List[float] = field(default_factory=list)
    cent_mzs: List[float] = field(default_factory=list)

    chunks: List[ProfileChunk] = field(default_factory=list)

    first_value: float = 0
    step: float = 0
    peak_count: int = 0  # number of chunks (peaks?)
    offset: int = 0
    a: float = 0
    b: float = 0
    c: float = 0  # pylint: disable=invalid-name
    count: int = 0  # centroided peak count
    peak_list_size: int = 0
    layout: int = 0
    n_param: int = 0
    n_bins: int = 0
    profile_size: int = 0
    time: float = 0
    index: int = 0
    filter: int = -1  # index of matching filter in RawFile.scan_filters


@dataclass
class ScanFilter:
    polarity: Polarity = Polarity.UNDEFINED
    scan_mode: ScanMode = ScanMode.UNDEFINED
    dependent_scans: bool = None
    ms_level: int = None
    scan_type: ScanType = ScanType.UNDEFINED
    analyzer: Analyzer = Analyzer.UNDEFINED
    precursor_mz: List[float] = field(default_factory=list)
    energy: List[float] = field(default_factory=list)
    low_mass: float = None
    high_mass: float = None
    filter_string: str = ""


@dataclass
class RawFile:
    data: bytes = b""
    instrument_model: str = ""
    file_name: str = ""
    file_format: str = "raw"
    scans: List[Scan] = field(default_factory=list)
    num_scans: int = 0
    scan_trailer_addr: int = 0
    first_scan_number: int = 0
    last_scan_number: int = 0
    low_mz: float = 0
    high_mz: float = 0
    start_time: float = 0
    end_time: float = 0
    scan_index_addr: int = 0
    scan_params_addr: int = 0
    version: int = 0
    n_controllers: int = 0
    data_addr: int = 0
    scan_filters: List[ScanFilter] = field(default_factory=list)
    timestamp: int = 0
