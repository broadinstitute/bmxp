"""A python implementation of an mzML file reader"""

import zlib
import base64
import struct
import warnings
from copy import copy
from xml.sax import parse
from xml.sax.handler import ContentHandler
import pynumpress
from rawfilereader import Analyzer, ScanMode, ScanType, Polarity, accessions
from rawfilereader.main_dataclasses import Scan, ScanFilter, RawFile

COMPRESSION_TYPES = {
    "no compression": lambda x: x,
    "zlib": zlib.decompress,
    "numpress linear": pynumpress.decode_linear,
    "numpress positive int": pynumpress.decode_pic,
    "numpress short logged float": pynumpress.decode_slof,
}


def parse_mzml(filename_or_stream):
    """Returns a Rawfile given a file string or stream
    :param filename_or_stream:string (filename) or file object
    :return: Rawfile
    """
    parser = MzmlParser()
    parse(filename_or_stream, parser)  # nosec
    return parser.rawfile


def parse_filter_string(filter_string, scan_filter=ScanFilter()):
    """
    Fills analyzer, scan type, scan mode, and dependent scan info from a filter string.
    :param filter_string: string, filter string saved from thermo raw files
    :param scan_filter: ScanFilter, a pre-existing ScanFilter to populate, defaults to
    new ScanFilter
    :return: ScanFilter, copy of the provided ScanFilter with additional info filled
    """
    filter_string = filter_string.lower()
    scan_filter = copy(scan_filter)
    for a in Analyzer:
        if a.name.lower() in filter_string:
            scan_filter.analyzer = a
            break
    for scan_type in ScanType:
        if scan_type.name.lower() in filter_string:
            scan_filter.scan_type = scan_type
    scan_filter.dependent_scans = " d " in filter_string

    if " p " in filter_string:
        scan_filter.scan_mode = ScanMode.PROFILE
    elif " c " in filter_string:
        scan_filter.scan_mode = ScanMode.CENTROID

    return scan_filter


class MzmlParser(ContentHandler):
    """mzML parser class to be used with xml.sax"""

    def __init__(self):
        super().__init__()
        self.stack = []
        self.current_function = None
        self.rawfile = RawFile()
        self.rawfile.low_mz = -1

        self.warn = False  # scan_mode warning should or shouldn't be raised

        self.referenceable_groups = {}
        self.instrument_config = {}

        self.start_method_map = {
            "mzML": self.set_version,
            "sourceFile": self.set_file_name,
            "cvParam": self.handle_cv_param,
            "spectrumList": self.set_spectrum_count,
            "spectrum": self.create_spectrum,
        }

        self.end_method_map = {
            "referenceableParamGroup": self.handle_referenceable_group,
            "instrumentConfiguration": self.handle_instrument_config,
            "spectrum": self.finalize_spectrum,
            "mzML": self.finalize_rawfile,
            "binary": self.handle_binary_data,
            "binaryDataArray": self.finalize_binary_data_array,
        }

    def startElement(self, name, attrs):
        # callbacks, if applicable
        if self.current_function:
            self.current_function(name, dict(attrs))  # pylint: disable=not-callable
        elif name in self.start_method_map:
            method = self.start_method_map[name]
            method(name, dict(attrs))

        if name in ("indexedmzML", "mzML"):
            return
        element = {"name": name, "attrs": dict(attrs), "value": "", "children": []}
        if len(self.stack) > 0:
            self.stack[-1]["children"].append(element)
        self.stack.append(element)

    def endElement(self, name):
        if name in self.end_method_map:
            self.end_method_map[name]()

        if name in ("indexedmzML", "mzML"):
            return
        self.stack.pop()

    def characters(self, content):
        if content.strip() != "":
            self.stack[-1]["value"] += repr(content)

    def handle_cv_param(self, _, attrs):
        """Callback to handle generic cvParam elements"""
        acc = attrs["accession"]
        if acc in accessions.INSTRUMENTS:
            self.rawfile.instrument_model = attrs["name"]

    def set_version(self, _, attrs):
        """Callback for mzML element. Sets file version."""
        self.rawfile.version = attrs["version"]

    def set_file_name(self, _, attrs):
        """Callback for sourceFile. Sets file name."""
        if len(self.stack[-1]["children"]) == 0:  # first file only
            self.rawfile.file_name = attrs["name"]

    def handle_referenceable_group(self):
        """Callback for referenceableParamGroup. Stores group contents in parser."""
        group = self.stack[-1]
        self.referenceable_groups[group["attrs"]["id"]] = group["children"]

    def handle_instrument_config(self):
        """Callback for instrumentConfiguration. Stores group contents in parser."""
        group = self.stack[-1]
        self.instrument_config[group["attrs"]["id"]] = group["children"]

    def set_spectrum_count(self, _, attrs):
        """Callback for spectrumList. Sets num_scans."""
        self.rawfile.num_scans = int(attrs["count"])

    def create_spectrum(self, _, attrs):
        """Callback for spectrum element. Creates new Scan."""
        # create and store a temporary ScanFilter within the Scan
        self.rawfile.scans.append(Scan(index=int(attrs["index"]), filter=ScanFilter()))
        self.current_function = self.handle_scan_cv_params

    def handle_scan_cv_params(self, name, attrs):  # pylint: disable=too-many-branches
        """Callback for elements within a spectrum. Sets various Scan fields."""
        if name == "binaryDataArray":
            self.current_function = self.handle_binary_data_array
        if name != "cvParam":
            return
        acc = attrs["accession"]
        unit = None if "unitAccession" not in attrs else attrs["unitAccession"]
        current_scan = self.rawfile.scans[-1]
        if acc == accessions.POSITIVE_SCAN:
            current_scan.filter.polarity = Polarity.POSITIVE
        elif acc == accessions.NEGATIVE_SCAN:
            current_scan.filter.polarity = Polarity.NEGATIVE
        elif acc == accessions.MS_LEVEL:
            current_scan.filter.ms_level = int(attrs["value"])
        elif acc == accessions.PROFILE:
            current_scan.filter.scan_mode = ScanMode.PROFILE
        elif acc == accessions.CENTROID:
            current_scan.filter.scan_mode = ScanMode.CENTROID
        elif acc == accessions.TIME:
            current_scan.time = float(attrs["value"])
        elif acc == accessions.FILTER_STRING:
            current_scan.filter.filter_string = attrs["value"]
        elif acc == accessions.LOW_MZ and unit == accessions.UNIT_MZ:
            current_scan.filter.low_mass = float(attrs["value"])
            if self.rawfile.low_mz < 0 or float(attrs["value"]) < self.rawfile.low_mz:
                self.rawfile.low_mz = float(attrs["value"])
        elif acc == accessions.HIGH_MZ and unit == accessions.UNIT_MZ:
            current_scan.filter.high_mass = float(attrs["value"])
            if float(attrs["value"]) > self.rawfile.high_mz:
                self.rawfile.high_mz = float(attrs["value"])
        # consider "selected ion m/z" MS:1000744 if isolation window is not reliable
        elif acc == accessions.ISOLATION_WINDOW and unit == accessions.UNIT_MZ:
            current_scan.filter.precursor_mz.append(float(attrs["value"]))
        elif acc == accessions.COLLISION_ENERGY:
            current_scan.filter.energy.append(float(attrs["value"]))

    def finalize_spectrum(self):
        """
        Callback for spectrum element. Finalizes ScanFilter and sets default values if
        needed.
        """
        self.current_function = None
        current_scan = self.rawfile.scans[-1]

        if len(current_scan.filter.precursor_mz) == 0:
            current_scan.filter.precursor_mz = [
                0 for _ in range(current_scan.filter.ms_level)
            ]
        if len(current_scan.filter.energy) == 0:
            current_scan.filter.energy = [
                0 for _ in range(current_scan.filter.ms_level)
            ]
        if len(current_scan.intensities) == 0:  # no profile data
            current_scan.intensities = current_scan.cent_intensities
            current_scan.mzs = current_scan.cent_mzs

        new_filter = parse_filter_string(
            current_scan.filter.filter_string, current_scan.filter
        )
        if (
            new_filter.scan_mode == ScanMode.PROFILE
            and current_scan.filter.scan_mode == ScanMode.CENTROID
        ):
            self.warn = True
        current_scan.filter = new_filter

        if current_scan.filter not in self.rawfile.scan_filters:
            self.rawfile.scan_filters.append(current_scan.filter)
        current_scan.filter = self.rawfile.scan_filters.index(current_scan.filter)

    def handle_binary_data_array(self, name, attrs):
        """
        Callback for binaryDataArray within a spectrum. Sets units, compression, and
        binary data.
        """
        if name != "cvParam":
            return

        acc = attrs["accession"]
        if acc in accessions.BINARY_UNITS:
            self.stack[-1]["unit"] = accessions.BINARY_UNITS[acc]
        elif acc in accessions.BINARY_COMPRESSION:
            self.stack[-1]["compression"] = accessions.BINARY_COMPRESSION[acc]
        elif acc == accessions.INTENSITY_ARRAY:
            self.stack[-1]["array type"] = "intensity"
        elif acc == accessions.MZ_ARRAY:
            self.stack[-1]["array type"] = "mz array"

    def handle_binary_data(self):
        """Callback for binary element"""
        if self.stack[-2]["name"] != "binaryDataArray":
            return
        # copy binary data to parent in stack
        self.stack[-2]["binary"] = self.stack[-1]["value"]

    def finalize_binary_data_array(self):
        """
        Callback for binaryDataArray. Decompresses and decodes m/z and intensity data.
        """
        # ignore chromatograms
        if len(self.stack) < 3 or self.stack[-3]["name"] != "spectrum":
            return
        self.current_function = self.handle_scan_cv_params  # back to spectrum
        data_array = self.stack[-1]
        if "array type" not in data_array:
            # is not intensity or mz, can be skipped
            return
        if "unit" not in data_array:
            raise RuntimeError("Binary unit missing or not supported.")
        if "compression" not in data_array:
            raise RuntimeError("Binary compression type missing or not supported.")

        binary = base64.b64decode(self.stack[-1]["binary"])
        for compression in data_array["compression"]:
            func = COMPRESSION_TYPES[compression]
            binary = func(binary)

        unit = data_array["unit"]
        if "numpress" in data_array["compression"][-1]:
            decoded = binary
        elif unit:
            format_string = f"<{int(len(binary) / unit[1])}{unit[0]}"
            decoded = struct.unpack(format_string, binary)
        else:
            decoded = binary.decode("ascii")

        if (
            data_array["array type"] == "intensity"
            and self.rawfile.scans[-1].filter.scan_mode == ScanMode.PROFILE
        ):
            self.rawfile.scans[-1].intensities = list(decoded)
        elif data_array["array type"] == "intensity":
            self.rawfile.scans[-1].cent_intensities = list(decoded)
        elif (
            data_array["array type"] == "mz array"
            and self.rawfile.scans[-1].filter.scan_mode == ScanMode.PROFILE
        ):
            self.rawfile.scans[-1].mzs = list(decoded)
        else:
            self.rawfile.scans[-1].cent_mzs = list(decoded)

    def finalize_rawfile(self):
        """Callback for close of mzML element"""
        self.rawfile.file_format = "mzML"
        self.rawfile.first_scan_number = self.rawfile.scans[0].index
        self.rawfile.last_scan_number = self.rawfile.scans[-1].index
        self.rawfile.start_time = self.rawfile.scans[0].time
        self.rawfile.end_time = self.rawfile.scans[-1].time

        if self.warn:
            warnings.warn(
                "One or more scans are labeled as being acquired in profile mode but "
                + "only contain centroid data."
            )
