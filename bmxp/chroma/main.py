# pylint: skip-file
"""A python implementation of the C finniganreader"""

import struct
from rawfilereader import Analyzer, ScanMode, ScanType, Polarity
from rawfilereader.main_dataclasses import Xic, ProfileChunk, Scan, ScanFilter, RawFile
from rawfilereader.mzml_parser import parse_mzml


def from_pascal(data, pos):
    """
    Decode a pascal string from bytes
    :param data: bytes
    :param pos: index of the beginning of the string
    :return: string
    """
    [length] = struct.unpack("I", data[pos : pos + 4])
    pos += 4
    pascal = []
    for _ in range(length):
        pascal += data[pos : pos + 1]
        pos += 2
    return bytes(pascal).decode("utf-8")


def move_pascal(data, pos, n):
    """
    Skips past n pascal strings, beginning at index=pos
    :param data: bytes
    :param pos: index of the beginning of the first string
    :param n: number of strings to skip over
    :return: position after skipping over n strings
    """
    for _ in range(n):
        [length] = struct.unpack("I", data[pos : pos + 4])
        pos += 4 + length * 2
    return pos


def fill_scan_events(rawfile):
    """
    Create Scans and populate with basic scan event info
    :param rawfile:
    """
    pos = rawfile.scan_trailer_addr + 4

    for _ in range(rawfile.num_scans):
        scan_i = Scan()
        scan_filter = ScanFilter()
        # scan_mode at pos+5, ms_level at pos+6, scan_type at pos+7, data_dependent at
        # pos+10, analyzer at pos+40
        [
            polarity,
            scan_mode,
            scan_filter.ms_level,
            scan_type,
            dependent_scans,
            analyzer,
        ] = struct.unpack("4B2xB29xB", rawfile.data[pos + 4 : pos + 41])

        scan_filter.polarity = Polarity(polarity)
        scan_filter.scan_mode = ScanMode(scan_mode)
        scan_filter.scan_type = ScanType(scan_type)
        scan_filter.analyzer = Analyzer(analyzer)
        scan_filter.dependent_scans = bool(dependent_scans)

        pos += 136
        [n_precursors] = struct.unpack("I", rawfile.data[pos : pos + 4])
        pos += 4

        # in non-ID-X data, MS1 data has n_precursors=1 and 56 bytes of unknowns to skip
        # in ID-X data, MS1 data has n_precursors=0 and 56 bytes of data/unknowns
        for _ in range(n_precursors):
            if scan_filter.ms_level > 1:
                scan_filter.precursor_mz.append(
                    round(struct.unpack("d", rawfile.data[pos : pos + 8])[0], 4)
                )
                scan_filter.energy.append(
                    struct.unpack("d", rawfile.data[pos + 16 : pos + 24])[0]
                )
            pos += 56

        [scan_filter.low_mass, scan_filter.high_mass, scan_i.n_param] = struct.unpack(
            "ddI", rawfile.data[pos + 4 : pos + 24]
        )
        pos += 24

        if scan_i.n_param > 0:
            [scan_i.a, scan_i.b, scan_i.c] = struct.unpack(
                "3d", rawfile.data[pos + 16 : pos + 40]
            )

        if rawfile.instrument_model == "Orbitrap Exploris 240":
            pos += 48
            pos = move_pascal(rawfile.data, pos, 1)
        elif scan_i.n_param == 7:  # Orbitrap ID-X
            pos += 68
        elif scan_i.n_param == 5:
            pos += 60
        elif scan_i.n_param == 0:
            pos += 12

        for j, filter_j in enumerate(rawfile.scan_filters):
            if filter_j == scan_filter:
                scan_i.filter = j
                break
        if scan_i.filter == -1:
            scan_i.filter = len(rawfile.scan_filters)
            rawfile.scan_filters.append(scan_filter)

        scan_i.filter = rawfile.scan_filters.index(scan_filter)

        rawfile.scans.append(scan_i)


def fill_scan_indices(rawfile):
    """
    Populate scans with index, time, offset
    :param rawfile:
    """
    pos = rawfile.scan_index_addr
    for i in range(rawfile.num_scans):
        # index at pos+4, time at pos+24, offset at pos+72
        [
            rawfile.scans[i].index,
            rawfile.scans[i].time,
            rawfile.scans[i].offset,
        ] = struct.unpack("=I16xd40xQ", rawfile.data[pos + 4 : pos + 80])
        pos += 88


def read_scan_data_packet(rawfile, index, scan, profile):
    """
    Populate scan with actual peaks/data
    :param rawfile:
    :param index: address/index where scan data begins (scan.offset)
    :param scan: scan to be read
    """
    [scan.profile_size, scan.peak_list_size, scan.layout] = struct.unpack(
        "3I", rawfile.data[index + 4 : index + 16]
    )
    index += 40
    total_peaks = 0
    if scan.profile_size > 0:
        [scan.first_value, scan.step, scan.peak_count, scan.n_bins] = struct.unpack(
            "ddII", rawfile.data[index : index + 24]
        )
        index += 24
        if not profile:
            for _ in range(scan.peak_count):
                chunk = ProfileChunk()
                [chunk.first_bin, chunk.n_bins] = struct.unpack(
                    "II", rawfile.data[index : index + 8]
                )
                index += 8
                if scan.layout > 0:
                    index += 4
                index += 4 * chunk.n_bins
        else:
            for _ in range(scan.peak_count):
                chunk = ProfileChunk()
                [chunk.first_bin, chunk.n_bins] = struct.unpack(
                    "II", rawfile.data[index : index + 8]
                )
                total_peaks += chunk.n_bins
                index += 8
                if scan.layout > 0:
                    [chunk.fudge] = struct.unpack("f", rawfile.data[index : index + 4])
                    index += 4
                chunk.signal = list(
                    struct.unpack(
                        f"{chunk.n_bins}f",
                        rawfile.data[index : index + 4 * chunk.n_bins],
                    )
                )
                index += 4 * chunk.n_bins
                scan.chunks.append(chunk)

    if scan.peak_list_size > 0:
        [scan.count] = struct.unpack("I", rawfile.data[index : index + 4])
        index += 4
        for _ in range(scan.count):
            [mz, intensity] = struct.unpack("ff", rawfile.data[index : index + 8])
            scan.cent_intensities.append(intensity)
            scan.cent_mzs.append(mz)
            index += 8

    scan.peak_total = total_peaks


def fill_scans(rawfile, profile, centroid):
    """
    Populate scans with usable MZs/intensities
    :param rawfile:
    """
    for i in range(rawfile.num_scans):
        scan = rawfile.scans[i]
        read_scan_data_packet(rawfile, scan.offset, scan, profile)
        if not profile:
            scan.intensities = []
            scan.mzs = []
            continue
        # convert profile readings into mz/intensities
        for peak_n in range(scan.peak_count):
            for bin_n in range(scan.chunks[peak_n].n_bins):
                # convert Hz values to m/z
                v = (
                    scan.first_value
                    + (scan.chunks[peak_n].first_bin + bin_n) * scan.step
                )
                if scan.n_param == 4:
                    scan.mzs.append(
                        scan.a + scan.b / v + scan.c / v / v + scan.chunks[peak_n].fudge
                    )
                elif scan.n_param in (5, 7):
                    scan.mzs.append(
                        scan.a
                        + scan.b / v / v
                        + scan.c / v / v / v / v
                        + scan.chunks[peak_n].fudge
                    )
                else:
                    scan.mzs.append(v + scan.chunks[peak_n].fudge)
            scan.intensities.extend(scan.chunks[peak_n].signal)
        scan.chunks = []

        # # put centroided data in mz/intensity lists if there's no profile data
        # if scan.peak_count == 0:
        #     scan.mzs = scan.cent_mzs.copy()
        #     scan.intensities = scan.cent_intensities.copy()


def initialize(rawfile, profile, centroid):
    """
    Fetch all relevant info from the provided rawfile
    :param rawfile:
    """
    # Skip to the version number in FileHeader
    pos = 36
    [rawfile.version] = struct.unpack("I", rawfile.data[pos : pos + 4])
    [rawfile.timestamp] = struct.unpack("I", rawfile.data[pos + 4 : pos + 12])
    print(f"version: {rawfile.version}")
    if rawfile.version < 66:
        raise ValueError("Versions under 66 are not supported.")
    pos = 1420  # start right after InjectionData in SequencerRow
    pos = move_pascal(rawfile.data, pos, 11)  # skip labels, comments, method files
    rawfile.file_name = from_pascal(rawfile.data, pos)
    pos = move_pascal(rawfile.data, pos, 5)  # skip filename, path, vial name
    pos += 4
    pos = move_pascal(rawfile.data, pos, 15)
    pos += 24  # AutoSamplerInfo
    pos = move_pascal(rawfile.data, pos, 1)
    pos += 28  # skip the methodfile and timestamp
    [rawfile.n_controllers] = struct.unpack("I", rawfile.data[pos : pos + 4])
    pos += 780  # skip other stuff in the preamble and the padding
    [rawfile.data_addr] = struct.unpack("Q", rawfile.data[pos : pos + 8])
    pos += 16
    run_header_addr = []
    for i in range(rawfile.n_controllers):
        run_header_addr.append(struct.unpack("Q", rawfile.data[pos : pos + 8])[0])
        pos += 16  # skip runHeaderAddr and the unknown 7
    rawfile.scan_trailer_addr = 0
    for i in range(rawfile.n_controllers):
        if rawfile.scan_trailer_addr != 0:
            break
        pos = run_header_addr[i]
        [rawfile.first_scan_number, rawfile.last_scan_number] = struct.unpack(
            "II", rawfile.data[pos + 8 : pos + 16]
        )
        [
            rawfile.low_mz,
            rawfile.high_mz,
            rawfile.start_time,
            rawfile.end_time,
        ] = struct.unpack("4d", rawfile.data[pos + 56 : pos + 88])
        [rawfile.scan_index_addr] = struct.unpack(
            "Q", rawfile.data[pos + 7408 : pos + 7416]
        )
        [rawfile.scan_trailer_addr, rawfile.scan_params_addr] = struct.unpack(
            "QQ", rawfile.data[pos + 7448 : pos + 7464]
        )
        pos += 7588
        pos = move_pascal(rawfile.data, pos, 1)
        rawfile.instrument_model = from_pascal(rawfile.data, pos)

    # pull the scan Events
    rawfile.num_scans = rawfile.last_scan_number - rawfile.first_scan_number + 1

    if rawfile.n_controllers == 0 or rawfile.num_scans == 0:
        raise ValueError("Invalid raw file.")

    fill_scan_events(rawfile)
    fill_scan_indices(rawfile)
    for j in range(rawfile.num_scans):
        rawfile.scans[j].offset += rawfile.data_addr
    fill_scans(rawfile, profile, centroid)


def find_first_ge(arr, target_rt):
    """
    Find the first RT that is greater than or equal to the target RT
    :param arr: list of scans
    :param target_rt: RT to search for
    :return: int, index of resulting RT
    """
    n = len(arr)
    if arr[0].time >= target_rt:
        return 0
    if n == 1:
        return 1
    midpoint = int(n / 2)
    if arr[midpoint].time >= target_rt:
        return find_first_ge(arr[:midpoint], target_rt)
    return midpoint + find_first_ge(arr[midpoint:], target_rt)


def find_first_ge_float(arr, target):
    """
    find_first_ge, but generalized to find any target in any list
    :param arr: list of values to search
    :param target: target value to find
    :return: int, index of resulting value
    """
    n = len(arr)
    if n == 0 or arr[0] >= target:
        return 0
    if n == 1:
        return 1
    midpoint = int(n / 2)
    if arr[midpoint] >= target:
        return find_first_ge_float(arr[:midpoint], target)
    return midpoint + find_first_ge_float(arr[midpoint:], target)


def pull_xic(
    rawfile, rt1, rt2, mz, ppm=20, scan_filter=None, centroid=True, pull_mzs=False
):
    """
    Generate RTs/intensities for a given MZ and given RT range
    :param rawfile:
    :param rt1: float, lower end of RT range
    :param rt2: float, upper end of RT range
    :param mz: float, target MZ
    :param ppm: float, allowed delta ppm from target MZ
    :param scan_filter: int, index of filter in RawFile.scan_filters list
    :param centroid: bool, use centroided peaks instead of profile data
    :return: Xic containing extracted data
    """
    start = find_first_ge(rawfile.scans, rt1)
    stop = find_first_ge(rawfile.scans, rt2) - 1
    length = stop - start + 1
    intensity = []
    rt = []
    return_mzs = []

    for i in range(length):
        if scan_filter is None or rawfile.scans[start + i].filter == scan_filter:
            min_mz = mz - (mz / 1000000) * ppm
            max_mz = mz + (mz / 1000000) * ppm

            rt.append(rawfile.scans[start + i].time)
            intensity.append(0)
            return_mzs.append(0)

            if centroid:
                mzs = rawfile.scans[start + i].cent_mzs
                intensities = rawfile.scans[start + i].cent_intensities
            else:
                mzs = rawfile.scans[start + i].mzs
                intensities = rawfile.scans[start + i].intensities

            peak_start = find_first_ge_float(mzs, min_mz)
            peak_stop = find_first_ge_float(mzs, max_mz) - 1

            if peak_start > peak_stop:
                continue
            avg_mz = 0
            for j in range(peak_stop - peak_start + 1):
                intensity[-1] += intensities[j + peak_start]
                if pull_mzs:
                    avg_mz += intensities[j + peak_start] * mzs[j + peak_start]
            if pull_mzs:
                try:
                    return_mzs[-1] = avg_mz / intensity[-1]
                except ZeroDivisionError:
                    return_mzs[-1] = 0

    return Xic(rt=rt, intensity=intensity, length=length, mz=return_mzs)


def from_bytes(data, profile=False, centroid=True):
    """
    Generate a rawfile from bytes
    :param data: bytes, rawfile contents
    :return: rawfile
    """
    rawfile = RawFile(data=data)
    initialize(rawfile, profile, centroid)
    rawfile.data = b""
    return rawfile


def open_rawfile(filename, profile=False, centroid=True):
    """
    Open a rawfile from a file name
    :param filename: string, filename
    :return: rawfile
    """
    with open(filename, mode="rb") as file:
        data = file.read()
    rawfile = from_bytes(data, profile, centroid)
    return rawfile


def open_mzml(filename_or_stream):
    """
    Open an mzML file from a file name. Wraps parse_mzml so it can be accessed from the
    main module.
    :param filename_or_stream: string (filename) or file object
    :return: rawfile
    """
    return parse_mzml(filename_or_stream)


if __name__ == "__main__":
    # for i in range(1):
    #     print(i)
    #     rawfile = open_rawfile(r"D:\Work\TFE Demo\PREFA05.raw")
    # exit()
    # rawfile = open_rawfile(r"D:\Work\TFE Demo\PREFA05.raw")

    rawfile = open_rawfile(r"D:\Work\TFE Demo\PREFA05.raw", profile=True, centroid=True)
    result = pull_xic(rawfile, 9.76, 9.8, 89.1073)
    print(result.intensity)
    print(result.mz)
    result = pull_xic(rawfile, 9.76, 9.8, 89.1073, pull_mzs=True)
    print(result.intensity)
    print(result.mz)
    result = pull_xic(rawfile, 9.76, 9.8, 89.1073, centroid=False)
    print(result.intensity)
    print(result.mz)
    result = pull_xic(rawfile, 9.76, 9.8, 89.1073, centroid=False, pull_mzs=True)
    print(result.intensity)
    print(result.mz)
