# pylint: disable-all
import numpy as np
from bmxp.chroma import RawFile

if __name__ == "__main__":
    folder = "D:\\Work\\Tests\\RawFileReader\\"

    print("**** Testing Orbitrap MS1 ****")
    # Test full scan ms1 data
    alanine_rt = 7.73
    alanine_mz = 90.0554
    rt_dif = 0.02
    ppm = 5
    rts = np.array(
        [
            7.7116547,
            7.716111,
            7.720567,
            7.7250257,
            7.7294817,
            7.7339406,
            7.7383966,
            7.7428546,
            7.747311,
        ]
    )
    intensities = (
        np.array(
            [
                27994996.0,
                32120916.0,
                36518720.0,
                38417432.0,
                38131832.0,
                45260688.0,
                40830152.0,
                37549100.0,
                34188064.0,
            ]
        ),
    )

    pr_intensities = np.array(
        [
            9.6551744e07,
            1.1090329e08,
            1.2637120e08,
            1.3285288e08,
            1.3186667e08,
            1.5642437e08,
            1.4103850e08,
            1.2982436e08,
            1.1825154e08,
        ]
    )

    # test full scan ms1
    ms1_files = [
        "QE-ihmp2-rx0020a-hp.raw",
        "32-Centroid-zlib-NLC-QE.mzML",
        "32-Centroid-NoCompression-QE.mzML",
        "64-Centroid-NoCompression.mzML",
        "64-Centroid-NLC-QE.mzML",
        "64-Centroid-SLFC-QE.mzML",
        "64-Centroid-zlib-NLC-QE.mzML",
        "64-Centroid-zlib-QE.mzML",
    ]
    for ms1_file in ms1_files:
        print(ms1_file, " ", end="")
        profile = False
        centroid = True
        rf = RawFile(folder + ms1_file, profile=profile, centroid=centroid)
        xic = rf.xic(
            alanine_mz, alanine_rt - rt_dif, alanine_rt + rt_dif, ppm, centroid=centroid
        )
        print(
            np.isclose(xic.rt, rts).all(),
            np.isclose(
                xic.intensity,
                intensities,
            ).all(),
        )

    ms1_files = [
        "64-Profile-NoCompression.mzML",
        "32-Profile-NoCompression.mzML",
    ]
    for ms1_file in ms1_files:
        print(ms1_file, " ", end="")
        profile = True
        centroid = False
        rf = RawFile(folder + ms1_file, profile=profile, centroid=centroid)
        xic = rf.xic(
            alanine_mz, alanine_rt - rt_dif, alanine_rt + rt_dif, ppm, centroid=centroid
        )
        print(
            np.isclose(xic.rt, rts).all(),
            np.isclose(
                xic.intensity,
                pr_intensities,
            ).all(),
        )

    print("**** Testing QQQ Transitions****")
    # test triple quad
    filename = "Sciex-QTRAP-6500.mzML"
    print(filename, " ", end="")
    rf = RawFile(folder + filename)
    xic = rf.xic(369.4, 8.5, 8.6, precursor=675.2)
    print(np.isclose(xic.rt[0], 8.512), np.isclose(xic.intensity[0], 15604))

    # # test triple quad
    filename = "Agilent-6495-QQQ.mzML"
    print(filename, " ", end="")
    rf = RawFile(folder + filename)
    xic = rf.xic(44.296, 7.8, 7.9, precursor=90.06)
    print(np.isclose(xic.rt[0], 7.8011165), np.isclose(xic.intensity[0], 52846.855))
