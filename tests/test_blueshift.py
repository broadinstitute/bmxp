# pylint: disable=redefined-outer-name, missing-function-docstring, consider-using-with
"""
Tests for blueshift
"""
import pickle
from pathlib import Path
import pytest
import pandas as pd
import numpy as np
from bmxp import blueshift as b


@pytest.fixture()
def path_dc_input_1():
    return Path(__file__).parent / "DCinput1.csv"


@pytest.fixture()
def path_sample_info_1():
    return Path(__file__).parent / "DCinfo1.csv"


@pytest.fixture()
def path_dc_input_2():
    return Path(__file__).parent / "DCinput2.csv"


@pytest.fixture()
def path_sample_info_2():
    return Path(__file__).parent / "DCinfo2.csv"


@pytest.fixture()
def df_dc_input_1(path_dc_input_1):
    return pd.read_csv(path_dc_input_1)


@pytest.fixture()
def df_sample_info_1(path_sample_info_1):
    return pd.read_csv(path_sample_info_1)


@pytest.fixture()
def df_dc_input_2(path_dc_input_2):
    return pd.read_csv(path_dc_input_2)


@pytest.fixture()
def df_sample_info_2(path_sample_info_2):
    return pd.read_csv(path_sample_info_2)


@pytest.fixture()
def pickled_results():
    return pickle.load(open(Path(__file__).parent / "blueshift.pickle", "rb"))


def test_data_validation(df_dc_input_1, df_sample_info_1):
    # missing required column in injection information
    info = df_sample_info_1.drop("injection_order", axis=1)
    with pytest.raises(ValueError) as e:
        b.DriftCorrection(df_dc_input_1, info)
    assert "injection_order" in str(e.value)

    # missing injection in data input
    data = df_dc_input_1.drop("B0005_COL_ExampleProject_CN-M36058078", axis=1)
    with pytest.raises(ValueError) as e:
        b.DriftCorrection(data, df_sample_info_1)
    assert "data sheet: B0005_COL_ExampleProject_CN-M36058078" in str(e.value)

    # no error when missing "not_used" injection in data input
    data = df_dc_input_1.drop("B0008_COL_ExampleProject_CN-M59244903", axis=1)
    b.DriftCorrection(data, df_sample_info_1)

    # duplicate injection order
    info = df_sample_info_1.copy()
    info.loc[14, "injection_order"] = info.loc[15, "injection_order"]
    with pytest.raises(ValueError) as e:
        b.DriftCorrection(df_dc_input_1, info)
    assert "duplicate values" in str(e.value)

    # duplicate injection id
    info = df_sample_info_1.copy()
    info.loc[14, "injection_id"] = info.loc[15, "injection_id"]
    with pytest.raises(ValueError) as e:
        b.DriftCorrection(df_dc_input_1, info)
    assert "duplicate injection_ids" in str(e.value)

    # out-of-order injection order
    info = df_sample_info_1.copy()
    info.loc[14, "injection_order"] = 700
    with pytest.raises(ValueError) as e:
        b.DriftCorrection(df_dc_input_1, info)
    assert "must be sorted" in str(e.value)

    # invalid label in batches column
    info = df_sample_info_1.copy()
    info.loc[13, "batches"] = "batch nd"
    with pytest.raises(ValueError) as e:
        b.DriftCorrection(df_dc_input_1, info)
    assert "invalid label" in str(e.value)

    # non-numeric character in data
    data = df_dc_input_1.copy()
    data.iloc[5, 5] = "f"
    with pytest.raises(TypeError) as e:
        b.DriftCorrection(data, df_sample_info_1)
    assert "non-numeric" in str(e.value)

    # data and sample are not in same order
    data = df_dc_input_1.copy()
    col_list = list(data.columns)
    col_list = col_list[:10] + col_list[11:] + col_list[10:11]
    data = data.loc[:, col_list]
    with pytest.raises(ValueError) as e:
        b.DriftCorrection(data, df_sample_info_1)
    assert "usable samples" in str(e.value)


def test_batch_generation(
    df_dc_input_1,
    df_sample_info_1,
    path_dc_input_2,
    path_sample_info_2,
    pickled_results,
):
    a = b.DriftCorrection(df_dc_input_1, df_sample_info_1)
    for batch, ref_batch in zip(a.batches["default"], pickled_results["default1"]):
        assert (batch.values == ref_batch.values).all()
    for batch, ref_batch in zip(a.batches["override"], pickled_results["override1"]):
        assert (batch.values == ref_batch.values).all()

    a = b.DriftCorrection(path_dc_input_2, path_sample_info_2)
    for batch, ref_batch in zip(a.batches["default"], pickled_results["default2"]):
        assert (batch.values == ref_batch.values).all()
    for batch, ref_batch in zip(a.batches["override"], pickled_results["override2"]):
        assert (batch.values == ref_batch.values).all()


def test_internal_standard_correction(
    path_dc_input_1,
    df_dc_input_1,
    path_sample_info_1,
    df_dc_input_2,
    df_sample_info_2,
    pickled_results,
):
    # one internal standard
    a = b.DriftCorrection(path_dc_input_1, path_sample_info_1)
    a.internal_standard_correct("Internal Standard 1")
    assert np.isclose(
        a.data.round().fillna(0),
        pickled_results["DCinput1_IS_InternalStandard1"].round().loc[:, a.data.columns],
        equal_nan=True,
    ).all()

    # one internal standard with nonquant duplicate
    nonquant_df = df_dc_input_1.copy()
    nonquant_df.loc[4, "Metabolite"] = "Internal Standard 1"
    nonquant_df.loc[4, "Non_Quant"] = True
    a = b.DriftCorrection(nonquant_df, path_sample_info_1)
    a.internal_standard_correct("Internal Standard 1")
    assert np.isclose(
        a.data.round().fillna(0),
        pickled_results["DCinput1_IS_InternalStandard1"].round().loc[:, a.data.columns],
        equal_nan=True,
    ).all()

    # nonquant "missing" internal standard
    nonquant_df = df_dc_input_1.copy()
    nonquant_df.loc[0, "Non_Quant"] = True
    a = b.DriftCorrection(nonquant_df, path_sample_info_1)
    with pytest.raises(ValueError) as e:
        a.internal_standard_correct("Internal Standard 1")
    assert "not found in" in str(e.value)

    # two internal standards
    a = b.DriftCorrection(df_dc_input_2, df_sample_info_2)
    a.internal_standard_correct(["15R-15-methyl-PGA2", "15R-15-methyl-PGF2a"])
    assert np.isclose(
        a.data.round(),
        pickled_results["DCinput2_IS_PGA2_PGF2a"].loc[:, a.data.columns],
        equal_nan=True,
    ).all()

    # missing IS value
    data = df_dc_input_2.copy()
    data.iloc[14, 50] = 0
    a = b.DriftCorrection(data, df_sample_info_2)
    with pytest.raises(ValueError) as e:
        a.internal_standard_correct("15S-15-methyl-PGD2")
    assert "missing values" in str(e.value)

    # wrong IS name
    with pytest.raises(ValueError) as e:
        a.internal_standard_correct("not_a_real_metabolite")
    assert "not found in" in str(e.value)


def test_pool_correction(
    path_dc_input_1,
    path_sample_info_1,
    path_dc_input_2,
    path_sample_info_2,
    pickled_results,
):
    # linear with override
    a = b.DriftCorrection(path_dc_input_1, path_sample_info_1)
    a.pool_correct(
        interpolation="linear", pool="PREFA", override=True, max_missing_percent=100
    )
    assert np.isclose(
        a.data.apply(np.floor),
        pickled_results["DCinput1_linear_PREFA_override"].loc[:, a.data.columns],
        equal_nan=True,
    ).all()

    # linear without override
    a = b.DriftCorrection(path_dc_input_1, path_sample_info_1)
    a.pool_correct(
        interpolation="linear", pool="PREFA", override=False, max_missing_percent=100
    )
    assert np.isclose(
        a.data.apply(np.floor),
        pickled_results["DCinput1_linear_PREFA"].loc[:, a.data.columns],
        equal_nan=True,
    ).all()

    # internal standard + NN
    a = b.DriftCorrection(path_dc_input_2, path_sample_info_2)
    a.internal_standard_correct("15R-15-methyl-PGA2")
    a.pool_correct(interpolation="NN", pool="PREFB", max_missing_percent=100)
    assert np.isclose(
        a.data.round(),
        pickled_results["DCinput2_IS_PGA2_NN_PREFB"].loc[:, a.data.columns],
        equal_nan=True,
    ).all()

    # linear with max_missing_percent=30
    a = b.DriftCorrection(path_dc_input_1, path_sample_info_1)
    a.pool_correct(
        interpolation="linear", pool="PREFA", override=True, max_missing_percent=30
    )
    assert np.isclose(
        a.data.apply(np.floor),
        pickled_results["DCinput1_linear_PREFA_override_maxmissing30"]
        .loc[:, a.data.columns]
        .apply(np.floor),
        equal_nan=True,
    ).all()


def test_cv_calculation(
    path_dc_input_2,
    path_sample_info_2,
    pickled_results,
):
    # CV calculation
    a = b.DriftCorrection(path_dc_input_2, path_sample_info_2)
    a.pool_correct(interpolation="linear", pool="PREFA", max_missing_percent=100)
    a.calculate_cvs()
    res = a.cvs.loc[:, ["CV" in col for col in a.cvs.columns]]
    assert np.isclose(
        res.fillna(0),
        pickled_results["DCinput2_linear_PREFA_CVs"].loc[:, res.columns].fillna(0),
    ).all()
