# pylint: disable=redefined-outer-name
"""
Tests for MSEclipse
"""
from pathlib import Path
import pickle
import numpy as np
import pandas as pd
import bmxp.eclipse as ms
import networkx as nx
import pytest


@pytest.fixture()
def filepath_1():
    """
    First file path for test
    """
    return Path(__file__).parent / "test1.csv"


@pytest.fixture()
def filepath_2():
    """
    Second file path for test
    """
    return Path(__file__).parent / "test2.csv"


@pytest.fixture()
def filepath_3():
    """
    Third file path for test
    """
    return Path(__file__).parent / "test3.csv"


@pytest.fixture()
def dataframe_1(filepath_1):
    """
    First dataframe for test
    """
    return pd.read_csv(filepath_1)


@pytest.fixture()
def dataframe_2(filepath_2):
    """
    Second dataframe for test
    """
    return pd.read_csv(filepath_2)


@pytest.fixture()
def dataframe_3(filepath_3):
    """
    Third dataframe for test
    """
    return pd.read_csv(filepath_3)


@pytest.fixture()
def pickled_ms():
    """
    Load pickled file
    """
    return pickle.load(
        open(
            Path(__file__).parent / "mseclipse.pickle", "rb"
        )  # pylint: disable=consider-using-with
    )


@pytest.fixture()
def integration_data():
    """
    Load pickled file
    """
    pickles = []
    datasets = []
    for f in [
        "full_cliques.pickle",
        "ranked_cliques.pickle",
        "all_cliques.pickle",
        "1_non_clique.pickle",
    ]:
        df = pickle.load(
            open(Path(__file__).parent / f, "rb")  # pylint: disable=consider-using-with
        )
        df = df.sort_values(by=["DS1", "DS2", "DS3", "DS4"]).reset_index(drop=True)
        pickles.append(df)
    for f in ["DS1.csv", "DS2.csv", "DS3.csv", "DS4.csv"]:
        datasets.append(pd.read_csv(Path(__file__).parent / f))
    return [pickles, datasets]


def test_initialize_add(filepath_1, filepath_2, dataframe_1, dataframe_2):
    """
    Test MSAligner initializes or throws errors
    """

    # example test: show that importing and utilizing the class is successful
    ms.MSAligner(filepath_1, filepath_2, names=["hp1", "hp2"])

    # Check that dataframes can be read
    ms.MSAligner(dataframe_1, dataframe_2, names=["hp1", "hp2"])

    # datasets require names
    with pytest.raises(ValueError) as e:
        ms.MSAligner(dataframe_1, dataframe_2)
    assert "must be provided" in str(e.value)

    # check mismatch in length
    with pytest.raises(ValueError) as e:
        ms.MSAligner(filepath_1, filepath_2, names=["hp1", "hp2", "hp3"])
    assert "The number of names" in str(e.value)

    # check duplicate names
    with pytest.raises(ValueError) as e:
        ms.MSAligner(filepath_1, filepath_2, names=["hp2", "hp2"])
    assert "duplicate" in str(e.value)


def test_add_dataset(filepath_1, filepath_2, dataframe_1, dataframe_2):
    """
    Test MSAligner add_dataset method
    """
    a = ms.MSAligner()

    # example test: show that importing and utilizing the class is successful
    a.add_dataset(filepath_1, "hp1")
    a.add_dataset(dataframe_2, "hp2")

    # check duplicate names
    with pytest.raises(ValueError) as e:
        a.add_dataset(filepath_2, "hp2")
    assert "duplicate" in str(e.value)

    # test dataframe can be added
    a.add_dataset(dataframe_1, "hp3")

    # dataframes require names
    with pytest.raises(ValueError) as e:
        a.add_dataset(dataframe_2)
    assert "must be provided" in str(e.value)


def test_alignment_methods(dataframe_1, dataframe_2, pickled_ms):
    """
    Test the methods associated with alignment
    """

    a = ms.MSAligner(dataframe_1, dataframe_2, names=["HP1", "HP2"])

    # Test anchors
    a.gen_anchors()
    assert a.anchors["HP1"] == pickled_ms.anchors["HP1"]
    assert a.anchors["HP2"] == pickled_ms.anchors["HP2"]

    # Test coarse matching
    a.gen_coarse()
    assert np.equal(
        a.coarse_matches["HP1"]["HP2"].values,
        pickled_ms.coarse_matches["HP1"]["HP2"].values,
    ).all()

    # Test scaler generation; RT, MZ, and intensity
    a.gen_scalers()
    hp1_scalers = a.scalers["HP1"]["HP2"]
    pickled_scalers = pickled_ms.scalers["HP1"]["HP2"]
    assert np.isclose(hp1_scalers["RT"]["y"], pickled_scalers["RT"]["y"]).all()
    assert np.isclose(hp1_scalers["RT"]["x"], pickled_scalers["RT"]["x"]).all()
    assert np.isclose(hp1_scalers["MZ"]["y"], pickled_scalers["MZ"]["y"]).all()
    assert np.isclose(hp1_scalers["MZ"]["x"], pickled_scalers["MZ"]["x"]).all()
    assert np.isclose(
        hp1_scalers["Intensity"]["x"], pickled_scalers["Intensity"]["x"]
    ).all()
    assert np.isclose(
        hp1_scalers["Intensity"]["y"], pickled_scalers["Intensity"]["y"]
    ).all()

    # Test scaling functions
    a.gen_scaled_values()
    scaled_values = a.scaled_values["HP1"]["HP2"]
    pickled_values = pickled_ms.scaled_values["HP1"]["HP2"]
    assert np.isclose(scaled_values["RT"], pickled_values["RT"]).all()
    assert np.isclose(scaled_values["MZ"], pickled_values["MZ"]).all()
    assert np.isclose(scaled_values["Intensity"], pickled_values["Intensity"]).all()

    # Test standard deviation calculations
    a.gen_stds()
    assert np.isclose(a.stds["HP1"]["HP2"]["RT"], pickled_ms.stds["HP1"]["HP2"]["RT"])
    assert np.isclose(a.stds["HP1"]["HP2"]["MZ"], pickled_ms.stds["HP1"]["HP2"]["MZ"])
    assert np.isclose(
        a.stds["HP1"]["HP2"]["Intensity"], pickled_ms.stds["HP1"]["HP2"]["Intensity"]
    )

    # Test matches
    a.gen_matches()
    assert a.matches["HP1"]["HP2"].equals(pickled_ms.matches["HP1"]["HP2"])

    # Check that align passes
    a = ms.MSAligner(dataframe_1, dataframe_2, names=["HP1", "HP2"])
    a.align()
    a.results()

    # Check that custom scalers are not overwritten
    rt = ms.calc_scalers([1, 2, 3], [1.1, 2.2, 3.3])
    mz = ms.calc_scalers([100, 200, 300], [100.0005, 200.0008, 300.0009])
    intensity = ms.calc_scalers([2, 3, 4], [3, 4, 5])

    a = ms.MSAligner(dataframe_1, dataframe_2, names=["HP1", "HP2"])
    a.set_scalers(
        {
            "HP1": {
                "HP2": {
                    "RT": zip(*rt),
                    "MZ": zip(*mz),
                    "Intensity": zip(*intensity),
                }
            }
        },
        rec=False,
    )
    a.align()
    a_scalers = a.scalers["HP1"]["HP2"]
    assert np.equal(a_scalers["RT"]["x"], rt[0]).all()
    assert np.equal(a_scalers["RT"]["y"], rt[1]).all()
    assert np.equal(a_scalers["MZ"]["x"], mz[0]).all()
    assert np.equal(a_scalers["MZ"]["y"], mz[1]).all()
    assert np.equal(a_scalers["Intensity"]["x"], intensity[0]).all()
    assert np.equal(a_scalers["Intensity"]["y"], intensity[1]).all()

    # check reciprocity for scaler generation works
    a = ms.MSAligner(dataframe_1, dataframe_2, names=["HP1", "HP2"])
    a.set_scalers(
        {
            "HP1": {
                "HP2": {
                    "RT": zip(*rt),
                    "MZ": zip(*mz),
                    "Intensity": zip(*intensity),
                }
            }
        },
        rec=True,
    )
    a.align()
    a_scalers = a.scalers["HP2"]["HP1"]
    assert np.equal(a_scalers["RT"]["x"], rt[0] + rt[1]).all()
    assert np.equal(a_scalers["RT"]["y"], -rt[1]).all()
    new_mzs = mz[0] + mz[0] * mz[1] / 1000000
    assert np.equal(a_scalers["MZ"]["x"], new_mzs).all()
    assert np.equal(a_scalers["MZ"]["y"], -mz[1] * mz[0] / new_mzs).all()
    assert np.equal(a_scalers["Intensity"]["x"], intensity[0] + intensity[1]).all()
    assert np.equal(a_scalers["Intensity"]["y"], -intensity[1]).all()


def test_calc_score():
    """
    Test score calculation
    """
    score = ms.calc_score(
        {"RT": 1, "MZ": 200, "Intensity": 100},
        {"RT": 2, "MZ": 200.0002, "Intensity": 1000},
        {"RT": "linear", "MZ": "ppm", "Intensity": "log10"},
        {"RT": 1, "MZ": 1, "Intensity": 1},
        {"RT": 1, "MZ": 1, "Intensity": 1},
    )
    assert np.isclose(score, 3)


@pytest.mark.filterwarnings("ignore:invalid value")
def test_calc_scalers():
    """
    Supplemental tests to calculating scalers
    """

    # Will return either nans, the original y_vals, or 1s.
    # Test it does not make an infinite loop
    x_vals = np.array([0, 0, 1])
    y_vals = np.array([1, 1, 2])
    results = ms.calc_scalers(x_vals, y_vals, smoothing="lowess", mode="linear")[1]
    all_nans = np.isnan(results).all()
    is_original = np.isclose(np.array([1, 1, 1]), results).all()
    assert all_nans or is_original


def test_adding_intensity_column(dataframe_3):
    """
    Tests to see if intensity is autocalculated when missing
    """
    # Create Dataframes without Intensity
    df_3_no_intensity = dataframe_3.drop(columns="Intensity")

    a = ms.MSAligner()
    a.add_dataset(df_3_no_intensity, "df3")
    # Test that _calc_intensity functions on its own
    assert np.isclose(a.datasets["df3"]["Intensity"], dataframe_3["Intensity"]).all()


def test_dataset_selection(dataframe_1, dataframe_2, dataframe_3):
    """
    Tests the wrappers which control the inner and outer dataframe for loops
    """
    # test all works
    names = ["HP1", "HP2", "HP3"]
    a = ms.MSAligner(dataframe_1, dataframe_2, dataframe_3, names=names)
    a.align()
    assert list(a.matches) == ["HP1", "HP2", "HP3"]
    assert list(a.matches["HP1"]) == ["HP2", "HP3"]
    assert list(a.matches["HP2"]) == ["HP1", "HP3"]
    assert list(a.matches["HP3"]) == ["HP1", "HP2"]

    # test list to all
    a = ms.MSAligner(dataframe_1, dataframe_2, dataframe_3, names=names)
    a.align(["HP1", "HP2", "HP3"])
    assert list(a.matches) == ["HP1", "HP2", "HP3"]
    assert list(a.matches["HP1"]) == ["HP2", "HP3"]
    assert list(a.matches["HP2"]) == ["HP1", "HP3"]
    assert list(a.matches["HP3"]) == ["HP1", "HP2"]

    # test one way to all
    a = ms.MSAligner(dataframe_1, dataframe_2, dataframe_3, names=names)
    a.align("HP1")
    assert list(a.matches) == ["HP1"]
    assert list(a.matches["HP1"]) == ["HP2", "HP3"]

    # test one to one
    a = ms.MSAligner(dataframe_1, dataframe_2, dataframe_3, names=names)
    a.align("HP1", "HP2")
    assert list(a.matches) == ["HP1"]
    assert list(a.matches["HP1"]) == ["HP2"]


def test_schema(dataframe_1, dataframe_2):
    dataframe_1 = dataframe_1.copy()
    dataframe_2 = dataframe_2.copy()
    schema_labels = {
        "Compound_ID": "feature_id",
        "RT": "rt (min)",
        "MZ": "m/z",
        "Intensity": "abundance",
    }

    a = ms.MSAligner(
        dataframe_1.rename(columns=schema_labels),
        dataframe_2.rename(columns=schema_labels),
        names=["DS1", "DS2"],
        schema_labels=schema_labels,
    )
    a.align()


def test_instance_params(dataframe_1, dataframe_2, dataframe_3):
    schema_labels = {"RT": "rt (min)"}
    dataframe_1 = dataframe_1.copy().rename(columns=schema_labels)
    dataframe_2 = dataframe_2.copy().rename(columns=schema_labels)
    dataframe_3 = dataframe_3.copy().rename(columns=schema_labels)

    a = ms.MSAligner(
        dataframe_1,
        dataframe_2,
        dataframe_3,
        names=["DS1", "DS2", "DS3"],
        schema_labels=schema_labels,
    )
    a.descriptors = {"rt (min)": "linear", "MZ": "ppm"}
    a.default_coarse_params["rt (min)"] = {"upper": 1.0, "lower": -1.0}
    a.default_cutoff = 10
    a.default_cutoffs = {"MZ": 6}
    a.align()

    a = ms.MSAligner(
        dataframe_1,
        dataframe_2,
        dataframe_3,
        names=["DS1", "DS2", "DS3"],
        schema_labels=schema_labels,
    )
    a.descriptors = {"rt (min)": "linear", "MZ": "ppm"}
    a.set_instance_defaults(
        {
            "coarse_params": {"rt (min)": {"upper": 1.0, "lower": -1.0}},
            "cutoff": 10,
            "cutoffs": {"MZ": 6},
        }
    )
    a.align()


@pytest.fixture()
def graph():
    """
    Digraph for testing compression
    """
    graph = nx.DiGraph()
    nodes = ["1__a", "2__a", "3__a", "4__a", "1__b", "2__b", "3__b", "4__b", "1__c"]
    graph.add_nodes_from(nodes)
    edges = [
        ("1__a", "2__a", {"score": 0.5, "rank": 0}),  # kept
        ("1__a", "3__a", {"score": 0.9, "rank": 0}),  # kept
        ("1__a", "4__a", {"score": 2, "rank": 0}),  # kept
        ("2__a", "1__a", {"score": 0.3, "rank": 0}),  # kept
        ("2__a", "3__a", {"score": 0.1, "rank": 0}),  # kept
        ("3__a", "1__a", {"score": 0.6, "rank": 0}),  # kept
        ("3__a", "2__a", {"score": 0.1, "rank": 0}),  # kept
        ("3__a", "4__a", {"score": 100, "rank": 0}),
        ("1__b", "2__a", {"score": 100, "rank": 0}),
        ("1__b", "3__b", {"score": 10, "rank": 0}),  # kept
        ("2__b", "1__b", {"score": 100, "rank": 0}),
        ("2__b", "3__a", {"score": 100, "rank": 0}),
        ("3__b", "1__b", {"score": 11, "rank": 0}),  # kept
        ("3__b", "2__b", {"score": 100, "rank": 0}),
        ("4__b", "2__b", {"score": 100, "rank": 0}),
        ("4__a", "1__a", {"score": 1, "rank": 0}),
        ("1__c", "4__a", {"score": 100, "rank": 0}),
    ]
    graph.add_edges_from(edges)
    return graph


def test_compress_graph(graph):
    """
    Test graph compression
    """

    compressed = ms.MSAligner.compress_graph(graph)

    # check nodes are present
    assert set(compressed.nodes) == set(graph.nodes)

    # check edges and weights
    edges_weights = set(
        frozenset([e[0], e[1], e[2]["score"]]) for e in compressed.edges(data=True)
    )
    assert edges_weights == {
        frozenset(["1__a", "2__a", 0.3 + 0.5]),
        frozenset(["2__a", "3__a", 0.1 + 0.1]),
        frozenset(["3__a", "1__a", 0.6 + 0.9]),
        frozenset(["4__a", "1__a", 1 + 2]),
        frozenset(["1__b", "3__b", 10 + 11]),
    }


@pytest.fixture()
def complex_compressed():
    """
    A complex feature graph
    """
    graph = nx.Graph()
    nodes = [
        ("1__a", {"dataset": 1}),
        ("2__a", {"dataset": 2}),
        ("3__a", {"dataset": 3}),
        ("4__a", {"dataset": 4}),
        ("5__a", {"dataset": 5}),
        ("1__b", {"dataset": 1}),
        ("2__b", {"dataset": 2}),
        ("3__b", {"dataset": 3}),
        ("4__b", {"dataset": 4}),
        ("1__c", {"dataset": 1}),
        ("2__c", {"dataset": 2}),
        ("3__d", {"dataset": 3}),
        ("4__c", {"dataset": 4}),
    ]
    graph.add_nodes_from(nodes)

    edges = [
        ("1__a", "2__a", {"score": 1}),
        ("3__a", "2__a", {"score": 1}),
        ("1__a", "3__a", {"score": 1}),
        ("4__a", "5__a", {"score": 0.7}),
        ("1__a", "5__a", {"score": 0.5}),
        ("1__a", "4__a", {"score": 0.3}),
        ("4__a", "2__b", {"score": 2}),
        ("4__a", "3__b", {"score": 2}),
        ("2__b", "3__b", {"score": 2}),
        ("4__b", "1__b", {"score": 0.5}),
        ("2__c", "1__b", {"score": 1}),
        ("3__d", "1__b", {"score": 1}),
        ("2__c", "3__d", {"score": 5}),
        ("3__b", "1__c", {"score": 1}),
        ("4__c", "1__c", {"score": 1}),
        ("5__b", "3__d", {"score": 0.1}),
        ("3__d", "4__b", {"score": 0.1}),
        ("2__c", "4__b", {"score": 0.1}),
        ("2__c", "5__b", {"score": 0.1}),
        ("1__c", "2__b", {"score": 0.1}),
    ]
    graph.add_edges_from(edges)
    return graph


@pytest.fixture()
def clique_compressed():
    """
    A simple clique graph
    """
    graph = nx.Graph()
    nodes = [
        ("1__a", {"dataset": 1}),
        ("2__a", {"dataset": 2}),
        ("3__a", {"dataset": 3}),
    ]
    graph.add_nodes_from(nodes)

    edges = [
        ("1__a", "2__a", {"score": 1}),
        ("3__a", "2__a", {"score": 1}),
        ("1__a", "3__a", {"score": 1}),
    ]
    graph.add_edges_from(edges)
    return graph


def test_deconvolution(complex_compressed, clique_compressed):
    """
    Test the deconvolution function
    """
    max_size = 3
    # test it returns a full clique
    result = ms.MSAligner.get_rank_groups(clique_compressed, max_size)
    assert len(result) == 1

    # test that it only reports 1 group
    result = ms.MSAligner.get_rank_groups(clique_compressed, max_size, 1, 1, 2)
    assert len(result) == 1

    max_size = 5
    # max cliques, 4, 3, >=1
    result = ms.MSAligner.get_rank_groups(complex_compressed, max_size)
    assert len(result) == 0
    result = ms.MSAligner.get_rank_groups(complex_compressed, max_size, 4, 4, 1)
    assert len(result) == 1
    result = ms.MSAligner.get_rank_groups(complex_compressed, max_size, 3, 3, 1)
    assert len(result) == 6
    result = ms.MSAligner.get_rank_groups(complex_compressed, max_size, 1, 1, 1)
    assert len(result) == 7

    # test the settings return the correct, sorted cliques
    expected = [
        ({"5__b", "2__c", "4__b", "3__d", "1__b"}, 5, 4, 8, 0.9874999999999998),
        ({"1__a", "3__a", "5__a", "2__a", "4__a"}, 5, 3, 6, 0.75),
        ({"1__a", "5__a", "3__b", "4__a", "2__b"}, 5, 3, 6, 1.25),
        ({"1__c", "2__b", "4__a", "3__b"}, 4, 3, 5, 1.42),
        ({"1__c", "4__c", "3__b", "2__b"}, 4, 3, 4, 1.025),
        ({"3__b", "1__b", "4__a", "2__b"}, 4, 3, 4, 1.75),
    ]
    results = ms.MSAligner.get_rank_groups(complex_compressed, max_size, 1, 1, 2)
    for i, r in enumerate(results):
        for x, y in zip(r, expected[i]):
            if isinstance(y, (float, int)):
                assert np.isclose(x, y)
            else:
                assert x == y


def test_integration(integration_data):
    """
    Tests alignment and 4 different aggregation methods
    """
    pickles = integration_data[0]
    datasets = integration_data[1]

    a = ms.MSAligner(*datasets, names=["DS1", "DS2", "DS3", "DS4"])
    a.align()

    results = a.results().loc[:, "DS1":"DS4"]
    results = results.sort_values(by=["DS1", "DS2", "DS3", "DS4"]).reset_index(
        drop=True
    )
    pickles[0].equals(results)

    results = a.results(c_size_or_loss=1, g_size_or_loss=1).loc[:, "DS1":"DS4"]
    results = results.sort_values(by=["DS1", "DS2", "DS3", "DS4"]).reset_index(
        drop=True
    )
    pickles[1].equals(results)

    results = a.results(c_size_or_loss=1, g_size_or_loss=1, remove_rerank=False).loc[
        :, "DS1":"DS4"
    ]
    results = results.sort_values(by=["DS1", "DS2", "DS3", "DS4"]).reset_index(
        drop=True
    )
    pickles[2].equals(results)

    results = a.results(c_size_or_loss=3, g_size_or_loss=3, diameter=2).loc[
        :, "DS1":"DS4"
    ]
    results = results.sort_values(by=["DS1", "DS2", "DS3", "DS4"]).reset_index(
        drop=True
    )
    pickles[3].equals(results)
