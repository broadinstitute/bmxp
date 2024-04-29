# pylint: disable=redefined-outer-name
"""
Tests for Gravity
"""
from pathlib import Path
import pytest
import pandas as pd
from bmxp import gravity


@pytest.fixture()
def filepath():
    """
    First file path for test
    """
    return Path(__file__).parent / "test_gravity.csv"


@pytest.fixture()
def dataframe(filepath):
    """
    First dataframe for test
    """
    return pd.read_csv(filepath, header=0, index_col="Compound_ID")


def test_initialize_add(dataframe):
    """
    Test MSAligner initializes or throws errors
    """
    gravity.cluster(dataframe)
