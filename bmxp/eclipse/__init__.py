"""
MSEclipse is a LCMS dataset alignment package, designed to combine datasets which have
been collected on the same method but perhaps on different instruments or columns.
"""

# pylint: disable=too-many-lines
from itertools import chain, combinations
import copy
import logging
import collections
import io
from typing import Optional, List, Union
from itertools import product
import numpy as np
import pandas as pd
import networkx as nx
import textwrap
from scipy import interpolate
import statsmodels.api as sm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from bmxp import FMDATA

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)

np.random.seed(0)

lowess = sm.nonparametric.lowess
__version__ = "0.2.7"


def dataset_loops(attr=None):
    """
    Wraps a function with a double For-loop for iterating over the datasets.
    It converts the supplied datasets, which may be None, 1 dataset, or multiple,
    to a list.
    Also handles creation of a dict key for an attribute and makes it if it doesn't
     exist. E.g. self.scalers[DS1] -- if that key doesn't exist, it will create it
     so that self.scalers[DS1][DS2] can be assigned.
    """

    def attr_handler(func):
        def loops(self, datasets_1=None, datasets_2=None, *args, **kwargs):
            # convert datasets to a list of datasets for the loop
            if datasets_1 is None:
                datasets_1 = list(self.datasets)
            if isinstance(datasets_1, str):
                datasets_1 = [datasets_1]

            if datasets_2 is None:
                datasets_2 = list(self.datasets)
            if isinstance(datasets_2, str):
                datasets_2 = [datasets_2]

            for ds1 in datasets_1:
                if attr == "coarse_matches":
                    LOGGER.info(f"Generating coarse matches for {ds1} to others...")
                if attr == "scalers":
                    LOGGER.info(f"Generating scalers for {ds1} to others...")
                if attr == "scaled_values":
                    LOGGER.info(f"Scaling {ds1} to others...")
                if attr == "stds":
                    LOGGER.info(f"Calculating variance from {ds1} to others...")
                if attr == "matches":
                    LOGGER.info(f"Generating {ds1} feature tables...")
                if attr is None:
                    LOGGER.info(f"Loading feature graph with {ds1} feature tables...")

                if attr and ds1 not in self.__getattribute__(attr):
                    self.__getattribute__(attr)[ds1] = {}
                for ds2 in datasets_2:
                    if ds1 == ds2:
                        continue
                    func(self, ds1, ds2, *args, **kwargs)

        return loops

    return attr_handler


def ds2_loop(func):
    """
    Wraps a function with an inner ds2 loop. Does not check for attribute
     dictionary keys (handled in functions)
    """

    def wrapper(self, ds1=None, ds2=None, *args, **kwargs):
        ds2 = list(self.datasets) if ds2 is None else ds2
        if not isinstance(ds2, str):
            for item in ds2:
                if item == ds1:
                    continue
                func(self, ds1, item, *args, **kwargs)
        else:
            func(self, ds1, ds2, *args, **kwargs)

    return wrapper


class MSAligner:  # pylint: disable=too-many-instance-attributes
    """
    MSAligner stores 2 or more MS Datasets and performs matching, scaling, and
    alignments between them.
    """

    # references FMDATA keys, not schema values
    # values are evaluated in __init__
    descriptors = {"RT": "linear", "MZ": "ppm", "Intensity": "log10"}
    mz_col = "MZ"
    default_cutoff = 6
    default_weight = 1
    default_cutoffs = {}
    default_weights = {}
    default_coarse_params = {
        "RT": {"upper": 0.5, "lower": -0.5},
        "MZ": {"upper": 15, "lower": -15},
        "Intensity": {"upper": 2, "lower": -2},
    }
    default_scaler_params = {
        "smoothing_method": "lowess",
        "smoothing_params": {"frac": 0.1},
    }
    feature_name = "Compound_ID"
    intensity_col = "Intensity"
    annotation_col = "Annotation_ID"

    def __init__(
        self,
        *args: Union[str, pd.DataFrame],
        names: Optional[List[str]] = None,
        schema_labels: Optional[dict[str, str]] = None,
    ):
        fmdata = FMDATA.copy()

        self.mz_col = (schema_labels or {}).get(
            MSAligner.mz_col, FMDATA[MSAligner.mz_col]
        )

        if schema_labels:
            fmdata.update(schema_labels)

        self.descriptors = {fmdata[k]: v for k, v in MSAligner.descriptors.items()}
        self.default_cutoff = MSAligner.default_cutoff
        self.default_weight = MSAligner.default_weight
        self.default_weights = {
            fmdata[k]: v for k, v in MSAligner.default_weights.items()
        }
        self.default_cutoffs = {
            fmdata[k]: v for k, v in MSAligner.default_cutoffs.items()
        }

        # for k in MSAligner.descriptors:
        #     if k not in MSAligner.default_weights:
        if MSAligner.default_weights is None:
            self.default_weights = {
                fmdata[k]: self.default_weight for k in MSAligner.descriptors
            }

        self.default_coarse_params = {
            fmdata[k]: v for k, v in MSAligner.default_coarse_params.items()
        }
        self.default_scaler_params = MSAligner.default_scaler_params
        self.feature_name = fmdata[MSAligner.feature_name]
        self.intensity_col = fmdata[MSAligner.intensity_col]
        self.annotation_col = fmdata[MSAligner.annotation_col]
        self.datasets = collections.OrderedDict()
        self.anchors = {}
        self.coarse_matches = {}
        self.scalers = {}
        self.scaled_values = {}
        self.stds = {}
        self.matches = {}
        self.scores = {}
        self.cutoffs = {}
        self.weights = {}
        self.coarse_params = {}
        self.scaler_params = {}
        self.multipliers = {}
        self.prescalers = {}
        self.graph = nx.DiGraph()
        self.remove_all = True
        # convert names to a list if not supplied
        if names is None:
            names = [None] * len(args)

        # Check there are the same number of datasets as names
        if len(names) != len(args):
            raise ValueError(
                "The number of names you provided does not match the "
                "number of datasets."
            )

        # Add datasets using the add_dataset method
        for arg, name in zip(args, names):
            self.add_dataset(arg, name)

    def add_dataset(self, dataset, name=None):
        """
        Adds a dataset to the existing alignment object.

        :param dataset: str, filepath of .csv dataset
        :param name: str, name of dataset
        """

        if name is None:
            # check if there are any dataframes present; if so, raise an exception
            if isinstance(dataset, pd.DataFrame):
                raise ValueError(
                    "Names must be provided if importing DataFrames directly."
                )
            name = dataset

        if "__" in name:
            raise ValueError(
                "Your dataset name (or filename) cannot contain double underscores."
            )

        if not isinstance(dataset, pd.DataFrame):
            dataset = pd.read_csv(dataset)

        # check that the name isn't already present
        if name in self.datasets.keys():
            raise ValueError("You have duplicate names for your datasets.")

        self.datasets[name] = dataset
        self.datasets[name].index = self.datasets[name].index.map(str)
        self._validate_data(name)

    def _validate_data(self, names=None):
        """
        Performs a check to see if a dataset has any compatibility errors.

        :param names: list, name of all datasets to check
        """

        # If names are not provided, check all
        if names is None:
            names = self.datasets.keys()

        # If one name is provided as a string, convert to a 1 member list
        elif isinstance(names, str):
            names = [names]

        # Check all datasets included in the names
        for name in names:
            ds = self.datasets[name]
            # calculate intensity if needed
            if (
                self.intensity_col not in ds.columns
                and self.intensity_col in self.descriptors
            ):
                gen_intensity(
                    ds,
                    self.intensity_col,
                    ignore=[
                        descr
                        for descr in self.descriptors
                        if descr != self.intensity_col
                    ]
                    + [self.feature_name],
                )

            # Check required columns
            required_columns = [*self.descriptors, self.feature_name]
            for req_col in required_columns:
                if req_col not in ds.columns:
                    raise ValueError(f"Column {req_col} not found in dataset {name}.")

    def set_instance_defaults(self, params):
        """
        Updates default values for weights and cutoffs.
        """

        for key, item in params.items():
            if key.lower() == "weight":
                self.default_weight = item
            elif key.lower() == "cutoff":
                self.default_cutoff = item
            elif key.lower() == "weights":
                self.default_weights.update(item)
            elif key.lower() == "cutoffs":
                self.default_cutoffs.update(item)
            elif key.lower() == "coarse_params":
                self.default_coarse_params.update(item)
            elif key.lower() == "scaler_params":
                self.default_scaler_params.update(item)
            else:
                raise KeyError(
                    f"Unknown value: {key}. Acceptable values are "
                    " 'weight', 'weights', 'cutoff', 'cutoffs',"
                    " 'coarse_params' and 'scaler_params'."
                )

    def set_params(self, params, rec=True):
        """
        Sets weights, cutoffs, and scalers for any dataset pairs
        For example:
        {"scalers": {
            "DS1": {"DS2": {"RT": 1}}
            },
        "cutoffs": { ...}
        }
        """
        for key, item in params.items():
            if key.lower() == "scalers":
                self.set_scalers(item, rec)
            elif key.lower() == "cutoffs":
                self.set_cutoffs(item, rec)
            elif key.lower() == "weights":
                self.set_weights(item, rec)
            elif key.lower() == "coarse_params":
                self.set_coarse_params(item, rec)
            elif key.lower() == "scaler_params":
                self.set_scaler_params(item, rec)
            else:
                raise KeyError(
                    f"Unknown value: {key}. Acceptable values are 'scalers', "
                    "'cutoffs', 'weights', 'coarse_params', and 'scaler_params'."
                )

    def set_scalers(self, scalers, rec=True):
        """
        Function to set custom scalers.
         If rec is true, will calculate the reciprocal scalers
         scalers is a nested dict with the dataframe names
         and scalers to be used. For example:
        {
            'DS1': {
                'DS2':{
                    'RT': {'x': [1, 2, 3], 'y': [.1, .2, .3]},
                    'MZ': {'x': [100, 200, 300], 'y' : [2, 4, 3]},
                    'Intensity': 1.001
                }
            }
        }
        While it seems clunky, it allows us to set multiple scalers at once.
         If rec is True, and recipricol scalers are specified, they will be overwritten.
        :param scalers: dict
        :param rec: bool
        """
        for ds1 in scalers:
            # make scaler dicts if they don't exist
            if ds1 not in self.scalers:
                self.scalers[ds1] = {}
            for ds2 in scalers[ds1]:
                if ds2 not in self.scalers[ds1]:
                    self.scalers[ds1][ds2] = {}
                if rec:
                    if ds2 not in self.scalers:
                        self.scalers[ds2] = {}
                    if ds1 not in self.scalers[ds2]:
                        self.scalers[ds2][ds1] = {}

                # Check if dataframes are in the aligner object
                if ds1 not in self.datasets or ds2 not in self.datasets:
                    raise IndexError(
                        f"Either {ds1} or {ds2} was not found in the datasets."
                    )

                # Check for typos in scaler names
                if not set(scalers[ds1][ds2].keys()).issubset(
                    set(self.descriptors.keys())
                ):
                    raise IndexError(
                        "The only acceptable columns for scalers are "
                        f"{self.descriptors.keys()}. "
                        f"You provided {list(scalers.keys())}."
                    )

                self._set_scalers(scalers[ds1][ds2], ds1, ds2, rec)

    def _set_scalers(self, scalers, ds1, ds2, rec):
        """
        Takes a list of scalers in a dictionary form
         and sets them as scalers in the object.
        For example:
        {
            'RT': {'x': [1, 2, 3], 'y': [.1, .2, .3]},
            'MZ': {'x': [100, 200, 300], 'y' : [2, 4, 3]},
            'Intensity': {'x': [3, 4, 5], 'y' : [.1, .3, .4]}
        }
        Scalers should already be in their proper space

        :param scalers: Dict(str: Dict(str:list),
         containing the scaler names as Keys and values
        :param ds1: str, name of Source dataset
        :param ds2: str, name of Target dataset
        :param rec: bool, Flag to create reciprocal or not (DS1->DS2, DS2->DS1)
            Convenience for manually specified scalers, not for use autogenerated
        """
        for k, v in scalers.items():
            values = list(zip(*v))
            df = pd.DataFrame({"x": values[0], "y": values[1]})
            self.scalers[ds1][ds2][k] = df
            if rec:
                if self.descriptors[k] in ["linear", "log10"]:
                    self.scalers[ds2][ds1][k] = pd.DataFrame(
                        {"x": df["x"] + df["y"], "y": -df["y"]}
                    )
                elif self.descriptors[k] == "ppm":
                    new_xs = df["x"] + df["x"] * df["y"] / 1000000
                    new_ys = -df["y"] * df["x"] / new_xs
                    self.scalers[ds2][ds1][k] = pd.DataFrame({"x": new_xs, "y": new_ys})
                else:
                    raise NotImplementedError(
                        "Not sure how to reverse scale the mode ", self.descriptors[k]
                    )

    def set_prescalers(self, prescalers):
        """
        Dict containing datasets and descriptors providing the dataset and
        the instructions to be applied.
        Five spread across the gradient should be plenty.
        This will shift DS1 RTs earlier, since the instructions are negative
        {
            'DS1': {
                'RT': [[1, 5, 10], [0, -0.2, -0.3]]
            }
        }
        """
        self.prescalers = {}
        for dataset in prescalers:
            self.prescalers[dataset] = {}
            for descriptor in prescalers[dataset]:
                if descriptor not in self.descriptors:
                    raise KeyError(
                        f"Your prescaler descriptor {descriptor} is not present as a "
                        "descriptor."
                    )
                self.prescalers[dataset][descriptor] = [
                    pd.Series(prescalers[dataset][descriptor][0]),
                    pd.Series(prescalers[dataset][descriptor][1]),
                ]

    def set_cutoffs(self, cutoffs, rec=True):
        """
        Function to set custom cutoffs. If rec is true, apply reciprocal cutoffs.
        {
            'DS1': {
                'DS2':{
                    'RT': 5.5,
                    'MZ': 5.5,
                    'Intensity': 5.0
                }
            }
        }
        :param cutoffs: dict
        :param rec: bool
        """
        cutoffs = cutoffs.copy()
        for ds1 in cutoffs:
            # make cutoff dicts if they don't exist
            if ds1 not in self.cutoffs:
                self.cutoffs[ds1] = {}
            for ds2 in cutoffs[ds1]:
                if ds2 not in self.cutoffs[ds1]:
                    self.cutoffs[ds1][ds2] = {}

                # Check if dataframes are in the aligner object
                if ds1 not in self.datasets or ds2 not in self.datasets:
                    raise IndexError(
                        f"Either {ds1} or {ds2} was not found in the datasets."
                    )

                # Check for typos in cutoff names
                if not set(cutoffs[ds1][ds2].keys()).issubset(set(self.descriptors)):
                    raise IndexError(
                        "The only acceptable columns for cutoffs are "
                        f"{self.descriptors.keys()}. You provided "
                        f"{cutoffs[ds1][ds2].keys()}."
                    )

                # replace 0's with inf, negating their use in cutoffs
                for key, value in cutoffs[ds1][ds2].items():
                    if value == 0:
                        cutoffs[ds1][ds2][key] = float("inf")

                self.cutoffs[ds1][ds2].update(cutoffs[ds1][ds2])
                if rec:
                    if ds2 not in self.cutoffs:
                        self.cutoffs[ds2] = {}
                    if ds1 not in self.cutoffs[ds2]:
                        self.cutoffs[ds2][ds1] = {}
                    self.cutoffs[ds2][ds1].update(cutoffs[ds1][ds2])

    def set_weights(self, weights, rec=True):
        """
        Function to set custom weights. If rec is true, apply reciprocal weights
        {
            'DS1': {
                'DS2':{
                    'RT': .5,
                    'MZ': 1,
                    'Intensity': 0
                }
            }
        }
        :param weights: dict
        :param rec: bool
        """
        weights = weights.copy()
        for ds1 in weights:
            # make weight dicts if they don't exist
            if ds1 not in self.weights:
                self.weights[ds1] = {}
            for ds2 in weights[ds1]:
                if ds2 not in self.weights[ds1]:
                    self.weights[ds1][ds2] = {}

                # Check if dataframes are in the aligner object
                if ds1 not in self.datasets or ds2 not in self.datasets:
                    raise IndexError(
                        f"Either {ds1} or {ds2} was not found in the datasets."
                    )

                # Check for typos in weight names
                if not set(weights[ds1][ds2].keys()).issubset(set(self.descriptors)):
                    raise IndexError(
                        "The only acceptable columns for weights are "
                        f"{self.descriptors.keys()}.  You provided "
                        f"{weights[ds1][ds2].keys()}."
                    )

                self.weights[ds1][ds2].update(weights[ds1][ds2])
                if rec:
                    if ds2 not in self.weights:
                        self.weights[ds2] = {}
                    if ds1 not in self.weights[ds2]:
                        self.weights[ds2][ds1] = {}
                    self.weights[ds2][ds1].update(weights[ds1][ds2])

    def set_coarse_params(
        self, coarse_params, rec=True
    ):  # pylint: disable=too-many-branches
        """
        Function to set custom cutoffs. If rec is true, apply reciprocal cutoffs
        {
            'DS1': {
                'DS2':{
                    'rt_minus': -.7,
                    'ppm_plus': 4
                }
            }
        }
        :param coarse_params: dict
        :param rec: bool
        """
        coarse_params = coarse_params.copy()
        for ds1 in coarse_params:
            # make coarse_params dicts if they don't exist
            if ds1 not in self.coarse_params:
                self.coarse_params[ds1] = {}
            for ds2 in coarse_params[ds1]:
                if ds2 not in self.coarse_params[ds1]:
                    self.coarse_params[ds1][ds2] = {}

                # Check if dataframes are in the aligner object
                if ds1 not in self.datasets or ds2 not in self.datasets:
                    raise IndexError(
                        f"Either {ds1} or {ds2} was not found in the datasets."
                    )

                # Check for typos in coarse_param names
                if not set(coarse_params[ds1][ds2].keys()).issubset(
                    set(key for key in self.descriptors)
                ):
                    raise IndexError(
                        f"The only acceptable columns for coarse_params are"
                        f" {set(key for key in self.descriptors)}. You"
                        f" provided {coarse_params[ds1][ds2].keys()}."
                    )

                self.coarse_params[ds1][ds2].update(coarse_params[ds1][ds2])
                if rec:
                    rec_dict = {}
                    for descriptor in coarse_params[ds1][ds2]:
                        new_feature = {}
                        if "upper" in coarse_params[ds1][ds2][descriptor]:
                            new_feature["lower"] = -coarse_params[ds1][ds2][descriptor][
                                "upper"
                            ]
                        if "lower" in coarse_params[ds1][ds2][descriptor]:
                            new_feature["upper"] = -coarse_params[ds1][ds2][descriptor][
                                "lower"
                            ]
                        rec_dict.update({descriptor: new_feature})

                    if ds2 not in self.coarse_params:
                        self.coarse_params[ds2] = {}
                    if ds1 not in self.coarse_params[ds2]:
                        self.coarse_params[ds2][ds1] = {}
                    self.coarse_params[ds2][ds1].update(rec_dict)

    def set_scaler_params(self, scaler_params, rec=True):
        """
        Function to set custom smoothing params for scaler generation. If rec is true,
        apply reciprocal params
        {
            'DS1': {
                'DS2':{
                    'smoothing_method': 'lowess',
                    'smoothing_params': {
                        'frac': 0.9
                    }
                }
            }
        }
        :param scaler_params: dict
        :param rec: bool
        """
        scaler_params = scaler_params.copy()
        for ds1 in scaler_params:
            # make scaler_params dicts if they don't exist
            if ds1 not in self.scaler_params:
                self.scaler_params[ds1] = {}
            for ds2 in scaler_params[ds1]:
                if ds2 not in self.scaler_params[ds1]:
                    self.scaler_params[ds1][ds2] = {}

                # Check if dataframes are in the aligner object
                if ds1 not in self.datasets or ds2 not in self.datasets:
                    raise IndexError(
                        f"Either {ds1} or {ds2} was not found in the datasets."
                    )

                self.scaler_params[ds1][ds2].update(scaler_params[ds1][ds2])
                if rec:
                    if ds2 not in self.scaler_params:
                        self.scaler_params[ds2] = {}
                    if ds1 not in self.scaler_params[ds2]:
                        self.scaler_params[ds2][ds1] = {}
                    self.scaler_params[ds2][ds1].update(scaler_params[ds1][ds2])

    def to_csv(
        self,
        *args,
        filepath=None,
        to_bytes=False,
        c_size_or_loss=0,
        g_size_or_loss=0,
        diameter=1,
        remove_rerank=True,
    ):
        """
        Condenses network into a flat CSV file
        :param filepath, str, location to be saved
        :param *args, str, names of datasets that you want in the export file
        :param, g_size_or_loss, int, minimum size of groups
            (positive, or #datasets - if 0 or negative)
        :param, c_size_or_loss, int, minimum size of clique
            (positive, or #datasets - if 0 or negative)
        :param, diameter, int, specifies maximum distance to be in a group
        :param, rank_and_remove, bool, True - removes nodes after adding to a group,
            False - produces all combinations that aren't subsets of larger groups.
        """

        if len(args) < 1:
            targ_datasets = list(self.datasets)
        else:
            targ_datasets = list(args)
        df = self.results(
            *targ_datasets,
            c_size_or_loss=c_size_or_loss,
            g_size_or_loss=g_size_or_loss,
            diameter=diameter,
            remove_rerank=remove_rerank,
        )
        if to_bytes:
            bytes_csv = df.to_csv().encode()
            return io.BytesIO(bytes_csv)
        if filepath is None:
            filepath = "results.csv"
        df.to_csv(filepath)
        return None

    def results(
        self,
        *args,
        c_size_or_loss=0,
        g_size_or_loss=0,
        diameter=1,
        remove_rerank=True,
    ):
        """
        Returns a results dataframe
        """
        if len(args) < 1:
            targ_datasets = list(self.datasets)
        else:
            targ_datasets = list(args)

        df = self.deconvolute(
            targ_datasets,
            remove_rerank,
            c_size_or_loss,
            g_size_or_loss,
            diameter,
        )
        for ds in df.columns:
            df = df.merge(
                self.datasets[ds],
                how="left",
                left_on=ds,
                right_index=True,
                suffixes=(None, f"_{ds}"),
            )
        return df

    @staticmethod
    def compress_graph(graph):
        """
        Compresses network into single highest rank non-directed graph
        """
        graph = graph.copy()
        selected_edges = [(u, v) for u, v, e in graph.edges(data=True) if e["rank"] > 0]
        for edge in selected_edges:
            graph.remove_edge(*edge)
        compressed = graph.to_undirected(reciprocal=True)

        # clear edge weights and update
        for u, v, d in compressed.edges(data=True):
            compressed[u][v]["score"] = 0
        for u, v, d in graph.edges(data=True):
            if u in compressed.neighbors(v):
                compressed[u][v]["score"] += d["score"]
        return compressed

    @staticmethod
    def get_rank_groups(graph, max_size, c_size=None, g_size=None, diameter=1):
        """
        Given a subgraph, ranks all possible groups based on provided info.
        Returns a list of ({group members}, clique_size, group_size,
        number of edges, avg edge weight)
        Ranking is accomplished by sorting by group size, clique size,
        number of edges, then edge weight (lowest to highest)
        """

        if c_size is None:
            c_size = max_size

        if g_size is None:
            g_size = max_size

        # some sanity checks
        if c_size > g_size:
            raise UserWarning(
                "You have indicated a clique size minimum greater than the group size m"
                "inimum. This is likely unintentional."
            )

        if max_size + 1 < c_size + diameter:
            raise UserWarning(
                "Your diameter and clique size setting is likely a mistake. If you use "
                "a diameter of 2, set 'c_size_or_loss' to -1, or any number less than t"
                "he number of datasets. For a diameter of 3, set 'c_size_or_loss' to -2"
                " or any number that is 2 below the number of datasets."
            )

        # use cliques as seeds
        cliques = set(
            frozenset(clique)
            for clique in nx.find_cliques(graph)
            if len(clique) >= c_size
        )

        potential_groups = set()

        if diameter == 1:
            potential_groups = cliques
        else:
            for clique in cliques:
                if diameter == 2:
                    # Try adding each nodes neighbors
                    neighbors = {
                        node: set(graph.neighbors(node)) - clique for node in clique
                    }
                    potential_groups.update(
                        set(
                            frozenset(clique.union(neighbors[node]))
                            for node in clique
                            if len(clique.union(neighbors[node])) >= g_size
                        )
                    )
                if diameter == 3:
                    # Initialize valid_nodes as a defaultdict of sets
                    valid_nodes = collections.defaultdict(set)

                    # Populate valid_nodes with neighbors grouped by partite
                    for clique_node in clique:
                        for neighbor_node in set(graph.neighbors(clique_node)) - clique:
                            partite = neighbor_node.split("__")[0]
                            valid_nodes[partite].add(neighbor_node)
                    target_g_size = sum(
                        1 for partite in valid_nodes.values() if partite
                    ) + len(clique)
                    if target_g_size < g_size:
                        continue
                    # add all frozensets to potential groups
                    non_empty = [nodes for nodes in valid_nodes.values() if nodes]
                    potential_groups.update(
                        set(
                            frozenset(combination).union(frozenset(clique))
                            for combination in product(*non_empty)
                        )
                    )

        # remove subsets
        filtered_groups = []
        for g1 in potential_groups:
            for g2 in potential_groups:
                if g1 == g2:
                    continue
                if g1.issubset(g2):
                    break
            else:
                filtered_groups.append(g1)

        # score and return
        final_groups = []
        for g in filtered_groups:
            sg = graph.subgraph(g)
            # score and record
            if (sg.number_of_edges()) == 0:
                weight = 0
            else:
                weight = (
                    sum(data["score"] for _, _, data in sg.edges(data=True))
                    / sg.number_of_edges()
                )
            max_clique_size = 0
            for clique in nx.find_cliques(sg):
                if len(clique) > max_clique_size:
                    max_clique_size = len(clique)
            final_groups.append(
                (set(g), len(sg), max_clique_size, len(sg.edges()), weight)
            )
        final_groups = sorted(final_groups, key=lambda x: (-x[1], -x[2], -x[3], x[4]))
        return final_groups

    def deconvolute(
        self,
        datasets,
        remove_rerank=True,
        c_size_or_loss=0,
        g_size_or_loss=0,
        diameter=1,
    ):
        """
        Makes a dataframe where each row is a feature, and each column is a dataset
        Entry is null if there is no matching clique
        """
        if g_size_or_loss <= 0:
            g_size_or_loss = len(self.datasets) + g_size_or_loss
        if c_size_or_loss <= 0:
            c_size_or_loss = len(self.datasets) + c_size_or_loss

        # compressed = self.compress_graph(self.graph)
        df_indices = {ds: [] for ds in self.datasets}
        results = []
        for dsg in nx.weakly_connected_components(self.graph):
            # skip graphs smaller than size
            if len(dsg) < g_size_or_loss:
                continue
            compressed = self.compress_graph(self.graph.subgraph(dsg))
            for sg in nx.connected_components(compressed):
                sg = compressed.subgraph(sg).copy()
                sg_results = []
                if remove_rerank:
                    while True:
                        new = self.get_rank_groups(
                            sg,
                            len(self.datasets),
                            c_size_or_loss,
                            g_size_or_loss,
                            diameter,
                        )
                        if len(new) == 0:
                            break
                        sg_results.append(new[0])
                        sg.remove_nodes_from(new[0][0])
                else:
                    sg_results = self.get_rank_groups(
                        sg,
                        len(self.datasets),
                        c_size_or_loss,
                        g_size_or_loss,
                        diameter,
                    )

                # remove those which don't have a dataset we're interested in
                ds_filtered = []
                for g in sg_results:
                    group_ds = [feature_name.split("__")[0] for feature_name in g[0]]
                    if len(set(datasets).intersection(group_ds)) >= 1:
                        ds_filtered.append(g)

                results.extend(ds_filtered)

        # parse results
        for row in results:
            # add an empty row fore everything, then fill in with features
            for df in df_indices.keys():
                df_indices[df].append(None)
            for feature in row[0]:
                ds_name, feature_name = feature.split("__")
                df_indices[ds_name][-1] = feature_name
        return pd.DataFrame.from_dict(df_indices)

    def align(self, ds1=None, ds2=None):
        """
        Performs all steps for generating the network
        """
        np.seterr(divide="ignore")
        self.gen_anchors(ds1, ds2)
        self.gen_coarse(ds1, ds2)
        self.gen_scalers(ds1, ds2)
        self.gen_scaled_values(ds1, ds2)
        self.gen_stds(ds1, ds2)
        self.gen_matches(ds1, ds2)
        self.gen_graph(ds1, ds2)
        np.seterr(divide="warn")

    def prescale(self):
        """
        Prescales all datasets
        """
        for ds in self.datasets:
            if ds not in self.prescalers:
                continue
            df = self.datasets[ds]
            for descriptor in self.descriptors:
                if descriptor not in self.prescalers[ds]:
                    continue
                prescalers = self.prescalers[ds][descriptor]
                try:
                    descriptor_col = df.columns.get_loc(descriptor)
                except KeyError as e:
                    raise KeyError(
                        f"The descriptor {descriptor} was not found in dataset {ds}."
                    ) from e
                if f"{descriptor}_Original" not in df.columns:
                    df.insert(
                        descriptor_col + 1, f"{descriptor}_Original", df[descriptor]
                    )
                df[descriptor] = scale_to(
                    df[f"{descriptor}_Original"],
                    prescalers[0],
                    prescalers[1],
                    self.descriptors[descriptor],
                )

    def gen_anchors(self, ds1=None, ds2=None):
        """
        Creates and sets anchors for all datasets specified.
        If either ds1 or ds2 is None, calculates anchors fall all datasets.
        """
        LOGGER.info("Generating anchors...")
        self.prescale()
        if ds1 is None or ds2 is None:
            datasets = self.datasets
        else:
            datasets = set()
            for ds in [ds1, ds2]:
                if isinstance(ds, str):
                    datasets.add(ds)
                else:
                    datasets.update(ds)
        coarse_params = copy.deepcopy(
            {k: self.default_coarse_params[k] for k in self.descriptors}
        )

        # copy the mode into the params
        for descr_param in coarse_params:
            coarse_params[descr_param].update({"mode": self.descriptors[descr_param]})

        for ds in datasets:
            if self.mz_col in set(self.descriptors.keys()):
                self.anchors[ds] = anchors_mz(
                    self.datasets[ds],
                    coarse_params,
                    self.mz_col,
                    self.remove_all,
                )
            else:
                self.anchors[ds] = anchors(
                    self.datasets[ds], coarse_params, self.remove_all
                )

    @dataset_loops("coarse_matches")
    def gen_coarse(self, ds1, ds2):
        """
        Generates all coarse anchors between datasets
        """
        coarse_params = self.default_coarse_params.copy()
        try:
            coarse_params.update(self.coarse_params[ds1][ds2])
        except KeyError:
            pass

        # remove coarse params which aren't in the descriptors
        coarse_params = copy.deepcopy({k: coarse_params[k] for k in self.descriptors})
        for descr_param in coarse_params:
            coarse_params[descr_param].update({"mode": self.descriptors[descr_param]})
        # extract anchors, match, and convert results to a dataframe
        if self.mz_col in self.descriptors.keys():
            res = simple_match_mz(
                self.datasets[ds1].loc[self.anchors[ds1], :],
                self.datasets[ds2].loc[self.anchors[ds2], :],
                coarse_params,
                self.mz_col,
            )
        else:
            res = simple_match(
                self.datasets[ds1].loc[self.anchors[ds1], :],
                self.datasets[ds2].loc[self.anchors[ds2], :],
                coarse_params,
            )
        self.coarse_matches[ds1][ds2] = pd.DataFrame(res)
        self.coarse_matches[ds1][ds2].columns = [ds1, ds2]

    @dataset_loops("scalers")
    def gen_scalers(self, ds1, ds2):
        """
        Creates and saves scalers for rt, mz, and intensity.
        These autogenerated scalers will not overwrite scalers already in place.
        """
        scaler_params = self.default_scaler_params.copy()
        try:
            scaler_params.update(self.scaler_params[ds1][ds2])
        except KeyError:
            pass
        coarse_match = self.coarse_matches[ds1][ds2]
        ds1_matches = self.datasets[ds1].loc[coarse_match[ds1], :]
        ds2_matches = self.datasets[ds2].loc[coarse_match[ds2], :]
        if ds2 not in self.scalers[ds1]:
            self.scalers[ds1][ds2] = {}

        scalers = {}
        for descr, mode in self.descriptors.items():
            if descr not in self.scalers[ds1][ds2]:
                scaled = calc_scalers(
                    ds1_matches[descr].values,
                    ds2_matches[descr].values,
                    smoothing=scaler_params["smoothing_method"],
                    mode=mode,
                    **scaler_params["smoothing_params"],
                )
                scalers[descr] = zip(*scaled)
        self._set_scalers(scalers, ds1, ds2, rec=False)

    @dataset_loops("scaled_values")
    def gen_scaled_values(self, ds1, ds2):
        """
        Scales the rt, mz, and intensity values for all datasets
        """
        # extract anchors, match, and convert results to a dataframe
        scaled = pd.DataFrame(index=self.datasets[ds1].index)
        scalers = self.scalers[ds1][ds2]
        for descr, mode in self.descriptors.items():
            scaled[descr] = scale_to(
                self.datasets[ds1][descr],
                scalers[descr]["x"],
                scalers[descr]["y"],
                mode=mode,
            )
        self.scaled_values[ds1][ds2] = scaled

    @dataset_loops("stds")
    def gen_stds(self, ds1, ds2):
        """
        Calculates standard deviations to be used for matching.
        """
        coarse_matches = self.coarse_matches[ds1][ds2]
        scalers = self.scalers[ds1][ds2]
        self.stds[ds1][ds2] = {}
        for descr, mode in self.descriptors.items():
            ds1_vals = self.datasets[ds1].loc[coarse_matches[ds1], descr]
            ds2_vals = self.datasets[ds2].loc[coarse_matches[ds2], descr]
            s_ds1_vals = scale_to(
                ds1_vals, scalers[descr]["x"], scalers[descr]["y"], mode=mode
            )
            if mode == "linear":
                res = (s_ds1_vals - ds2_vals).values
            if mode == "ppm":
                res = ((s_ds1_vals - ds2_vals) / s_ds1_vals * 1_000_000).values
            if mode == "log10":
                res = (np.log10(s_ds1_vals) - np.log10(ds2_vals)).values

            self.stds[ds1][ds2][descr] = np.std(res[~is_outlier(res)])

    @dataset_loops("matches")
    def gen_matches(self, ds1, ds2):
        """
        Generates best matches for each dataset
        """
        # generate the weights and cutoffs
        # if no "custom" cutoffs or weights are present, use default
        cutoffs = {}
        weights = {}
        for descr in self.descriptors:
            try:
                cutoffs[descr] = self.cutoffs[ds1][ds2][descr]
            except KeyError:
                try:
                    cutoffs[descr] = self.default_cutoffs[descr]
                except KeyError:
                    cutoffs[descr] = self.default_cutoff
            try:
                weights[descr] = self.weights[ds1][ds2][descr]

            except KeyError:
                try:
                    weights[descr] = self.default_weights[descr]
                except KeyError:
                    weights[descr] = self.default_weight
        if ds1 not in self.scores:
            self.scores[ds1] = {}

        if self.mz_col in self.descriptors.keys():
            self.matches[ds1][ds2], self.scores[ds1][ds2] = score_match_mz(
                self.scaled_values[ds1][ds2],
                self.datasets[ds2],
                self.stds[ds1][ds2],
                self.descriptors,
                cutoffs,
                weights,
                self.mz_col,
            )
        else:
            self.matches[ds1][ds2], self.scores[ds1][ds2] = score_match(
                self.scaled_values[ds1][ds2],
                self.datasets[ds2],
                self.stds[ds1][ds2],
                self.descriptors,
                cutoffs,
                weights,
            )

    @dataset_loops()
    def gen_graph(self, ds1, ds2):
        """
        Generates the network graph to be used for feature selection
        """
        nodes = [f"{ds1}__{str(i)}" for i in self.datasets[ds1].index]
        self.graph.add_nodes_from(nodes)
        matches = self.matches[ds1][ds2]
        scores = self.scores[ds1][ds2]
        edges = []
        for col in matches:
            this_match = matches[col][~pd.isnull(matches[col])]
            edge_scores = scores[col][~pd.isnull(matches[col])].values
            sources = this_match.index
            sources = [ds1 + "__" + str(name) for name in sources]
            targets = this_match.values
            targets = [ds2 + "__" + str(name) for name in targets]
            edges.extend(
                [
                    tuple([s, t, {"rank": col, "score": score}])
                    for s, t, score in zip(sources, targets, edge_scores)
                ]
            )
        self.graph.add_edges_from(edges)

    def report(
        self,
        datasets_1=None,
        datasets_2=None,
        filepath="report.pdf",
        to_bytes=False,
    ):
        """
        Creates a report of the alignment parameters.
        """
        if datasets_1 is None:
            datasets_1 = list(self.datasets)
        if isinstance(datasets_1, str):
            datasets_1 = [datasets_1]

        if datasets_2 is None:
            datasets_2 = list(self.datasets)
        if isinstance(datasets_2, str):
            datasets_2 = [datasets_2]
        if to_bytes:
            file_handle = io.BytesIO()
        else:
            file_handle = filepath

        # format plot text
        (
            title_wrap_width,
            x_label_wrap_width,
            y_label_wrap_width,
            plot_padding,
        ) = (
            55,
            65,
            40,
            1,
        )

        # create scaling report
        with PdfPages(file_handle) as pdf:
            for ds1 in datasets_1:
                for ds2 in datasets_2:
                    if (
                        ds1 == ds2
                        or ds1 not in self.scalers
                        or ds2 not in self.scalers[ds1]
                    ):
                        continue
                    for results in [False, True]:
                        if results and (
                            ds1 not in self.matches or ds2 not in self.matches[ds1]
                        ):
                            continue
                        plt.figure(figsize=(24, 14))
                        eval_chart(self, ds1, ds2, show=False, results=results)
                        for ax in plt.gcf().get_axes():
                            ax.set_title(
                                "\n".join(
                                    textwrap.wrap(ax.get_title(), title_wrap_width)
                                ),
                                fontsize=10,
                            )
                            ax.set_xlabel(
                                "\n".join(
                                    textwrap.wrap(ax.get_xlabel(), x_label_wrap_width)
                                ),
                                fontsize=8,
                                labelpad=5,
                            )
                            ax.set_ylabel(
                                "\n".join(
                                    textwrap.wrap(ax.get_ylabel(), y_label_wrap_width)
                                ),
                                fontsize=8,
                                labelpad=5,
                            )
                        plt.tight_layout(pad=plot_padding)
                        pdf.savefig()
                        plt.close("all")
        if to_bytes:
            file_handle.seek(0)
            return file_handle
        return None

    def explain(
        self,
        feature_tuples=None,
        annotation_ids=None,
        show_plot=True,
        label_edges=False,
        compress=True,
    ):
        """
        Creates a directed subgraph of feature's weakly interconnected nodes

        """
        combined_subgraph = nx.DiGraph()
        target_nodes = []
        if feature_tuples:
            if len(feature_tuples) == 2 and isinstance(feature_tuples[0], str):
                feature_tuples = [feature_tuples]
            for dataset, compound_id in feature_tuples:
                df = self.datasets[dataset]
                cmp_idx = df.loc[
                    df[self.feature_name] == compound_id, self.feature_name
                ].index[0]
                node_name = f"{dataset}__{cmp_idx}"

                # Find the target node
                target_node = None
                for node, data in self.graph.nodes(data=True):
                    if node == node_name:
                        target_node = node
                        target_nodes.append(target_node)
                        break
                if target_node is None:
                    continue

                component_subgraph = self.get_component(self.graph, target_node)
                if compress:
                    component_subgraph = self.compress_graph(component_subgraph)
                component_subgraph = component_subgraph.to_directed()
                combined_subgraph = nx.compose(combined_subgraph, component_subgraph)

        if annotation_ids:
            if isinstance(annotation_ids, str):
                annotation_ids = [annotation_ids]
            for dataset_name, df in self.datasets.items():
                matching_indices = df[
                    df[self.annotation_col].isin(annotation_ids)
                ].index
                for idx in matching_indices:
                    node_name = f"{dataset_name}__{idx}"
                    if node_name in self.graph:
                        target_nodes.append(node_name)
                        component_subgraph = self.get_component(self.graph, node_name)
                        if compress:
                            component_subgraph = self.compress_graph(component_subgraph)

                        # Convert to directed graph if not already
                        component_subgraph = component_subgraph.to_directed()
                        combined_subgraph = nx.compose(
                            combined_subgraph, component_subgraph
                        )
        node_labels = {}
        for node in combined_subgraph:
            [n_dataset, n_idx] = node.split("__")
            try:
                n_annotation = self.datasets[n_dataset].loc[n_idx, self.annotation_col]
            except:
                n_annotation = ""

            if pd.isnull(n_annotation):
                n_annotation = ""
            n_cmp_id = self.datasets[n_dataset].loc[n_idx, self.feature_name]

            node_labels[node] = f"{n_dataset}-{n_cmp_id}\n{n_annotation}"

        self.plot_subgraph(
            combined_subgraph,
            show_plot=show_plot,
            label_edges=label_edges,
            node_labels=node_labels,
            selected_nodes=target_nodes,
        )
        return combined_subgraph

    def get_component(self, graph, target_node):
        """
        Returns a subgraph of a target node's weakly interconnected nodes
        """
        visited_nodes = set()
        # Traverse outgoing edges
        for u, v in nx.bfs_edges(graph, source=target_node):
            visited_nodes.add(u)
            visited_nodes.add(v)
        # Traverse incoming edges
        for u, v in nx.bfs_edges(graph.reverse(), source=target_node):
            visited_nodes.add(u)
            visited_nodes.add(v)
        component_subgraph = graph.subgraph(visited_nodes).copy()
        return component_subgraph

    def plot_subgraph(
        self,
        subgraph,
        filepath="subgraph_plot.pdf",
        to_bytes=False,
        label_edges=True,
        show_plot=True,
        node_labels=None,
        selected_nodes=None,
    ):
        """
        Plots the subgraph.
        """
        if to_bytes:
            file_handle = io.BytesIO()
        else:
            file_handle = filepath

        plt.figure(figsize=(8, 6))
        pos = nx.spring_layout(subgraph)

        shape_cycle = [
            "o",  # Circle
            "^",  # Triangle
            "s",  # Square
            "p",  # Pentagon
            "*",  # Star
            "h",  # Hexagon
            "D",  # Diamond
            "d",  # Thin diamond
            "P",  # Plus
            "X",  # Cross
            "8",  # Figure 8
        ]
        color_cycle = [
            "red",
            "green",
            "blue",
            "orange",
            "purple",
        ]
        dataset_markers = {}

        for node in subgraph.nodes():
            dataset, _ = node.split("__")
            if dataset not in dataset_markers:
                dataset_markers[dataset] = [
                    shape_cycle[len(dataset_markers) % len(shape_cycle)],
                    color_cycle[len(dataset_markers) % len(color_cycle)],
                ]

        selected_node_list = []
        non_selected_node_list = []
        node_styles = {}

        for node in subgraph.nodes():
            dataset, _ = node.split("__")
            shape, color = dataset_markers[dataset]
            node_styles[node] = {"shape": shape, "color": color}

            if selected_nodes and node in selected_nodes:
                selected_node_list.append(node)
            else:
                non_selected_node_list.append(node)

        # Draw non-selected nodes
        for dataset, (shape, color) in dataset_markers.items():
            nodes_with_style = [
                node
                for node in non_selected_node_list
                if node.split("__")[0] == dataset
            ]
            if nodes_with_style:
                nx.draw_networkx_nodes(
                    subgraph,
                    pos,
                    nodelist=nodes_with_style,
                    node_color="white",
                    edgecolors=color,
                    node_shape=shape,
                    node_size=300,
                    linewidths=1.5,
                )

        # Draw selected nodes
        for dataset, (shape, color) in dataset_markers.items():
            nodes_with_style = [
                node for node in selected_node_list if node.split("__")[0] == dataset
            ]
            if nodes_with_style:
                nx.draw_networkx_nodes(
                    subgraph,
                    pos,
                    nodelist=nodes_with_style,
                    node_color=color,
                    edgecolors=color,
                    node_shape=shape,
                    node_size=300,
                    linewidths=1.5,
                )
        # Draw edges
        bidirectional_edges = []
        monodirectional_edges = []

        for u, v in subgraph.edges():
            if subgraph.has_edge(v, u):
                bidirectional_edges.append((u, v))
            else:
                monodirectional_edges.append((u, v))
        nx.draw_networkx_edges(
            subgraph,
            pos,
            edgelist=bidirectional_edges,
            edge_color="black",
            arrowstyle="->",
        )
        nx.draw_networkx_edges(
            subgraph,
            pos,
            edgelist=monodirectional_edges,
            edge_color="lightgray",
            style="dashed",
            arrowstyle="->",
        )
        # Draw labels
        label_pos = {node: (x + 0.02, y + 0.02) for node, (x, y) in pos.items()}
        nx.draw_networkx_labels(subgraph, label_pos, labels=node_labels, font_size=8)

        # Draw edge labels
        if label_edges:
            edge_labels = nx.get_edge_attributes(subgraph, "score")
            edge_labels = {w: f"{edge_labels[w]:.2f}" for w in edge_labels}
            nx.draw_networkx_edge_labels(
                subgraph, pos, edge_labels=edge_labels, font_size=6
            )
        # Create legend
        legend_elements = [
            plt.Line2D(
                [0],
                [0],
                marker=style[0],
                color=style[1],
                label=dataset,
                markerfacecolor="white",
                markeredgecolor=style[1],
                markeredgewidth=1.5,
                markersize=8,
            )
            for dataset, style in dataset_markers.items()
        ] + [
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="black",
                label="Selected Nodes",
                markerfacecolor="black",
                markeredgewidth=1.5,
                markersize=8,
            )
        ]
        plt.legend(
            handles=legend_elements,
            title="Datasets",
            loc="best",
            bbox_to_anchor=(1.0, 0.5),
            fontsize="small",
            title_fontsize="small",
        )
        plt.title("Directed Subgraph")
        plt.tight_layout()
        plt.savefig(file_handle, format="pdf", bbox_inches="tight")
        if show_plot:
            plt.show()
        plt.close()
        if to_bytes:
            file_handle.seek(0)
            return file_handle
        return None


def calc_scalers(x1, x2, smoothing=None, mode="linear", **kwargs):
    """
    Creates and returns scalers from X1 and X2 coordinates
    i.e. Dataset 1 RTs - [7.1, 8.43, 9.5]
         Dataset 2 RTs - [7.24, 8.74, 9.7]
    Returns scalers in the proper space, e.g. log10 space for intensity, ppm for mz
         Dataset 1 RTs -     [7.1, 8.43, 9.5]
         Dataset 2 Scalers - [0.14, 0.31, 0.2]

    :param x1: Array or List
    :param x2: Array or List
    :param smoothing: str, Smoothing mode (None or "lowess" currently)
    :param mode: str, "linear" or "ppm"
    :return: Numpy array
    """
    x1 = np.array(x1)
    x2 = np.array(x2)

    if mode == "linear":
        y = x2 - x1
    elif mode == "ppm":
        y = (x2 - x1) / x1 * 1000000
    elif mode == "log10":
        y = np.log10(x2) - np.log10(x1)
        x1 = np.log10(x1)
    else:
        raise NotImplementedError(f"We do not have a scaling method for {mode}.")
    if smoothing is None:
        return (x1, y)

    # +/- <0.5 ppt random noise to break ties
    np.random.seed(0)
    noise = (np.random.random(size=len(x1)) - 0.5) * x1 / 10**12

    # allow us to override frac via kwargs
    temp_kwargs = kwargs.copy()
    if "frac" not in temp_kwargs:
        temp_kwargs["frac"] = 0.1

    delta = 0.01 * (max(x1) - min(x1))
    scalers = lowess(y, x1 + noise, return_sorted=False, delta=delta, **temp_kwargs)

    # On failure, depending on statsmodels version, either returns NaNs, or
    # the residuals. Try adding more until it works or reaches 1. On reaching 1,
    # it should be obvious that smoothing failed or it is over smoothed.
    while np.isnan(scalers).all() or np.array_equal(scalers, y):
        temp_kwargs["frac"] += 0.01
        if temp_kwargs["frac"] >= 1.0:
            break
        scalers = lowess(
            y,
            x1 + noise,
            return_sorted=False,
            delta=delta,
            **temp_kwargs,
        )
    return (x1, scalers)


def score_match(df1, df2, stds, descriptors, cutoffs, weights, top_n=1):
    """
    Returns scored best matches for a dataframe to another
    :param df1: DataFrame, first dataset
    :param df2: DataFrame, second dataset
    :param stds: Dict, standard deviations for "RT", "MT", and "Intensity"
    :param cutoffs: Dict, Multipliers for finding cutoffs
    :param weights: Dict, Multipliers for weights
    :param top_n: Int, number of potential matches to record
    :return: Dataframe, df1 index as index, and matching df2 index in top_n columns
    """
    results_indices = []
    results_scores = []

    np1 = df1[list(descriptors)].values
    np2 = df2[list(descriptors)].values

    for row in np1:
        potential_hits = np.full(len(np2), True, dtype=bool)
        ds1_values = {}

        # generate values and bounds
        for i, (descr, mode) in enumerate(descriptors.items()):
            ds1_values[descr] = row[i]
            if mode == "linear":
                bounds = (
                    ds1_values[descr] - cutoffs[descr] * stds[descr],
                    ds1_values[descr] + cutoffs[descr] * stds[descr],
                )
            elif mode == "ppm":
                bounds = (
                    ds1_values[descr]
                    - (ds1_values[descr] / 1000000 * cutoffs[descr]) * stds[descr],
                    ds1_values[descr]
                    + (ds1_values[descr] / 1000000 * cutoffs[descr]) * stds[descr],
                )
            elif mode == "log10":
                bounds = (
                    10 ** (np.log10(ds1_values[descr]) - cutoffs[descr] * stds[descr]),
                    10 ** (np.log10(ds1_values[descr]) + cutoffs[descr] * stds[descr]),
                )
            else:
                raise NotImplementedError("Not sure how to handle the mode ", mode)

            # whittle hits down to list
            potential_hits = potential_hits & (
                (np2[:, i] >= bounds[0]) & (np2[:, i] <= bounds[1])
            )

        hits_index = df2.index[potential_hits]

        # score and sort results
        scores = []
        for hit in np2[potential_hits]:
            ds2_values = {descr: hit[i] for i, descr in enumerate(list(descriptors))}
            scores.append(
                calc_score(ds1_values, ds2_values, descriptors, stds, weights)
            )

        # append the top n
        this_hit = []
        this_score = []
        for score, hit in sorted(zip(scores, hits_index))[:top_n]:
            this_score.append(score)
            this_hit.append(hit)
        results_scores.append(this_score)
        results_indices.append(this_hit)

    return pd.DataFrame(results_indices, index=df1.index), pd.DataFrame(
        results_scores, index=df1.index
    )


def anchors_mz(ds, match_params, mz_key, remove_all=True, priority="Intensity"):
    """
    Shortcut for anchors when a descriminating column like MZ is present.
    Allows for binary search
    :param ds: DataFrame
    :param match_params: dict, parameters for matching
    :param mz_key: str, descriptor for binary search
    :param remove_all: bool, flag to indicate removing all confusable feature
    :param priority: str, the column indicating the which features are highest quality
        i.e. higher intensity is a strong feature
    :return: list, the index names of the anchors
    """
    thresholds = {k: {} for k in match_params if k != mz_key}
    if not remove_all:
        columns = list(set(match_params.keys()) | {priority})
        ds = ds[columns].copy()
    else:
        ds = ds[[*match_params.keys()]].copy()
    ds["index"] = ds.index
    ds = ds.sort_values(by=mz_key).reset_index(drop=True)

    for descriptor in thresholds:
        thresholds[descriptor]["values"] = ds[descriptor].to_numpy()
        tol = (
            abs(match_params[descriptor]["upper"] - match_params[descriptor]["lower"])
            / 2
        )

        if match_params[descriptor]["mode"] == "log10":
            thresholds[descriptor]["lower"] = thresholds[descriptor]["values"] / 10**tol
            thresholds[descriptor]["upper"] = thresholds[descriptor]["values"] * 10**tol
        if match_params[descriptor]["mode"] == "linear":
            thresholds[descriptor]["lower"] = thresholds[descriptor]["values"] - tol
            thresholds[descriptor]["upper"] = thresholds[descriptor]["values"] + tol
        if match_params[descriptor]["mode"] == "ppm":
            thresholds[descriptor]["lower"] = thresholds[descriptor]["values"] - (
                thresholds[descriptor]["values"] / 1e6 * tol
            )
            thresholds[descriptor]["upper"] = thresholds[descriptor]["values"] + (
                thresholds[descriptor]["values"] / 1e6 * tol
            )

    mz_array = ds[mz_key].to_numpy()
    orig_index = ds["index"].to_numpy()
    ppm_tol = abs(match_params[mz_key]["upper"] - match_params[mz_key]["lower"]) / 2
    mz_upper = mz_array + mz_array * ppm_tol / 1e6
    mz_lower = mz_array - mz_array * ppm_tol / 1e6
    graph = nx.Graph()
    graph.add_nodes_from(orig_index)

    for i in range(len(mz_array)):
        j_start = np.searchsorted(mz_array, mz_lower[i], side="left")
        j_end = np.searchsorted(mz_array, mz_upper[i], side="right")
        if j_start >= j_end:
            continue
        hits = np.array([True] * (j_end - j_start))
        for descriptor in thresholds:
            candidates = thresholds[descriptor]["values"][j_start:j_end]
            hits = (
                hits
                & (candidates > thresholds[descriptor]["lower"][i])
                & (candidates < thresholds[descriptor]["upper"][i])
            )
        for j in np.flatnonzero(hits) + j_start:
            if j == i:
                continue
            graph.add_edge(orig_index[i], orig_index[j])
    to_keep = []
    for component in nx.connected_components(graph):
        if len(component) == 1:
            to_keep.append(next(iter(component)))
        elif not remove_all:
            component_list = list(component)
            max_idx = ds.loc[component_list, "Intensity"].idxmax()
            to_keep.append(max_idx)
    to_keep.sort()
    return to_keep


def score_match_mz(df1, df2, stds, descriptors, cutoffs, weights, mz_key, top_n=1):
    """
    Shortcut for score_match when a discriminating column, like "MZ" is present.
    Allows for a fast binary search
    Returns scored best matches for a dataframe to another
    :param df1: DataFrame, first dataset
    :param df2: DataFrame, second dataset
    :param stds: Dict, standard deviations for "RT", "MT", and "Intensity"
    :param cutoffs: Dict, Multipliers for finding cutoffs
    :param weights: Dict, Multipliers for weights
    :param mz_key: Str, descriptor label for binary search
    :param top_n: Int, number of potential matches to record
    :return: Dataframe, df1 index as index, and matching df2 index in top_n columns
    """
    results_indices = []
    results_scores = []

    other_keys = [k for k in descriptors if k != mz_key]
    df1 = df1[[mz_key] + other_keys].copy()
    df2 = df2[[mz_key] + other_keys].copy()
    df2 = df2.sort_values(by=mz_key)

    mz1 = df1[mz_key].to_numpy()
    mz2 = df2[mz_key].to_numpy()
    index2 = df2.index.to_numpy()

    mz_tol = cutoffs[mz_key] * stds[mz_key]
    mz_lower = mz1 - mz1 * mz_tol / 1e6
    mz_upper = mz1 + mz1 * mz_tol / 1e6

    # Prepare descriptor values and bounds
    vals1 = {}
    vals2 = {}
    for k in other_keys:
        v1 = df1[k].to_numpy()
        v2 = df2[k].to_numpy()
        tol = cutoffs[k] * stds[k]
        mode = descriptors[k]
        if mode == "log10":
            lower = v1 / 10**tol
            upper = v1 * 10**tol
        elif mode == "linear":
            lower = v1 - tol
            upper = v1 + tol
        elif mode == "ppm":
            lower = v1 - (v1 * tol / 1e6)
            upper = v1 + (v1 * tol / 1e6)
        else:
            raise ValueError(f"Unknown mode for {k}")
        vals1[k] = {"v": v1, "lo": lower, "hi": upper}
        vals2[k] = v2

    for i in range(len(mz1)):
        mz_i = mz1[i]
        j_start = np.searchsorted(mz2, mz_lower[i], side="left")
        j_end = np.searchsorted(mz2, mz_upper[i], side="right")

        if j_start >= j_end:
            results_indices.append([None] * top_n)
            results_scores.append([None] * top_n)
            continue

        # Apply descriptor filters
        hits = np.ones(j_end - j_start, dtype=bool)
        for k in other_keys:
            candidates = vals2[k][j_start:j_end]
            hits &= (candidates >= vals1[k]["lo"][i]) & (
                candidates <= vals1[k]["hi"][i]
            )

        if not np.any(hits):
            results_indices.append([None] * top_n)
            results_scores.append([None] * top_n)
            continue

        candidates_mz = mz2[j_start:j_end][hits]
        idx = index2[j_start:j_end][hits]
        score = ((((mz_i - candidates_mz) / mz_i) * 1e6) / stds[mz_key]) ** 2 * weights[
            mz_key
        ]
        for k in other_keys:
            v1 = vals1[k]["v"][i]
            v2 = vals2[k][j_start:j_end][hits]
            if descriptors[k] == "linear":
                score += ((v1 - v2) / stds[k]) ** 2 * weights[k]
            elif descriptors[k] == "log10":
                score += (np.log10(v1 / v2) / stds[k]) ** 2 * weights[k]
            elif descriptors[k] == "ppm":
                score += ((v1 - v2) / v1 * 1e6 / stds[k]) ** 2 * weights[k]

        top_idx = np.argsort(score)[:top_n]
        score_top = score[top_idx]
        index_top = idx[top_idx]

        pad_len = top_n - len(score_top)
        results_indices.append(index_top.tolist() + [None] * pad_len)
        results_scores.append(score_top.tolist() + [None] * pad_len)

    return (
        pd.DataFrame(results_indices, index=df1.index),
        pd.DataFrame(results_scores, index=df1.index),
    )


def simple_match_mz(df1, df2, match_params, mz_key):
    """
    Binary search match based on a discriminating mz_key with
    additional descriptor filtering.

    :param df1: DataFrame, source
    :param df2: DataFrame, target
    :param match_params: dict, format:
        {
            "MZ": {"lower": -5, "upper": 5, "mode": "ppm"},
            "RT": {"lower": -0.2, "upper": 0.2, "mode": "linear"},
            ...
        }
    :param mz_key: str, primary key used for binary search
    :return: list of [df1_index, ds2_index]
    """
    other_keys = [k for k in match_params if k != mz_key]
    df1 = df1[[mz_key] + other_keys].copy()
    df2 = df2[[mz_key] + other_keys].copy()
    df2 = df2.sort_values(by=mz_key)

    mz1 = df1[mz_key].to_numpy()
    mz2 = df2[mz_key].to_numpy()

    mz_lower = mz1 + mz1 * match_params[mz_key]["lower"] / 1e6
    mz_upper = mz1 + mz1 * match_params[mz_key]["upper"] / 1e6

    vals1 = {}
    vals2 = {}
    for k in other_keys:
        v1 = df1[k].to_numpy()
        v2 = df2[k].to_numpy()
        mode = match_params[k]["mode"]
        if mode == "log10":
            lower = v1 * 10 ** match_params[k]["lower"]
            upper = v1 * 10 ** match_params[k]["upper"]
        elif mode == "linear":
            lower = v1 + match_params[k]["lower"]
            upper = v1 + match_params[k]["upper"]
        elif mode == "ppm":
            lower = v1 + (v1 * match_params[k]["lower"] / 1e6)
            upper = v1 + (v1 * match_params[k]["upper"] / 1e6)
        else:
            raise ValueError(f"Unknown mode for {k}")
        vals1[k] = {"v": v1, "lo": lower, "hi": upper}
        vals2[k] = v2

    results = []
    for i in range(len(mz1)):
        j_start = np.searchsorted(mz2, mz_lower[i], side="left")
        j_end = np.searchsorted(mz2, mz_upper[i], side="right")
        if j_start >= j_end:
            continue

        hits = np.ones(j_end - j_start, dtype=bool)
        for k in other_keys:
            candidates = vals2[k][j_start:j_end]
            hits &= (candidates > vals1[k]["lo"][i]) & (candidates < vals1[k]["hi"][i])

        for m in df2.index[j_start:j_end][hits]:
            results.append([df1.index[i], m])
    return results


def anchors(ds, match_params, remove_all=True, priority="Intensity"):
    """
    Creates a list of anchors for a given dataset.
    :param ds: DataFrame
    :param match_params: dict, parameters for matching
    :param remove_all: bool, flag to indicate removing all confusable feature
    :param priority: str, the column indicating the which features are highest quality
        i.e. higher intensity is a strong feature
    :return: list, the index names of the anchors
    """
    if not remove_all:
        columns = list(set(match_params.keys()) | {priority})
        ds = ds[columns]
        ds = ds.sort_values(priority, ascending=False)
    else:
        ds = ds[[*match_params.keys()]]

    # convert to numpy for performance
    index = ds.index
    ds = ds.values

    inclusion = set()
    exclusion = set()
    for i, row in enumerate(ds):
        # if a feature has been marked for removal, skip it.
        # important for non-remove all algorithm
        if not remove_all and index[i] in exclusion:
            continue
        bool_array = np.full(len(index), True, dtype=bool)
        for j, param in enumerate(match_params.values()):
            descr_val = row[j]
            tol = abs(param["upper"] - param["lower"]) / 2
            if param["mode"] == "linear":
                lower = descr_val - tol
                upper = descr_val + tol
            elif param["mode"] == "ppm":
                lower = descr_val - (descr_val / 1_000_000 * tol)
                upper = descr_val + (descr_val / 1_000_000 * tol)
            elif param["mode"] == "log10":
                descr_val = float(descr_val)
                lower = descr_val / 10**tol
                upper = descr_val * 10**tol
            bool_array = bool_array & (ds[:, j] > lower) & (ds[:, j] < upper)

        results = ds[bool_array]
        result_indices = index[bool_array]
        if remove_all:
            if len(results) == 1:
                inclusion.add(index[i])
            elif len(results) > 1:
                # if there is more than one result, exclude all of them
                exclusion.update(set(result_indices))
        else:
            if len(results) == 1:
                # if the only result is itself add it and remove the others
                inclusion.add(index[i])
            elif len(results) > 1:
                if result_indices[0] == index[i]:
                    inclusion.add(index[i])
                    not_max = np.delete(result_indices, 0, 0)
                    exclusion.update(set(not_max))

    ds_anchors = list(inclusion - exclusion)
    ds_anchors.sort()
    return ds_anchors


def simple_match(df1, df2, match_params):
    """
    Simply iterates through rows and finds features that fall within a
     threshold in a reference dataset.
    In the future this function will be much more generalized
    :param df1: DataFrame, first dataset
    :param df2: DataFrame, second dataset
    :param match_params: dict, match parameters
    :return: list, 2D Nested containing indices for df1 and df2 matches
    """
    results = []

    np1 = df1[[*match_params.keys()]].values
    np2 = df2[[*match_params.keys()]].values
    for i, row in enumerate(np1):
        bool_array = np.full(len(df2), True, dtype=bool)
        for j, param in enumerate(match_params.values()):
            descr_val = row[j]
            if param["mode"] == "linear":
                lower = descr_val + param["lower"]
                upper = descr_val + param["upper"]
            elif param["mode"] == "ppm":
                lower = descr_val + (descr_val / 1_000_000 * param["lower"])
                upper = descr_val + (descr_val / 1_000_000 * param["upper"])
            elif param["mode"] == "log10":
                lower = descr_val * 10 ** param["lower"]
                upper = descr_val * 10 ** param["upper"]
            bool_array = bool_array & (np2[:, j] > lower) & (np2[:, j] < upper)

        for result in df2.index[bool_array]:
            results.append([df1.index[i], result])

    return results


def gen_intensity(dataset, intensity_col="Intensity", ignore=None):
    """
    calculate an intensity column, using all columns that are not named '
     RT','MZ','Compound_ID'
    overwrites the existing dataframe with another one where Intensity is at index 1

    :param dataset: pandas Dataframe where there are intensities to calculate
    """
    if ignore:
        calc_df = dataset.drop(columns=ignore)
    calc_df.fillna(0, inplace=True)
    dataset.insert(1, intensity_col, calc_df.mean(axis=1, numeric_only=True))


def scale_to(x, x_scalers, y_scalers, mode="linear"):
    """
    Adjusts the x values to the scalers provided.

    :param x: Numpy Array, unscaled values
    :param x_scalers: Numpy Array, x values for the scalers
    :param y_scalers: Numpy Array, y values for the scalers
    :return: Numpy Array, scaled values for x
    """
    if mode == "linear":
        pass
    elif mode == "ppm":
        y_scalers = y_scalers * x_scalers / 1_000_000
    elif mode == "log10":
        x = np.log10(x)
    else:
        raise NotImplementedError(f"We do not have a scaling method for {mode}.")

    x_min_y = y_scalers.iloc[np.argmin(x_scalers)]
    x_max_y = y_scalers.iloc[np.argmax(x_scalers)]
    scaled = x + interpolate.interp1d(
        x_scalers, y_scalers, bounds_error=False, fill_value=(x_min_y, x_max_y)
    )(x)
    if isinstance(scaled, pd.Series):
        scaled = scaled.values
    if mode == "log10":
        scaled = 10**scaled
    return scaled


def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False
    otherwise.
    Borrowed from https://stackoverflow.com/questions/22354094/pythonic-way-of-detecting
    -outliers-in-one-dimensional-observation-data

    Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.

    :param points: An numobservations by numdimensions array of observations
    :param thresh: The modified z-score to use as a threshold. Observations with
    :return: Numpy Array, A numobservations-length boolean array.
    """
    try:
        points = points.values
    except AttributeError:
        pass
    median = np.median(points, axis=0)
    diff = (points - median) ** 2
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)
    modified_z_score = 0.6745 * diff / med_abs_deviation
    return modified_z_score > thresh


def calc_score(v1, v2, descriptors, stds, weights):
    """
    Returns the score of two features, provided their parameters, standard deviations,
     and scoring weights
    :return: float, score
    """
    penalty = 0
    for descr, mode in descriptors.items():
        if mode == "linear":
            penalty += ((v1[descr] - v2[descr]) / stds[descr]) ** 2 * weights[descr]
        elif mode == "ppm":
            penalty += (
                (((v1[descr] - v2[descr]) / v1[descr]) * 1_000_000) / stds[descr]
            ) ** 2 * weights[descr]
        elif mode == "log10":
            penalty += (np.log10(v1[descr] / v2[descr]) / stds[descr]) ** 2 * weights[
                descr
            ]
        else:
            raise NotImplementedError("Not sure how to score mode ", mode)
    return penalty


def eval_chart(self, key1, key2, show=True, results=False):
    """
    Description
    """

    ds1 = self.datasets[key1]
    ds2 = self.datasets[key2]

    if not results:
        coarse_df = self.coarse_matches[key1][key2]
        ds1 = ds1.loc[coarse_df[key1]]
        ds2 = ds2.loc[coarse_df[key2]]
        title = "Coarse"
    else:
        # filter our datasets to just the matched features
        has_match = pd.notnull(self.matches[key1][key2][0])
        matches = self.matches[key1][key2].loc[has_match]
        ds1_matches = matches.index
        ds2_matches = matches[0]
        ds1 = ds1.loc[ds1_matches]
        ds2 = ds2.loc[ds2_matches]
        title = "Matched"

    plt.figure(figsize=(20, 12))
    for i, (descr, mode) in enumerate(self.descriptors.items()):
        descr_1 = ds1[descr].values
        descr_2 = ds2[descr].values
        scalers = self.scalers[key1][key2][descr].sort_values(by="x")

        # calculate coarse scalers, or pull matched scaled values
        if not results:
            s_descr_1 = scale_to(descr_1, scalers["x"], scalers["y"], mode=mode)
        else:
            s_descr_1 = self.scaled_values[key1][key2].loc[ds1_matches][descr]
        if mode == "linear":
            y_vals = descr_2 - descr_1
            s_y_vals = descr_2 - s_descr_1
        if mode == "ppm":
            y_vals = 1e6 * (descr_2 - descr_1) / descr_1
            s_y_vals = 1e6 * (descr_2 - s_descr_1) / s_descr_1
        if mode == "log10":
            y_vals = np.log10(descr_2) - np.log10(descr_1)
            s_y_vals = np.log10(descr_2) - np.log10(s_descr_1)
            descr_1 = np.log10(descr_1)

        max_length = 50
        if len(key1) > max_length:
            label1 = f"...{key1[-(max_length - 3):]}"
        else:
            label1 = key1
        if len(key2) > max_length:
            label2 = f"...{key2[-(max_length - 3):]}"
        else:
            label2 = key2

        # plot the deltas
        plt.subplot(len(self.descriptors), 3, 3 * i + 1)
        plt.scatter(
            descr_1,
            y_vals,
            alpha=0.4,
            s=3,
            rasterized=True,
        )
        plt.title(f"{descr} {label1} vs {label2}")
        plt.xlabel(f"{title} {descr} {label1}")
        plt.ylabel(f"{title} {descr} ({label2} - {label1}) : {mode}")
        plt.plot(scalers["x"], scalers["y"], color="red", rasterized=True)

        # plot the scaled
        plt.subplot(len(self.descriptors), 3, 3 * i + 2)
        plt.scatter(
            descr_1,
            s_y_vals,
            alpha=0.4,
            s=3,
            rasterized=True,
        )
        plt.title(f"{descr} {label1} vs {label2}")
        plt.xlabel(f"{title} {descr} {label1}")
        plt.ylabel(f"{title} {descr} ({label2} - {label1}) : {mode}")
        plt.axhline(y=0, color="red", rasterized=True)

        # plot histogram
        plt.subplot(len(self.descriptors), 3, 3 * i + 3)
        plt.hist(
            x=(s_y_vals),
            bins="auto",
            alpha=0.7,
            rwidth=0.85,
        )
        plt.title(f"{descr} Histograms")
        plt.xlabel(f"{title} {descr} ({label2} - {label1})")
        plt.tight_layout()

    if show:
        plt.show()
