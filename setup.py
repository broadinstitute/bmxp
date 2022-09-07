"""Setup for Eclipse"""
from setuptools import setup, find_namespace_packages

setup(
    name="bmxp",
    version="0.0.1",
    description="A package to perform calculations and match between "
    "nontargeted LCMS Datasets",
    packages=find_namespace_packages(include=["bmxp.*"]),
)
