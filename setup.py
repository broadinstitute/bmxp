"""Setup for Eclipse"""
import os
from setuptools import setup, find_namespace_packages

HERE = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(HERE, "README.md"), encoding="utf-8") as f:
    LONG_DESCRIPTION = f.read()

with open(os.path.join(HERE, "requirements.txt"), encoding="utf-8") as f:
    REQUIREMENTS = f.read()


setup(
    name="bmxp",
    install_requires=REQUIREMENTS.split("\n"),
    version="0.0.3",
    description="A package to perform calculations and match between "
    "nontargeted LCMS Datasets",
    packages=find_namespace_packages(include=["bmxp.*"]),
    license="MIT",
    author="Daniel S. Hitchcock",
    author_email="daniel.s.hitchcock@gmail.com",
    long_description_content_type="text/markdown",
    long_description=LONG_DESCRIPTION,
    keywords=["LCMS", "Alignment", "Processing", "Metabolomics"],
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
    url="https://github.com/broadinstitute/bmxp",
)
