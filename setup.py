"""Setup for Eclipse"""
import os
from setuptools import setup, find_namespace_packages
from distutils.core import Extension

HERE = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(HERE, "README.md"), encoding="utf-8") as f:
    LONG_DESCRIPTION = f.read()

with open(os.path.join(HERE, "requirements.txt"), encoding="utf-8") as f:
    REQUIREMENTS = f.read()

module1 = Extension
setup(
    name="bmxp",
    install_requires=REQUIREMENTS.split("\n"),
    version="0.0.17",
    description="LCMS Processing tools used by the Metabolomics Platform at the Broad"
    " Institute.",
    packages=find_namespace_packages(include=["bmxp.*"]),
    license="MIT",
    author="Daniel S. Hitchcock",
    author_email="daniel.s.hitchcock@gmail.com",
    long_description_content_type="text/markdown",
    long_description=LONG_DESCRIPTION,
    keywords=["LCMS", "Alignment", "Processing", "Metabolomics", "Clustering"],
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
    package_data={"": ["*.dll", "*.so"]},
    url="https://github.com/broadinstitute/bmxp",
)
