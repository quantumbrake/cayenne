#!/usr/bin/env python
"""The setup script."""

from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy as np

with open("README.md") as readme_file:
    readme = readme_file.read()

requirements = ["numpy", "Cython", "matplotlib"]

setup_requirements = ["Cython"]

test_requirements = ["pytest", "pytest-runner", "pytest-benchmark"]

ext_modules = [
    Extension(
        "cayenne.algorithms.direct",
        ["cayenne/algorithms/direct.pyx"],
        include_dirs=[np.get_include()],
    ),
    Extension(
        "cayenne.algorithms.tau_leaping",
        ["cayenne/algorithms/tau_leaping.pyx"],
        include_dirs=[np.get_include()],
    ),
    Extension(
        "cayenne.algorithms.tau_adaptive",
        ["cayenne/algorithms/tau_adaptive.pyx"],
        include_dirs=[np.get_include()],
    ),
    Extension("cayenne.utils", ["cayenne/utils.pyx"], include_dirs=[np.get_include()]),
]

setup(
    author="Dileep Kishore, Srikiran Chandrasekaran",
    author_email="k.dileep1994@gmail.com",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
    description="Python package for stochastic simulations",
    install_requires=requirements,
    license="Apache Software License 2.0",
    long_description=readme + "\n\n",
    long_description_content_type="text/markdown",
    include_package_data=True,
    keywords="cayenne stochastic gillepsie simulation",
    name="cayenne",
    packages=find_packages(exclude=["tests"]),
    setup_requires=setup_requirements,
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/Heuro-labs/cayenne",
    version="1.0.0",
    zip_safe=False,
    ext_modules=cythonize(
        ext_modules,
        annotate=True,
        compiler_directives={"binding": True, "linetrace": False, "profile": False},
    ),
)
