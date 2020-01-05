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
        "pyssa.algorithms.direct",
        ["pyssa/algorithms/direct.pyx"],
        define_macros=[("CYTHON_TRACE", "1")],
        include_dirs=[np.get_include()],
    ),
    Extension(
        "pyssa.algorithms.tau_leaping",
        ["pyssa/algorithms/tau_leaping.pyx"],
        define_macros=[("CYTHON_TRACE", "1")],
        include_dirs=[np.get_include()],
    ),
    Extension(
        "pyssa.algorithms.tau_adaptive",
        ["pyssa/algorithms/tau_adaptive.pyx"],
        define_macros=[("CYTHON_TRACE", "1")],
        include_dirs=[np.get_include()],
    ),
    Extension(
        "pyssa.utils",
        ["pyssa/utils.pyx"],
        define_macros=[("CYTHON_TRACE", "1")],
        include_dirs=[np.get_include()],
    ),
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
    keywords="pyssa stochastic gillepsie simulation",
    name="pyssa",
    packages=find_packages(exclude=["tests"]),
    setup_requires=setup_requirements,
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/Heuro-labs/pyssa",
    version="0.9.0",
    zip_safe=False,
    ext_modules=cythonize(
        ext_modules,
        annotate=True,
        compiler_directives={"binding": True, "linetrace": False, "profile": False},
    ),
)
