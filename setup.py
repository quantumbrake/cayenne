#!/usr/bin/env python
"""The setup script."""

from setuptools import setup, find_packages

with open("README.md") as readme_file:
    readme = readme_file.read()

requirements = ["numpy", "numba", "matplotlib"]

setup_requirements = []

test_requirements = ["pytest", "pytest-runner", "pytest-benchmark"]

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
    keywords="pyssa stochastic gillepsie simulation numba",
    name="pyssa",
    packages=find_packages(include=["pyssa"]),
    setup_requires=setup_requirements,
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/Heuro-labs/pyssa",
    version="0.7.1",
    zip_safe=False,
)
