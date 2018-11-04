#!/usr/bin/env python
"""The setup script."""

from setuptools import setup, find_packages
from setuptools.extension import Extension

from Cython.Build import cythonize, build_ext
import numpy

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'numpy',
    'cython',
]

setup_requirements = [
    'numpy',
    'numba',
    'cython',
]

test_requirements = [
    'pytest',
    'pytest-runner',
    'pytest-benchmark',
]

extensions = [
    Extension(
        "pyssa.pyssa_cython",
        ["pyssa/pyssa_cython.pyx"],
        include_dirs=[numpy.get_include()],
    ),
]

setup(
    author="Dileep Kishore, Srikiran Chandrasekaran",
    author_email='k.dileep1994@gmail.com',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="Python package for stochastic simulations",
    install_requires=requirements,
    license="Apache Software License 2.0",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='pyssa',
    name='pyssa',
    packages=find_packages(include=['pyssa']),
    ext_modules=cythonize(extensions),
    cmdclass={'build_ext': build_ext},
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/dileep-kishore/pyssa',
    version='0.1.0',
    zip_safe=False,
)
