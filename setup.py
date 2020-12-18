import os
import re
import sys
import sysconfig
import platform
import subprocess

from distutils.version import LooseVersion
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from setuptools.command.test import test as TestCommand


class PyTest(TestCommand):
    def run_tests(self):
        import shlex

        # import here, cause outside the eggs aren't loaded
        import pytest

        errno = pytest.main(["tests/"])
        sys.exit(errno)


setup(
    name="pyclouds",
    version="0.1.1",
    author="Leif Denby",
    author_email="l.c.denby@leeds.ac.uk",
    description="Cloud model and thermodynamics",
    long_description="",
    zip_safe=False,
    install_requires=open("requirements.txt").readlines(),
    packages=["pyclouds", "pyclouds.ccfm", "pyclouds.ccfm.ccfmpython"],
    # add custom build_ext command
    cmdclass=dict(test=PyTest),
    tests_require=["pytest"],
)
