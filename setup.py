import os
import re
import sys
import sysconfig
import platform
import subprocess

import versioneer
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


cmdclass = versioneer.get_cmdclass()
# add custom build_ext command
cmdclass["test"] = PyTest

setup(
    name="pyclouds",
    version=versioneer.get_version(),
    author="Leif Denby",
    author_email="l.c.denby@leeds.ac.uk",
    description="Cloud model and thermodynamics",
    long_description="",
    zip_safe=False,
    install_requires=open("requirements.txt").readlines(),
    packages=["pyclouds", "pyclouds.ccfm", "pyclouds.ccfm.ccfmpython"],
    cmdclass=cmdclass,
    tests_require=["pytest"],
)
