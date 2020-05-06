RELEASE = True

import codecs
from setuptools import setup, find_packages
import sys, os

classifiers = """\
Development Status :: 5 - Production/Stable
Environment :: Console
Intended Audience :: Developers
Intended Audience :: Science/Research
License :: OSI Approved :: ISC License (ISCL)
Operating System :: OS Independent
Programming Language :: Python
Topic :: Scientific/Engineering
Topic :: Software Development :: Libraries :: Python Modules
"""

version = '0.1'

HERE = os.path.abspath(os.path.dirname(__file__))

def read(*parts):
    """
    Build an absolute path from *parts* and and return the contents of the
    resulting file.  Assume UTF-8 encoding.
    """
    with codecs.open(os.path.join(HERE, *parts), "rb", "utf-8") as f:
        return f.read()

setup(
        name='hpx',
        version=version,
        description="Python wrapper module for HPx simulations",
        packages=["hpx"],
        long_description=read("README.md"),
        classifiers=filter(None, classifiers.split("\n")),
        keywords='hydrological models',
        author='Joel Rahman',
        author_email='joel@flowmatters.com.au',
        url='https://github.com/flowmatters/py-hpx',
        #download_url = "http://cheeseshop.python.org/packages/source/p/Puppy/Puppy-%s.tar.gz" % version,
        license='ISC',
        py_modules=['hpx'],
        include_package_data=True,
        zip_safe=True,
        test_suite = 'nose.collector',
        install_requires=[
            'numpy',
            'pandas'
        ],
        extras_require={
            'test': ['nose'],
        },
)
