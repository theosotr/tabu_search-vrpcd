from setuptools import setup, find_packages
from vrpcd import __version__, __author__

# Package info
PACKAGE_NAME = "tabu_vrpcd"
SHORT_DESCRIPTION = ('Tabu Search Algorithm for solving Vehicle Routing'
                     'Problem with Cross-Docking')

PACKAGES_ROOT = '.'
PACKAGES = find_packages(PACKAGES_ROOT)

# Package meta
CLASSIFIERS = []

# Package requirements
INSTALL_REQUIRES = ['networkx']

EXTRAS_REQUIRES = {}

TESTS_REQUIRES = []


setup(
    name=PACKAGE_NAME,
    version=__version__,
    author=__author__,
    author_email='theosotr@windowslive.com',
    licence='Apache v2',
    description=SHORT_DESCRIPTION,
    classifiers=CLASSIFIERS,
    packages=PACKAGES,
    package_dir={'': PACKAGES_ROOT},
    include_package_data=True,
    zip_safe=False,
    install_requires=INSTALL_REQUIRES,
    extras_require=EXTRAS_REQUIRES,
    tests_require=TESTS_REQUIRES,
)
