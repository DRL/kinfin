import pip
from setuptools import setup, find_packages

__version__ = '1.1'

# Get the long description from the README file
with open('README.md', 'r') as readme:
    long_description = readme.read()

# get the dependencies and installs
with open('requirements.txt', 'r') as requirements:
    reqs = requirements.read().splitlines()

setup(
    name='kinfin',
    version=__version__,
    description='Taxon-aware analysis of clustered protein data',
    long_description=long_description,
    url='https://github.com/DRL/kinfin',
    download_url='https://github.com/DRL/kinfin/tarball/' + __version__,
    license='GnuGPL3',
    classifiers=[
      'Development Status :: 3 - Alpha',
      'Intended Audience :: Developers',
      'Programming Language :: Python :: 3',
    ],
    keywords='Comparative genomics',
    packages=find_packages(exclude=['docs', 'tests*']),
    include_package_data=True,
    author='Dominik R Laetsch',
    entry_points={
        'console_scripts': [
            "kinfin=src.kinfin:main",
            ],
        },
    author_email='dominik.laetsch@gmail.com'
)
