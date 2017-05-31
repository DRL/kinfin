import pip
from setuptools import setup, find_packages
from setuptools.command.install import install
from pip.req import parse_requirements
from codecs import open
from os import path

__version__ = '0.9'

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# get the dependencies and installs
install_reqs = parse_requirements(path.join(here, 'requirements.txt'), session=False)
reqs = [str(ir.req) for ir in install_reqs]

class OverrideInstall(install):

    """
    Emulate sequential install of pip install -r requirements.txt
    To fix numpy bug in scipy, scikit in py2
    """

    def run(self):
        for req in reqs:
            pip.main(["install", req])



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
      'Programming Language :: Python :: 2.7',
    ],
    keywords='Comparative genomics',
    packages=find_packages(exclude=['docs', 'tests*']),
    include_package_data=True,
    author='Dominik R Laetsch',
    cmdclass={'install': OverrideInstall},
    author_email='dominik.laetsch@gmail.com'
)
