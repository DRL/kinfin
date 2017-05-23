from setuptools import setup, find_packages
from codecs import open
from os import path

__version__ = '0.9'

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# get the dependencies and installs
with open(path.join(here, 'requirements.txt'), encoding='utf-8') as f:
    all_reqs = f.read().split('\n')

try:
    import matplotlib
except ImportError:  # We do not have numpy installed
    build_requires = ['matplotlib>=2.0.2']
else:
    # If we're building a wheel, assume there already exist numpy wheels
    # for this platform, so it is safe to add numpy to build requirements.
    # See gh-5184.
    build_requires = (['matplotlib>=2.0.2'] if 'bdist_wheel' in sys.argv[1:] else [])

install_requires = [x.strip() for x in all_reqs if 'git+' not in x]
dependency_links = [x.strip().replace('git+', '') for x in all_reqs if x.startswith('git+')]


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
    build_requires=build_requires,
    install_requires=install_requires,
    dependency_links=dependency_links,
    author_email='dominik.laetsch@gmail.com'
)
