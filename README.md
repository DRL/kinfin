<img src="https://cloud.githubusercontent.com/assets/167909/26763490/8f07758a-494b-11e7-8fb7-83b8153f4691.png" width="128"> 

Citing KinFin
-------------

> [Laetsch DR and Blaxter ML, 2017. KinFin: Software for Taxon-Aware Analysis of Clustered Protein Sequences. G3: Genes, Genomes, Genetics. Doi:10.1534/g3.117.300233](https://doi.org/10.1534/g3.117.300233)

Dependencies
------------
- UNIX system (bash, wget, tar, gunzip) 
- Python 2.7
- ```pip```

Installation
------------

Create Conda environment

    $ conda create -n kinfin -c conda-forge docopt==0.6.2 scipy==0.19.0 matplotlib networkx==1.11 ete3
    $ conda activate blobtools
    $ conda install -c anaconda matplotlib docopt tqdm wget pyyaml git
    $ conda install -c bioconda pysam --update-deps

    $ git clone https://github.com/DRL/kinfin.git


Option A: Create Conda environment

conda create -n blobtools
conda activate blobtools
conda install -c anaconda matplotlib docopt tqdm wget pyyaml git
conda install -c bioconda pysam --update-deps

Tip: Check if samtools exists by executing the command 'samtools' in the commandline. If samtools complains about dependencies, simply run the pysam install twice.

Option B: Install dependencies via PIP

python setup.py install --user

    $ ./install

Usage
-----

    $ ./kinfin -h

Example
-------

    $ ./test

Documentation
-------------

[kinfin.readme.io](https://kinfin.readme.io)
