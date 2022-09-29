<img src="https://cloud.githubusercontent.com/assets/167909/26763490/8f07758a-494b-11e7-8fb7-83b8153f4691.png" width="128"> 

Citing KinFin
-------------

> [Laetsch DR and Blaxter ML, 2017. KinFin: Software for Taxon-Aware Analysis of Clustered Protein Sequences. G3: Genes, Genomes, Genetics. Doi:10.1534/g3.117.300233](https://doi.org/10.1534/g3.117.300233)

Dependencies
------------
- UNIX system 

Installation
------------

- Create [Conda](https://conda.io/en/latest/miniconda.html) environment
    
```
$ conda create -n kinfin -c conda-forge docopt==0.6.2 scipy==0.19.0 matplotlib networkx==1.11 ete3
$ conda activate kinfin
```

- Clone github repo

```
$ git clone https://github.com/DRL/kinfin.git
```

- Run install script for fetching databases

```
$ cd kinfin
$ ./install
```

- Test kinfin 

```
$ ./test
```

Usage
-----

    $ ./kinfin -h

Documentation
-------------

[kinfin.readme.io](https://kinfin.readme.io)
