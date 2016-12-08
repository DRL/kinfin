### kinFin
Taxon-aware analysis of clustered proteomes

----------

### Installation

#### Requirements
- Python 2.7+
- Python libraries:
    - Docopt (pip install docopt)
    - Matplotlib (pip install matplotlib)
    - Seaborn (pip install seaborn)
    - SciPy (pip install scipy)
    - ETE3 ToolKit (pip install ete3)

#### Analyse clusters
```
src/kinfin.py \
    -p test_data/test.SpeciesIDs.txt \
    -g test_data/test.OrthologousGroups.txt \
    -c test_data/test.SpeciesClassification.txt \
    -o test \
    -s test_data/test.SequenceIDs.txt \
    --functional_annotation test_data/test.functional_annotation.txt \
    -a test_data/ \
    -t test_data/test.tree.nwk
```

