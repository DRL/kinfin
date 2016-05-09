### clusterbuster
Analysis of clustered proteome data of many species 

----------

#### Analyse clusters
```
src/kinfin.py data/nematoda/SpeciesIDs.txt \
              data/nematoda/SpeciesClassification.txt \
              data/nematoda/OrthologousGroups_I1.5.txt \
              data/nematoda/OrthologousGroups_I2.0.txt
```
#### Hierarchality between clusters at two inflation values (I1.5 and I2.0): i.e. "Jaccard-magic"
```
src/calculate_jaccard_distances.py cluster_stats.txt I1.5 I2.0
```
