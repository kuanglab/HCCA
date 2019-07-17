# Hierarchical Canonical Correlation Analysis (HCCA)


### Dataset

Dataset used in the paper can be download [here](http://compbio.cs.umn.edu/HCCA/Data.zip).

Please make sure to unzip the `Data` folder to the same directory as the script files.

### Package

HCCA can be run using the following command:

```
c = 4;
perc = 0.85;
[Us,~] = HCCA({Xexp;Xmut;Xmet;Climate}, c, perc);
X = Us{end};
```
Where Xexp, Xmut, Xmet and Climate are, respectively, gene expression, mutation profile, DNA methylation and Geoclimate datasets. `c` indicates the condition number and `perc` is the percentage of eigenvalues to select

Please refer to `script_example_prediction.m` and `script_example_correlation_analysis.m` for more examples of how to run the package and baselines.
