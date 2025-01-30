## get_maps

Takes a list of Kegg KOs as input in a file and return a list of Kegg pathways with the number of KO including in each pathway.

```
get_maps.py -i FILE [-o FOLDER] [-l FOLDER]
```
  -i: file with KOs  
  -o: output folder [.\]  
  -l: log folder [.\]  

## build_pathways_matrix

The script takes as input pathway counts file generated with get_maps.py and create a count matrix.

```
build_pathways_matrix.py -f FOLDER [-o FILE]
```
  -f FOLDER   folder containing pathways count files  
  -o FILE     output matrix name [KO_matrix.csv]  

## rf_functions

Utility functions for random forests using the randomForest package in R.

```
rf_overfitting_test <- function(table, target, mtry, ntree, index_split, iteration, type, maxnode = NULL, print=TRUE)
```
Compute the difference in Mean Squared Error between the model and predictions. Usefull for determining if a model is overfitting or not.

```
best_rf <- function(table, target, mintree, maxtree, step)
```
Find the best combination of mtry and ntree values in order to obtain the best model.

```
remove_least_imp <- function(table, model, proportion)
```
Only keep most important predictors. Usefull for decomplexing model, lowering overfitting and calculation time.

## gbk_commands

Different functions to use with similar Genbank files. Useful for sequence vizualisation with tools like Clinker.

```
gbk_commands.py find -i [FOLDER]
```
Find genes common to all Genbank files. Duplicated genes are ignored.

  -i: Folder with all Genbank files

```
gbk_commands.py order [-h] -i FOLDER -o FOLDER -t STR
```
Re-order Genbank files considering a gene in common as origin.

  -i: Folder with all Genbank files
  -o: Output folder for all reordered Genbank files
  -t: Target gene considered as origin

```
gbk_commands.py extract [-h] -i FOLDER -o FOLDER -t STR -r NUM
```
Extract region surrounding a specified gene.

  -i: Folder with all Genbank files
  -o: Output folder for all extracted regions
  -t: Target gene
  -r: Sequence length to extract before and after the target

