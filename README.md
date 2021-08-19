# TDA suspense
Analysis code for the master's thesis: Revealing brain network dynamics during the emotional state of suspense using topological data analysis

## Dataset

- Please note that the large size of the original dataset, the sharing-   restrictions, and the size limitations for repositories prevent us from uploading the full data here. You can obtain the original data from the CamCan Consortium or request the extracted data from Prof <a href="https://github.com/rschmaelzle"> Schmaelzle</a> or Prof. <a href="https://github.com/claregrall"> Grall </a>.
- Extracted ... are in the folder Data </li>


## Code
The notebooks to reproduce the analyses are in the Code folder. 
There are some considerations:
- The file `requirements.txt` specifies the Python and R packages used in the project.
- Most functions used in the project are in the `files _functions-py.py` and `functions-R.R` . This plain scripts were automatically generated and linked by jupytext from their respectives Jupyter notebooks.
- The `Mapper` package in R can be installed via [devtools](https://github.com/r-lib/devtools) package: 

```R
require("devtools")
devtools::install_github("aaolaveh/Mapper")
```
This personal Mapper's version is a forked project from the original [Mapper R package ](https://github.com/peekxc/Mapper) of Matt Piekenbrock. 

In addtion, Mapper construction relies on [simplextree](https://github.com/peekxc/simplextree). The current version is not yet compatible with Mapper, thus, install or maintain if asked the older version
```R
devtools::install_github("peekxc/simplextree@6e34926")
```








