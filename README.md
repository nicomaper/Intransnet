# Intransnet
Command line-based Python 3.7 script for the automatization of the integration of protein-protein interaction (PPI) data with RNA-Seq data.

I did this Python script during the developing of my Master Final Project. The objective was to automatize the integration of protein protein interaction data from a network with expression data derived from RNA-seq experiments. The purpose is to simplify a given network according to a transcriptomic profilling experiment, by eliminating from the network the proteins which genes are not being expressed under certain experimental conditions, tissues, cell types, etc.

## Getting Started
This script has three mandatory input arguments:
* A PPI network in the form of an edge list with only two columns (txt file)
* A data frame with the expression values of the network proteins genes from an RNA-seq experiment as FPKM, RPKM, TPM, or similar (txt file)
* A cut-off expression value (float)

The genes which expression is under the provided cut-off value will be deleted from the original given network, giving birth to a new version that contains only the nodes that are expressed under certain experimental conditions, tissues, cell types, etc. Instead of providing the cut-off value as a float, there is also the possibility of getting its automatic computation by providing the whole RNA-Seq dataset as a data frame (it should have the same structure as the previous one, but with all the genes and not only the network genes).

The ouput of the script is a data frame with the value of 18 network parameters for all the different networks generated.

Other optional arguments are:

* Option to generate the edge list of all the generated networks.
* Option to perform statistical comparison of all the networks among them, with the original and with a random network and and visualize the results.

### Prerequisites

To use the script it is required to [install Python 3.7](https://tecadmin.net/install-python-3-7-on-ubuntu-linuxmint/).

For the script to work properly, the following Python packages should be installed:
* [Pandas](https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html)
* [Numpy](https://docs.scipy.org/doc/numpy/user/install.html)
* [NetworkX](https://networkx.github.io/documentation/networkx-1.1/install.html)
* [Scipy](https://scipy.org/install.html)
* [Scikit-learn](https://scikit-learn.org/stable/install.html)
* [Matplotlib](https://matplotlib.org/users/installing.html)
* [Seaborn](https://seaborn.pydata.org/installing.html)

## Tutorial

For the script to work, it is required to have on the same directory the python file intransnet.py, the network edgelist and the expression values data frame in a txt format. For automatic cut-off value computation, it should also be added the whole RNA-Seq experiment Count Matrix. Then proceed to navigate to the said directory.

To get help on what are the program options type on the command line:

```
python intransnet.py -h
```

____________________

To use the example network and data set included in the directory 'Example' providing a pre-obtained cut-off value (in this case a random integer), copy and paste on the command line:

```
python intransnet.py edge_list.txt net_count_matrix.txt -cf 8
```

A new directory will be created called 'intransnet_results', which contains the data frame with the results. **Important**: before running the program again please change the name of the directory, otherwise it will be rewritten and the results will be lost.

Please note that if a cut-off value is provided, its value cannot be higher than any of the expression values contained in the Count Matrix. Otherwise an AssertionError will rise.

______________________

To use the example network and data set included in the directory 'Example' but with the automatic computation of the cut-off value, copy and paste on the command line:

```
python intransnet.py edge_list.txt net_count_matrix.txt -cfc whole_count_matrix.txt
```

**Important**: Note that the -cf and -cfc optiones are mutually exclusive. This means that they cannot be used together, but for the program to function at least one of them must be provided.

____________

Optinal arguments:

* -el, --edgelists: option for the program to generate the edge lists of the generated networks
* -p, --plotting: option to generate plots with statistical analyses (pval_hm, par_hm, pca, violin)

Copy and paste the following in the command line to used the mentioned above options:

```
python intransnet.py edge_list.txt net_count_matrix.txt -cfc whole_count_matrix.txt -el -p pval_hm par_hm pca violin
```

The edge lists as well as the plots in a png format are stored in the 'intransnet_results' directory.

## Built With

* [Python 3.7](https://www.python.org/) - The programming language used to write this program

## Authors

* **Nicolás Manosalva Pérez** - contact information: nicolas.manosalva@estudiante.uam.es

## License

No license
