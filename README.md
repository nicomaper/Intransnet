# Intransnet
Command line-based Python 3.7 script for the automatization of the integration of protein-protein interaction (PPI) data with RNA-Seq data.

I did this Python script during the developing of my Master Final Project. The objective was to automatize the integration of protein protein interaction data from a network with expression data derived from RNA-seq experiments. The objective is to simplify a given network according to a transcriptomic profilling experiment, by eliminating from the network the proteins which genes are not being expressed under certain experimental conditions, tissues, cell types, etc.

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

To use the script it is required to have [Python 3.7](https://tecadmin.net/install-python-3-7-on-ubuntu-linuxmint/) installed.

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
Python intransnet.py -h
```


## Built With

* [Python 3.7](https://www.python.org/) - The programming language used to write this program

## Authors

* **Nicolás Manosalva Pérez** - contact information: nicolas.manosalva@estudiante.uam.es

## License

No license
