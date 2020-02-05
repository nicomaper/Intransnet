# Intransnet
Command line-based Python 3.7 script for the automatization of the integration of protein-protein interaction (PPI) data with RNA-Seq data.

I did this Python script during the developing of my Master Final Project. The objective was to automatize the integration of protein protein interaction data from a network with expression data derived from RNA-seq experiments. The objective is to simplify a given network according to a transcriptomic profilling experiment, by eliminating from the network the proteins which genes are not being expressed under certain experimental conditionds, tissues, cell types, etc.

## Getting Started
This script has three mandatory input arguments:
* A PPI network in the form of an edge list with only two columns (txt file)
* A data frame with the expression values of the network proteins genes from an RNA-seq experiment as FPKM, RPKM, TPM, or similar (txt file)
* A cut-off expression value (float)

The genes which expression is under the provided cut-off value will be deleted from the original given network, giving birth to a new version that contains only the nodes that are expressed under certain experimental conditions, tissues, cell types, etc. Instead of providing the cut-off value as a float, there is also the possibility of getting its automatic computation by providing the whole RNA-Seq dataset as a data frame (has the same structure as the previous one, but with all the genes and not only the network genes).

The ouput of the script is a data frame with the value of 18 network parameters for all the different networks generated.

Another optional arguments are:

* Option to generate the edge list of all the generated networks.
* Option to perform statistical comparison of all the networks among them, with the original and with a random network and and visualize the results.

### Prerequisites

To use the script it is required to have Python 3.7 installed. For further information on how to do this visit the following [link](https://tecadmin.net/install-python-3-7-on-ubuntu-linuxmint/)

Also it is important that the following Python packages are installed:
* Pandas
* Numpy
* NetworkX
* Scipy
* Scikit-learn
* Matplotlib
* Seaborn

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
