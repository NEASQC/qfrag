### About KraChem

KraChem is a Python library designed for working with chemical molecules, fragmentation, based on chemical files input and mathematical tools from graph theory. Therefore, this package has two main functionalities, and they are

1. Intermolecular Fragmentation: given a chemical compound, where we may have a collection of disconnected molecules we want to find
2. Intramolecular Fragmentation: given a molecule defined by its atoms and connections, we want to brake this molecule into smaller parts, so we have a collection of subsets of atoms.

There are two main packages where we develop our framework for fragmenting molecules, based on our approach of analysing chemical molecules as graphs and certain arrays. In this way, we use the following packages

1. numpy: a library for the python programming language, adding support for large, multi-dimensional arrays and matrices, along with a large collection of high-level mathematical functions to operate on these arrays

2. networkx : a python library for studying graphs and networks, where we encode our molecules as graphs and perform several routines and algorithms in our graphs 

Moreover, we implement some routines for plotting graphs and certain analysis, as well as packages for creating subsets and manipulating paths and directories, and they are
3. matplotlib: a plotting library for the python programming language and its numerical mathematics extension numpy, where we have some features for plotting figures
4. itertools: this module implements a number of iterator building blocks inspired by constructs from APL, Haskell, and SML. We use this package for combining parts of a chemical compound or atoms in a molecule into a collection of subsets
4. pathlib: this module offers classes representing filesystem paths with semantics appropriate for different operating systems, where we can access and manipulate files in folders

### Requirements

Numpy
Networkx
Openbabel
Matplotlib

### Tutorial

We have a [quickstart tutorial](https://git.ichec.ie/neasqc-frag/qfrag/-/blob/master/quickstart_guide.ipynb) if you want to know the basic features in this package.

For documentation [click here](https://www.overleaf.com/project/6123ad698e1e62f3766d1b08)
