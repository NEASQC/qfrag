# About KraChem
KraChem is Python library that focuses on fragmentation schemes for chemical molecules by using graph partitioning methods. Fragmentation of molecules and complexes is a very important topic in chemistry, particularly in computational chemistry.

Currently available as a beta version, KraChem implements two fragmentation schemes:

1. Intermolecular fragmentation: Separation of two or mode independent interacting moieties (typically two molecules).
2. Intramolecular fragmentation: Separation of two or more moieties which belong to the same molecule, chain or system.

While KraChem is intended to aid computational chemistry as a whole, it was developed in the context of quantum computational chemistry use-cases in the [NEASQC](https://www.neasqc.eu) project.

## KraChem user guide
To start using KraChem, download this repository. The [krachem](krachem) directory contains the library modules and functionalities, while the [input_molecules](input_molecules) directory contains sample input files which define the chemical molecules. The [results](results) directory is pre-populated with the fragmentation outputs for some of the examples using KraChem.

A user guide for KraChem is provided [here](krachem_user_guide.ipynb) as a Jupyter Notebook.

<!--## KraChem developer guide
An outline of the functionalities and modules implemented in the KraChem library is provided as the developer guide [here](krachem_developer_guide.md).-->

## Using and citing KraChem
KraChem is licensed under the terms of [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/).

To cite this work please use, `Bruno Chagas, Goar Sanchez-Sanz, Pablo Lauret Martínez de Rituerto, Venkatesh Kannan, "KraChem: Chemical Molecular Fragmentation", Irish Centre for High-End Computing (ICHEC), URL: https://github.com/NEASQC/krachem`.

## Contact
For any support with using or contributing to KraChem, write to support@ichec.ie.


<!-- KraChem is a beta-version Python library designed for working with chemical molecules, fragmentation, based on chemical files input and mathematical tools from graph theory. Therefore, this package has two main functionalities, and they are

1. Intermolecular Fragmentation: given a chemical compound, where we may have a collection of disconnected molecules we want to find
2. Intramolecular Fragmentation: given a molecule defined by its atoms and connections, we want to brake this molecule into smaller parts, so we have a collection of subsets of atoms.

There are two main packages where we develop our framework for fragmenting molecules, based on our approach of analysing chemical molecules as graphs and certain arrays. In this way, we use the following packages

1. numpy: a library for the python programming language, adding support for large, multi-dimensional arrays and matrices, along with a large collection of high-level mathematical functions to operate on these arrays

2. networkx : a python library for studying graphs and networks, where we encode our molecules as graphs and perform several routines and algorithms in our graphs 

Moreover, we implement some routines for plotting graphs and certain analysis, as well as packages for creating subsets and manipulating paths and directories, and they are

3. matplotlib: a plotting library for the python programming language and its numerical mathematics extension numpy, where we have some features for plotting figures

4. itertools: this module implements a number of iterator building blocks inspired by constructs from APL, Haskell, and SML. We use this package for combining parts of a chemical compound or atoms in a molecule into a collection of subsets

4. pathlib: this module offers classes representing filesystem paths with semantics appropriate for different operating systems, where we can access and manipulate files in folders


### Tutorial

We have a [quickstart tutorial](quickstart_guide.ipynb) if you want to know the basic features in this package.

For documentation [click here](https://www.overleaf.com/project/6123ad698e1e62f3766d1b08) -->
