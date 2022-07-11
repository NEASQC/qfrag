"""
Python Package:
    KraChem (version 1.0, 2021)
Authors:
    B. Chagas, G. Sánchez, V. Kannan
Organization:
    Irish Centre for High-End Computing (ICHEC)
    National University of Ireland Galway (NUIG)
"""

from itertools import combinations
from pathlib import Path
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import os
import shutil
from openbabel import openbabel


"""
Dictionary for coloring atoms as vertices in a graph
TODO: Add more atoms and colors
"""
node_color = {'H' : '#9ee1f7' ,'C' : '#cfc0c3' , 'O' : '#ffa1b2', 
              'N' : '#ce85ed' ,'P' : '#d1bb64' , 'Cl': '#6ccc87',
              'S' : '#b87163' ,'Na': '#e3c76d' , 'F' : '#d4ae31',
              'Br': '#c79202' ,'Mg': '#4dbda8' , 'Ca': '#857dab'}

class Molecule:
    """
    The Molecule class is responsible for transforming a chemical file
    into a graph and fragmenting between molecules in a compound or
    atoms which are considered leaves in a graph.
    
    
    Attributes
    ----------
        molecule_number: int
            indicates couting molecule objects in this class
        name: str
            name of the molecule given by the user
        path: str
            path of a chemical file
        atoms: list(str)
            contains strings which indicate the name of atoms
        edges: list(tuples(int))
            connections/bonds between each pair of atoms
        bonds: list(int)
            type of bonding (single, double or triple)
        graph: Graph from networkx
            a graph in which we have a converted a molecule/compound into this class
        inter_components: list(list(int))
            contains intermolecular fragments from a chemical compound
            
    
    Methods
    -------
        subgraph2xyz:
            extracts xyz coordinates of a subgraph given a graph of a molecule        
        xyz2list:
            converts an xyz chemical file into a list with atoms and their coordinates     
        sdf2graph:
            converts a sdf chemical file into a graph
        molecule2graph:
            creates a graph using networkx based on the molecule attributes
        draw:
            plots a graph with networkx
        IntermolecularFragmentation:
            finds the fragments (molecules) of a chemical compound
        fragments2file:
            creates files for the fragments and all possible combinations of molecules from a chemical compound
        
    """
    
    molecule_number = -1

    def __init__(self, path, name = 'molecule'):
        """
        Parameters
        ----------
            path: string
                path indicating (sub)folders and chemical file name
            name: string, optional
                given name for a molecule/compound (default is 'molecule')
        """
        
        Molecule.molecule_number          += 1
        self.name                          = name + '_' + str(Molecule.molecule_number)
        self.path                          = path
        self.atoms, self.edges, self.bonds = self.sdf2graph()
        self.graph                         = self.molecule2graph()
        self.xyzlist                       = self.xyz2listconverter()
        self.inter_components              = self.IntermolecularFragmentation()
        self.intra_components              = self.degreesearch(1)
        self.intra_edges                   = self.leaves_connections()

# ~~~~~~~> First Part: Intermolecular Fragmentation <~~~~~~~~~
    @staticmethod
    def subgraph2xyz(xyz, subgraph):
        """ Gets a list of xyz coordinates and extracts the coordinates
            from a given subgraph
        
        Parameters
        ----------
        xyz: list
            xyz coordinates of a molecule/compound
        subgraph: list
            list of vertices in a subgraph
        
        Returns
        ----------
        list
            list of coordinates of a subgraph from a xyz list           
        """
        
        subxyz = []
        for vertex in subgraph:
            subxyz.append(xyz[vertex])
        return subxyz
    
    @staticmethod
    def xyz2list(path_file):
        """ Converts an xyz chemical file to a list with
            atoms names and their coordinates
        
        Parameters
        ----------
        path_file:
            gives the location of an xyz chemical file
        
        Returns
        ----------
        list
            list of atoms and xyz coordinates        
        """
        
        mol_handle = open(path_file)
        xyz_mol = mol_handle.readlines()
        xyz_mol.pop(0)
        xyz_mol.pop(0)
        xyz_mol.pop(0)
        xyz_mol.pop(0)

        return xyz_mol
  
    def sdf2graph(self):
        """ Extracts graph properties of a molecule/compound from an sdf file
        
        Returns
        ----------
        list, list, list
            atoms, edges and bonds as lists        
        """
        sdf_file = open(self.path)
        sdf_lines = sdf_file.readlines()
        num_atoms = int(sdf_lines[3][0:3])

        atoms = []
        for a in range(num_atoms):
            atoms.append((sdf_lines[4+a][30:33]).strip(' '))

        edges = []
        bonds = []
        pos = num_atoms + 4
        start_node = (sdf_lines[pos][0:3]).strip(' ')
        while(start_node != 'M'):
            end_node = (sdf_lines[pos][3:6]).strip(' ')
            edges.append( (int(start_node)-1, int(end_node)-1) )
            bonds.append( int((sdf_lines[pos][6:9]).strip(' ')))

            pos += 1
            start_node = (sdf_lines[pos][0:3]).strip(' ')

        return atoms, edges, bonds
    
    def molecule2graph(self):
        """ Gets atoms, edges, and bonds, and transforms into a graph
            using networkx object
        ----------
        
        Returns
        ----------
        networkx
            given a molecule, a graph in networkx        
        """
        G = nx.Graph()
        d = {}
        for x in range(len(self.atoms)):
            d[x] = self.atoms[x]
            G.add_node(x, color = node_color[self.atoms[x]])


        for x in self.edges:
            G.add_edges_from([x])

        return G

        
    def IntermolecularFragmentation(self):
        """ Searches for fragments as molecules in a chemical compound
            using connected_components function from networkx
        
        Returns
        ----------
        list
            fragments as lists in a chemical compound
        
        """
        graph_components = [list(x) for x in list(nx.connected_components(self.graph))]
        return graph_components
    

    def fragments2file(self,folder_name='fragmentation_folder'):
        """ Writes files for each fragments and for every combination
            of fragments using combinations function from itertools
        
        Returns
        ----------
        xyz
            files and folders for every fragment and their combinations
        
        """
        number_of_components = len(self.inter_components)
        xyz_mol = self.xyz2list(self.path)

        #--- creating the main directory
        p = Path(f'results/{folder_name}/') 
        try:
            
            p.mkdir()
        except FileExistsError as exc:
            pass

        #--- looping over every combination
        for c in range(number_of_components):
            frag_comb = list(combinations(np.arange(number_of_components), c + 1))

            p = Path(f'results/{folder_name}/combination{c+1}/') 
            try:
                p.mkdir()
            except FileExistsError as exc:
                pass


            for idx, value in enumerate(frag_comb):
                open(f'results/{folder_name}/combination{c+1}/compound_{idx}.xyz','w+')

                with open(f'results/{folder_name}/combination{c+1}/compound_{idx}.xyz','w+') as file:

                    for x in value:
                        molecule_content = self.subgraph2xyz(xyz_mol, self.inter_components[x])
                        [file.write(x) for x in molecule_content]

                        
# ~~~~~~~> Second Part: Intramolecular Fragmentation <~~~~~~~~~
    def degreesearch(self, d):
        G = self.graph
        intra_components = [x for x in G.nodes() if G.degree(x) == d]
        return intra_components
    
    def breakingleaves(self, numberOfLeaves = 1):
        leavesCombination = list(combinations(self.intra_components, numberOfLeaves))
        components = []
        
        for x in leavesCombination:
            list_aux = list(range(self.graph.number_of_nodes()))
            aux = []
            for y in list(x):
                list_aux.remove(y)
                aux.append([y])
                
            aux.append(list_aux)
            components.append(aux)
            
        return components
   
    def xyz2listconverter(self):
        """ Converts an xyz chemical file to a list with
            atoms names and their coordinates
        
        Parameters
        ----------
        path_file:
            gives the location of an xyz chemical file
        
        Returns
        ----------
        list
            list of atoms and xyz coordinates        
        """
        
        mol_handle = open(self.path)
        xyz_mol = mol_handle.readlines()
        xyz_mol.pop(0)
        xyz_mol.pop(0)
        xyz_mol.pop(0)
        xyz_mol.pop(0)
        
        xyz_mol = xyz_mol[0:self.graph.number_of_nodes()]
        
        xyz = []
        for x in xyz_mol:
            aux = x.split()
            xyz.append([aux[3], aux[0], aux[1], aux[2]])
            
        return xyz
    
    def bridgesearch(self):
        G = self.graph
        return list(nx.bridges(G))
 

    def components2xyz(self, components, folder_name='leaves_folder'):
        
        numberOfComponents = len(components)
        
        p = Path(f'results/{folder_name}/') 
        
       
        try:
            shutil.rmtree(p)
        except:
            pass
        
        try:
            p.mkdir()
        except FileExistsError as exc:
            pass
        
        for idx in range(numberOfComponents):
            p = Path(f'results/{folder_name}/component{idx}/') 
            try:
                p.mkdir()
            except FileExistsError as exc:
                pass
            
            for c in range(len(components[idx])):
                open(f'results/{folder_name}/component{idx}/molecule{c}.xyz','w+')
             
                with open(f'results/{folder_name}/component{idx}/molecule{c}.xyz', "w") as file:
                    file.write(f'{len(components[idx][c])}\n')
                    file.write('#\n')
                    [file.write(' '.join(map(str, self.xyzlist[x])) +'\n') for x in components[idx][c]]
        
        return numberOfComponents

    def breakingbridges(self, bridges, d = 1):
        components = []
        C = list(combinations(bridges,d))
        F = self.graph.copy()

        return C

        
      
    def leaves_connections(self):
        """ extracts the connection between a leaf and other vertex

        Returns
        ----------
        list
            returns a list of lists in which we have pairs of vertices,
            and one them is a leaf in the graph

        """
        G = self.graph
        leaves = self.intra_components
        leaves_edges = []
        for vertex in leaves:
            leaves_edges.append([vertex, list(G[vertex])[0]])

        return leaves_edges
   

    def molecule_leaf_remainder(self):
        """ This methods creates files and folders for the
            entire molecule, a leaf, and the rest
            
        Returns
        ----------
        files
            for each leaf we have three files in a specific folder,
            containing a file for the entire molecule, other file
            for the leaf, and another for the rest, deleting the leaf

        """
        mol_name = self.name
        leaves_g = self.intra_components
        xyz_mol  = self.xyz2list(self.path)

        p = Path(f'{mol_name}/') 
        try:
            p.mkdir()
        except FileExistsError as exc:
            pass


        for idx, leaf in enumerate(leaves_g):
            p = Path(f'{mol_name}/leaf{idx}/') 
            try:
                p.mkdir()
            except FileExistsError as exc:
                pass

            open(f'{mol_name}/leaf{idx}/leaf_molecule.xyz','w+')
            open(f'{mol_name}/leaf{idx}/molecule.xyz','w+')
            open(f'{mol_name}/leaf{idx}/remainder_molecule.xyz','w+')


        for idx, x in enumerate(leaves_g):
            #--- molecule
            with open(f'{mol_name}/leaf{idx}/molecule.xyz', "w") as file:
                file.write(f'{len(xyz_mol)}\n')
                file.write('#\n')
                [file.write(x) for x in xyz_mol]

            #----- leaf
            molecule_content = self.subgraph2xyz(xyz_mol, [leaves_g[idx]])

            with open(f'{mol_name}/leaf{idx}/leaf_molecule.xyz', "w") as file:
                file.write('1\n')
                file.write('#\n')
                [file.write(x) for x in molecule_content]

            #----- remainder
            not_leaf = list(range(len(xyz_mol)))
            not_leaf.remove(leaves_g[idx])

            molecule_content = self.subgraph2xyz(xyz_mol, not_leaf)
            with open(f'{mol_name}/leaf{idx}/remainder_molecule.xyz', "w") as file:
                file.write(f'{len(not_leaf)}\n')
                file.write('#\n')
                [file.write(x) for x in molecule_content]