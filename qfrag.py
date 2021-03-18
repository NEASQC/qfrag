from openbabel import openbabel
import networkx	as nx
import numpy as np
import matplotlib.pyplot as plt
import math
from mendeleev import element

# --- Threshold for bonds
bond_thresh = 1.2

# --- Dictionary for covalent radii with mendeleev
cov_rads = {}
for x in range(1,119):
    cov_rads[element(x).symbol] = element(x).covalent_radius/100

# --- Color for the nodes of molecule graph
node_color = {'H':'#a1d4d6' ,'C':'#cfc0c3', 'O':'#ffa1b2', 
           'N':'#c487d6', 'P':'#d1bb64', 'Cl':'#6ccc87',
           'S':'#b87163', 'Na':'#e3c76d'}

# --- Color for the edges of molecule graphs based on the bonds
edge_color = {1:'k',2:'r',3:'b'}


# --- General converter between chemical molecular formats
def molecule_converter(molecule_input,input_format, output_format):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(input_format, output_format)
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, molecule_input)
    
    
    molecule_output = output_format + '/' + molecule_input.split('/')[1].split('.')[0] + '.' +output_format
    obConversion.WriteFile(mol, molecule_output)
    mol_file = open(molecule_output)
    
    return mol_file.read()

#--- SDF File to Graph Components: atoms (nodes), bonding (edges) and type of bonding
def sdf2graph(path: str):
    sdf_file = open(path)
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

#--- Graph Structures to Graph using Networkx
def graph2plot(atoms, edges, bonds):
    G = nx.Graph()
    d = {}
    for x in range(len(atoms)):
        d[x] = atoms[x]
        G.add_node(x, color = node_color[atoms[x]])
        
    for x in range(len(edges)):
        G.add_edge(edges[x][0],edges[x][1], color = edge_color[bonds[x]])     
    
    colors_n = nx.get_node_attributes(G,'color').values()
    colors_e = nx.get_edge_attributes(G,'color').values()
    nx.draw_kamada_kawai(G,labels = d, node_color = colors_n, edge_color = colors_e,font_weight = 'bold')
    plt.show()
    return G 

#---------------------------------------------------------------------------------
#-------------------------------XYZ CODE------------------------------------------
#---------------------------------------------------------------------------------
def get_file_string_array(file_name):
    try:
        file = open('molecules_xyz/' + file_name + '.xyz', "r")
    except IOError:
        print('Error: file (%s) not found!\n' % (file_name))
        sys.exit()
    lines = file.readlines()
    file.close()
    array = []
    for line in lines:
        array.append(line.split())
    return array

def get_geom(xyz_file_name):
    xyz_array = get_file_string_array(xyz_file_name)
    n_atoms = int(xyz_array[0][0])
    at_types = ['' for i in range(n_atoms)]
    coords = [[0.0 for j in range(3)] for i in range(n_atoms)]
    for i in range(n_atoms):
        at_types[i] = xyz_array[i+2][0]
        for j in range(3):
            coords[i][j] = float(xyz_array[i+2][j+1])
    geom = [at_types, coords]
    return geom

def get_r12(coords1, coords2):
    r2 = 0.0
    for p in range(3):
        r2 += (coords2[p] - coords1[p])**2
    r = math.sqrt(r2)
    return r

def get_bond_graph(geom):
    at_types, coords = geom[0:2]
    n_atoms = len(at_types)
    bond_graph = [[] for i in range(n_atoms)]
    for i in range(n_atoms):
        #covrad1 = (element(at_types[i]).covalent_radius)/100
        covrad1 = cov_rads[at_types[i]]
        for j in range(i+1, n_atoms):
            #covrad2 = (element(at_types[j]).covalent_radius)/100
            covrad2 = cov_rads[at_types[j]]
            thresh = bond_thresh * (covrad1 + covrad2)
            r12 = get_r12(coords[i], coords[j])
            if (r12 < thresh):
                bond_graph[i].append(j)
                bond_graph[j].append(i)
    return bond_graph

def mol2graph(mol_geom, connections):
    atoms = mol_geom[0]
    G = nx.Graph()
    for x in range(len(atoms)):
    	G.add_node(x, color = node_color[atoms[x]], 
                   pos = tuple(mol_geom[1][x][0:2]))
        

    for x in range(len(connections)):
        for y in range(len(connections[x])):
            G.add_edges_from([(x, connections[x][y])])
    return G	

'''
def graph2plot(mol_geom, G, draw = []):
    atoms = mol_geom[0]
    d = {}
    plt.figure(figsize=(4 + len(atoms)//10,4 + len(atoms)//10))
    for x in range(len(atoms)):
        d[x] = atoms[x]
        
    pos=nx.get_node_attributes(G,'pos')
    colors_n = nx.get_node_attributes(G,'color').values()
    if draw == 'kamada':
        nx.draw_kamada_kawai(G,node_color = colors_n,labels = d,font_weight = 'bold')
    else:
        nx.draw(G,pos,node_color = colors_n,labels = d,font_weight = 'bold')
    plt.show()
'''


