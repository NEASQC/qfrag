import networkx as nx
import matplotlib.pyplot as plt

def nudge(pos, x_shift, y_shift):
    """ 
    Shifts the positions of the atom number in the graph, in order to plot them above the nodes
    
    Parameters
    ----------
    pos : dic(float)
        dictionary with atom number and x-y positions based on the graph structure
    x_shift : float
        amount of shift to perform along x-axis
    y_shift : float
        amount of shift to perform along y-axis

    Returns
    ----------
    D: dict
        dictionary with the shifted layout
        
    max(pos_x), min(pos_x), max(pos_y), min(pos_y) : float
        maximum an mininum values on each axis

    """
    D = {n:(x + x_shift, y + y_shift) for n,(x,y) in pos.items()}
    
    pos_x = [D[x][0] for x in D.keys()]
    pos_y = [D[x][1] for x in D.keys()]
    
    return D, max(pos_x), min(pos_x), max(pos_y), min(pos_y)


def draw(molecule):
    """ 
    Draws a graph using its properties
    
    Parameters
    ----------
    molecule : krachem molecule class
    
    Returns
    ----------
    plot : matplotlib
        draws a graph based on its parameters

    """
    d = {}
    for x in range(len(molecule.atoms)):
        d[x] = molecule.atoms[x]
        
    layout = nx.kamada_kawai_layout(molecule.graph)
    s_layout, max_x, min_x, max_y, min_y = nudge(layout, 0.05, 0.1)



    colors_n = nx.get_node_attributes(molecule.graph,'color').values()
    
    plt.figure(dpi = 70)
    nx.draw_kamada_kawai(molecule.graph,labels = d, node_color = colors_n,font_weight = 'bold')
    
    nx.draw_networkx_labels(molecule.graph, pos=s_layout)
    
   
    plt.ylim([min_y-.5, max_y*1.1])
    plt.xlim([min_x-.3, max_x*1.1])
    plt.show()

def fragmentdraw(molecule, subgraph):
    """ 
    Draws a subgraph using its properties
    
    Parameters
    ----------
    molecule : krachem molecule class
    
    Returns
    ----------
    plot : matplotlib
        draws a subgraph based on its parameters

    """
    G = molecule.graph.copy()
    S = G.subgraph(subgraph)

    d = {}
    for x in subgraph:
        d[x] = molecule.atoms[x]
        
    layout = nx.kamada_kawai_layout(S)
    s_layout, max_x, min_x, max_y, min_y = nudge(layout, 0.05, 0.1)
    
    colors_n = nx.get_node_attributes(S,'color').values() 
    nx.draw_kamada_kawai(S,labels = d, node_color = list(colors_n), font_weight = 'bold')
    nx.draw_networkx_labels(S, pos=s_layout)
   
    plt.ylim([min_y-.3, max_y*1.1])
    plt.xlim([min_x-.3, max_x*1.1])

    plt.show()
    
def breakingdraw(molecule, components):
    """ 
    Draws a graph and highlights a connection that is broken based on a selection of components
    
    Parameters
    ----------
    molecule : krachem molecule class
    components : list(int)
        position of atoms to be broken
    
    Returns
    ----------
    plot : matplotlib
        draws a graph based on its breaking bonds
    """
    G = molecule.graph.copy()
    edges = G.edges()
    d = {}
    for x in range(len(molecule.atoms)):
        d[x] = molecule.atoms[x]

    for e in edges:
        nx.set_edge_attributes(G, {e: {"color": 'k'}})
        nx.set_edge_attributes(G, {e: {"weight": 1.0}})

    for x in range(len(components) - 1):
        nx.set_edge_attributes(G, {list(G.edges(components[x]))[0]: {"color": 'r'}})
        nx.set_edge_attributes(G, {list(G.edges(components[x]))[0]: {"weight": 1.2}})

    layout = nx.kamada_kawai_layout(G)
    s_layout, max_x, min_x, max_y, min_y = nudge(layout, 0.05, 0.1)

    colors_n = nx.get_node_attributes(G,'color').values()   
    colors_e = nx.get_edge_attributes(G,'color').values()
    w = nx.get_edge_attributes(G,'weight').values()

    nx.draw_kamada_kawai(G,labels = d, node_color = list(colors_n), edge_color=colors_e, font_weight = 'bold',width=list(w))
    nx.draw_networkx_labels(G, pos=s_layout)
   
    plt.ylim([min_y-.3, max_y*1.1])
    plt.xlim([min_x-.3, max_x*1.1])

    plt.show()

def draw_intra(molecule):
    """ Draws a graph based on its leaves and the rest,
        represented by different colors
        
    Parameters
    ----------
    molecule : krachem molecule class
    
    Returns
    ----------
    plot
        plots a graph based on leaves in green and the rest in brown
    """

    atoms = molecule.atoms
    edges = molecule.edges
    G = nx.Graph()
    d= {}
    plt.figure(figsize=(4 + len(atoms)//10,4 + len(atoms)//10))

    for x in range(len(atoms)):
        d[x] = atoms[x]
        G.add_node(x, color = '#946c5c')     

    for x in range(len(edges)):
        G.add_edge(edges[x][0],edges[x][1], color = '#535c5b')     


    for x in range(len(atoms)):

        if(G.degree()[x]==1):
            G.add_node(x, color = '#61ab95')

    layout = nx.kamada_kawai_layout(G)
    s_layout, max_x, min_x, max_y, min_y = nudge(layout, 0.05, 0.1)
    
    colors_n = nx.get_node_attributes(G,'color').values()
    colors_e = nx.get_edge_attributes(G,'color').values()   

    nx.draw_kamada_kawai(G,labels = d, node_color = colors_n,font_weight = 'bold')
    nx.draw_networkx_labels(G, pos=s_layout)
   
    plt.ylim([min_y-.3, max_y*1.1])
    plt.xlim([min_x-.3, max_x*1.1])
    
    plt.show()