import streamlit as st
import streamlit.components.v1 as components
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network
from functions import *
 
#Network(notebook=True)
st.title('Hello Pyvis')

# make Network show itself with repr_html
pathways_name=pd.read_csv("data/pathways.tsv", sep='\t')["pathway_name"]
last_selection=list(pathways_name.values)[25]
#def net_repr_html(self):
#  nodes, edges, height, width, options = self.get_network_data()
#  html = self.template.render(height=height, width=width, nodes=nodes, edges=edges, options=options)
#  return html

#Network._repr_html_ = net_repr_html
st.sidebar.title('Choose your favorite Graph')
option=st.sidebar.selectbox('select graph',pathways_name)
physics=st.sidebar.checkbox('add physics interactivity?')


skip_calcs=False
pathway_edges=read_pathway(option)
if (len(pathway_edges)==0):
    option=last_selection
    skip_calcs=True
else:
    skip_calcs=False
    last_selection=option
pathway_edges=read_pathway(option)
if skip_calcs==False:
    adj_matrix,nodes_renamed,inv_nodes_renamed=build_adj(pathway_edges)
    G = nx.from_numpy_matrix(adj_matrix)
    triad_cliques=get_triad(G)
    weighted_edges=calculate_weighted_edges(triad_cliques, adj_matrix,inv_nodes_renamed)
    to_remove=[]
    for x in weighted_edges.items():
        zeros=0
        ones=0
        minus=0
        for z in x[1]:
            if (z[1]==0):
                zeros+=1
            elif (z[1]==1):
                ones+=1
            else:
                minus+=1
        if (ones==0):
            if (minus==0):
                to_remove.append(x[0])
            else:
                if (zeros/minus>1):
                    to_remove.append(x[0])

    relabel={}
    for e,node in enumerate( G.nodes()):
        relabel[e]=str(inv_nodes_renamed[node])
    net=Network(height="700px",notebook=True,directed=True,width="100%", bgcolor='#222222', font_color='white')
    for i,node in relabel.items():
        net.add_node(str(node))

    for edge in pathway_edges.values:
            if(edge[2]==-1):
                net.add_edge(str(edge[0]), str(edge[1]), color="yellow")
            else:
                net.add_edge(str(edge[0]), str(edge[1]))
    for triad in triad_cliques:
        for i,x in enumerate(triad):
            for j,y in enumerate(triad):
                if ((str(inv_nodes_renamed[triad[i]])+","+str(inv_nodes_renamed[triad[j]])) in to_remove) or ((str(inv_nodes_renamed[triad[j]])+","+str(inv_nodes_renamed[triad[i]])) in to_remove):
                    color="red"
                    size=10
                else:
                    color="green"
                    size=3
                weight=pathway_edges[(pathway_edges[0]==inv_nodes_renamed[triad[i]]) & (pathway_edges[1]==inv_nodes_renamed[triad[j]])]
                if (weight.empty):
                    continue
                weight=int(weight[2].values)
                if (weight==1):
                    net.add_edge(str(inv_nodes_renamed[triad[i]]), str(inv_nodes_renamed[triad[j]]), color=color, width=size,title="Express")
                else:
                    net.add_edge(str(inv_nodes_renamed[triad[i]]), str(inv_nodes_renamed[triad[j]]), color=color, width=size,title="Suppress")
    net.hrepulsion(node_distance=120, central_gravity=0.0, spring_length=100, spring_strength=0, damping=0.09)
    net.show("data/graph.html")
HtmlFile = open("data/graph.html", 'r', encoding='utf-8')
source_code = HtmlFile.read() 
components.html(source_code, height = 900,width=900)

