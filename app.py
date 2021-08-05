import streamlit as st
import streamlit.components.v1 as components
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network
from functions import *
import PIL
#Network(notebook=True)

# make Network show itself with repr_html
pathways_name=pd.read_csv("data/pathways.tsv", sep='\t')["pathway_name"]
tmp=pathways_name.iloc[0]
pathways_name.iloc[0]=pathways_name.iloc[25]
pathways_name.iloc[25]=tmp
#def net_repr_html(self):
#  nodes, edges, height, width, options = self.get_network_data()
#  html = self.template.render(height=height, width=width, nodes=nodes, edges=edges, options=options)
#  return html
st.set_page_config(layout="wide")

#Network._repr_html_ = net_repr_html
st.sidebar.title('Choose a pathway')
option=st.sidebar.selectbox('',pathways_name)
st.sidebar.text("Edge legend:")
w=25
h=5
red = np.zeros((h, w, 3), dtype=np.uint8)
red[0:]=[255, 0, 0]
green = np.zeros((h, w, 3), dtype=np.uint8)
green[0:]=[0, 128, 0]
blue = np.zeros((h, w, 3), dtype=np.uint8)
blue[0:]=[30,144,255]
yellow = np.zeros((h, w, 3), dtype=np.uint8)
yellow[0:]=[255, 255,0]
st.sidebar.image(blue, caption='Expression edges')
st.sidebar.image(yellow, caption='Suppression edges')
st.sidebar.image(green, caption='Part of triad edges')
st.sidebar.image(red, caption='Removed edges')

skip_calcs=False
pathway_edges=read_pathway(option)
if (len(pathway_edges)==0):
    #option=last_selection
    skip_calcs=True
    st.error("Edges not found, try another pathway!")
else:
    skip_calcs=False
    #last_selection=option
    pathway_edges=read_pathway(option)
    adj_matrix,nodes_renamed,inv_nodes_renamed=build_adj(pathway_edges)
    G = nx.from_numpy_matrix(adj_matrix)
    triad_cliques=get_triad(G)
    weighted_edges=calculate_weighted_edges(triad_cliques, adj_matrix,inv_nodes_renamed)
    to_remove=[]
    signify_values={}
    essential_edges=[]
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
                m=(zeros+minus)/2
                if ((minus+zeros)/(zeros*minus+1)*zeros/(minus+zeros)>((minus+zeros)/(m*m+1))*zeros/(minus+zeros)):
                    to_remove.append(x[0])
        else:
            essential_edges.append(x[0])
        if (ones==0):
            signify_values[x[0]]=(minus+zeros)/(zeros*minus+1)*zeros/(minus+zeros)
        else:
            signify_values[x[0]]=0

    relabel={}
    for e,node in enumerate( G.nodes()):
        relabel[e]=str(inv_nodes_renamed[node])
    net=Network(height="825px",notebook=False,directed=True,width="1800px", bgcolor='#222222', font_color='white')
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
                value=""
                isessential=""
                min_node=min(triad[i],triad[j])
                max_node=max(triad[i],triad[j])
                weight=pathway_edges[(pathway_edges[0]==inv_nodes_renamed[min_node]) & (pathway_edges[1]==inv_nodes_renamed[max_node])]
                if (weight.empty):
                    continue
                if ((str(inv_nodes_renamed[min_node])+","+str(inv_nodes_renamed[max_node])) in to_remove):           
                    color="red"
                    size=10
                    value+=", significativity:  "+str(signify_values[str(inv_nodes_renamed[min_node])+","+str(inv_nodes_renamed[max_node])])
                else:
                    color="green"
                    size=3
                    value+=", significativity:  "+str(signify_values[str(inv_nodes_renamed[min_node])+","+str(inv_nodes_renamed[max_node])])
                if ((str(inv_nodes_renamed[min_node])+","+str(inv_nodes_renamed[max_node])) in essential_edges):   
                    isessential="Essential "
                weight=int(weight[2].values)
                if (weight==1):
                    net.add_edge(str(inv_nodes_renamed[min_node]), str(inv_nodes_renamed[max_node]), color=color, width=size,title=isessential+"Expression edge"+value)
                else:
                    net.add_edge(str(inv_nodes_renamed[min_node]), str(inv_nodes_renamed[max_node]), color=color, width=size,title=isessential+"Suppression edge"+value)
    net.hrepulsion(node_distance=120, central_gravity=0.0, spring_length=100, spring_strength=0, damping=0.09)
    net.show("data/graph.html")
    HtmlFile = open("data/graph.html", 'r', encoding='utf-8')
    source_code = HtmlFile.read() 
    components.html(source_code, height = 850,width=1850)

