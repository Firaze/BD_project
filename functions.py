import semopy
from semopy import Model
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import networkx as nx
from matplotlib import pyplot as plt

pathways=pd.read_csv("data/pathways.tsv", sep='\t')
esets=pd.read_csv("data/controls_counts_norm.csv")
esets.index=esets["Unnamed: 0"].values
esets=esets.drop(columns="Unnamed: 0")
gene_edges=pd.read_csv("data/gene_edges.tsv", sep='\t')
#for e,x in enumerate(gene_edges.values):
#    if (x[2]==-1):
#        gene_edges.iloc[e]=x[1],x[0],1
def read_pathway(pathway_name):
    pathway_id=np.where(pathways["pathway_name"]==pathway_name)[0]
    apathway=pathways["nodes"].loc[int(pathway_id)]
    apathway=apathway.split(";")
    apathway=[int(x) for x in apathway]
    pathway_edges=pd.DataFrame([list(x) for x in gene_edges.values if x[0] in apathway])
    pathway_edges=pd.DataFrame([x for x in pathway_edges.values if str(x[0]) in esets.index])
    pathway_edges=pd.DataFrame([x for x in pathway_edges.values if str(x[1]) in esets.index])
    return pathway_edges

def build_adj(pathway_edges):
    pathway_edges_0=pathway_edges[0].unique()
    pathway_edges_1=pathway_edges[1].unique()
    nodes=list(np.hstack((pathway_edges_0,pathway_edges_1)))
    nodes_renamed={}
    inv_nodes_renamed={}
    for e,x in enumerate(nodes):
        nodes_renamed[x]=e
        inv_nodes_renamed[e]=x
    nodes=len(nodes)
    adj_matrix=np.zeros((nodes,nodes))
    for x in pathway_edges.values:
        adj_matrix[nodes_renamed[x[0]]][nodes_renamed[x[1]]]=x[2]
    return (adj_matrix,nodes_renamed,inv_nodes_renamed)


def print_triad():
    temp=nx.DiGraph()
    
    for i,x in enumerate({k:x for k,x in inv_nodes_renamed.items() if k in triad_cliques[triad]}.items()):
        relabel[i]=x[1]
    relabel

    G_triad=nx.from_numpy_matrix(triad_matrix,create_using=temp,parallel_edges=True)
    G_triad=nx.relabel_nodes(G_triad,relabel)
    #nx.draw(G_triad, with_labels=True, font_weight='bold')
    pos=nx.random_layout(G_triad)
    nx.draw_networkx_nodes(G_triad, pos, node_color = 'r', node_size = 200, alpha = 0)

    nx.draw_networkx_labels(G_triad, pos)

    ax = plt.gca()
    ax.set_facecolor("yellow")
    ax.set_facecolor('silver')
    plt.rcParams['axes.facecolor']='red'
    m.inspect()["Estimate"][0]
    for i,e in enumerate(G_triad.edges):
        print(e)
        ax.annotate("",
                    xy=pos[e[0]], xycoords='data',
                    xytext=pos[e[1]], textcoords="data",
                    arrowprops=dict(arrowstyle="->", color="0.2",
                                    shrinkA=20, shrinkB=20,
                                    patchA=None, patchB=None,
                                    ),
                    )
        text_pos=np.abs(pos[e[0]]+pos[e[1]])/2
        if (i==2):
            continue
        plt.text(x=text_pos[0],y=text_pos[1],s='%.3f'%(m.inspect()["Estimate"][i]))


    plt.axis('off')
    plt.show()
    
def get_triad(G):
    all_cliques= nx.enumerate_all_cliques(G)
    triad_cliques=[x for x in all_cliques if len(x)==3 ]
    return triad_cliques
    
def calculate_weighted_edges(triad_cliques, adj_matrix,inv_nodes_renamed):
    weighted_edges={}
    new_triad_cliques=[]
    first_label=""
    second_label=""
    third_label=""
    mod = """ 
              y ~ x1 + x2
              """
    for triad in range(len(triad_cliques)):
        triad_matrix=np.zeros((3,3))
        for i,x in enumerate(triad_cliques[triad]):
            for j,y in enumerate(triad_cliques[triad]):
                triad_matrix[i][j]=adj_matrix[x][y]
      #  if list(triad_matrix[0])== [0, 1, 1]:
        zeros_count=np.array([len(np.where(x==0)[0]) for i,x in enumerate(triad_matrix) ])
        if (sum(zeros_count)==6):
            new_triad_cliques.append(triad_cliques[triad])
        else:
            continue
        first_index=int(np.where(zeros_count==1)[0])
        second_index=int(np.where(zeros_count==2)[0])
        third_index=int(np.where(zeros_count==3)[0])
        first_label=str(inv_nodes_renamed[triad_cliques[triad][first_index]])
        second_label=str(inv_nodes_renamed[triad_cliques[triad][second_index]])
        third_label=str(inv_nodes_renamed[triad_cliques[triad][third_index]])
        first_gene=(list(esets.loc[first_label,:].values),0)
        second_gene=(list(esets.loc[second_label,:].values),1)
        third_gene=(list(esets.loc[third_label,:].values),2)
       # elif list(triad_matrix[0])== [0, -1, -1]:
          #  first_gene=(list(esets.loc[str(inv_nodes_renamed[triad_cliques[triad][2]]),:].values),2)
            #second_gene=(list(esets.loc[str(inv_nodes_renamed[triad_cliques[triad][1]]),:].values),1)
            #third_gene=(list(esets.loc[str(inv_nodes_renamed[triad_cliques[triad][0]]),:].values),0)
            #first_label=str(inv_nodes_renamed[triad_cliques[triad][2]])
            #second_label=str(inv_nodes_renamed[triad_cliques[triad][1]])
            #third_label=str(inv_nodes_renamed[triad_cliques[triad][0]])

        y=third_gene
        x1=first_gene
        x2=second_gene
        to_df={"y":y[0],"x1":x1[0],"x2":x2[0]}
        data=pd.DataFrame(to_df).replace(np.inf, np.nan).replace(-np.inf, np.nan).dropna()
        m = Model(mod)
        r = m.fit(data)
        fac_sum=np.abs(r.x[0]+r.x[1])
        if (np.abs(r.x[0])<fac_sum*0.1):
            if (first_label+","+third_label in weighted_edges):
                weighted_edges[first_label+","+third_label].append((r.x[0],0))

            else:
                weighted_edges[first_label+","+third_label]=[(r.x[0],0)]
            if (second_label+","+third_label in weighted_edges) :
                weighted_edges[second_label+","+third_label].append((r.x[1],1))
            else:
                weighted_edges[second_label+","+third_label]=[(r.x[1],1)]

            if(first_label+","+second_label in weighted_edges) :
                weighted_edges[first_label+","+second_label].append((r.x[1],1))
            else:
                weighted_edges[first_label+","+second_label]=[(r.x[1],1)]
        elif(np.abs(r.x[1])<fac_sum*0.1):
            if (first_label+","+third_label in weighted_edges):
                weighted_edges[first_label+","+third_label].append((r.x[0],1))

            else:
                weighted_edges[first_label+","+third_label]=[(r.x[0],1)]
            if (second_label+","+third_label in weighted_edges):  
                weighted_edges[second_label+","+third_label].append((r.x[1],0))
            else:
                weighted_edges[second_label+","+third_label]=[(r.x[1],0)]
            if(first_label+","+second_label in weighted_edges) :
                weighted_edges[first_label+","+second_label].append((r.x[1],0))
            else:
                weighted_edges[first_label+","+second_label]=[(r.x[1],0)]
        else:
            if (first_label+","+third_label in weighted_edges):
                weighted_edges[first_label+","+third_label].append((r.x[0],-1))

            else:
                weighted_edges[first_label+","+third_label]=[(r.x[0],-1)]
            if (second_label+","+third_label in weighted_edges):  
                weighted_edges[second_label+","+third_label].append((r.x[1],-1))
            else:
                weighted_edges[second_label+","+third_label]=[(r.x[1],-1)]
            if(first_label+","+second_label in weighted_edges) :
                weighted_edges[first_label+","+second_label].append((r.x[1],-1))
            else:
                weighted_edges[first_label+","+second_label]=[(r.x[1],-1)]
    return weighted_edges, new_triad_cliques



#import statsmodels.api as sm
#import statsmodels.genmod.families.links as links
#from statsmodels.stats.mediation import Mediation
#probit = links.probit
#outcome_model = sm.GLM.from_formula("y ~ x1 + x2",datas[3],family=sm.families.Binomial())
#mediator_model = sm.OLS.from_formula("x2 ~ x1", datas[3]) 
#med = Mediation(outcome_model, mediator_model,exposure="x1").fit()
#med.summary()