# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 10:55:20 2022

@author: Grégoire

Ce programme permet à partir d'un fichier fasta comprenant plusieurs séquences de réaliser un alignement Pairwise entre chacune d'elles puis de présenter le résultat sous forme de graphe .
Les nœuds représenteront les protéines avec leur identifiant comme label.
Le diamètre des nœuds Sera proportionnel à la moyenne des scores d'alignement de la protéine.
Si la moyenne des scores est négative alors le nœud aura une taille minimale . 
La couleur des nœuds reflète également la moyenne et permet de voir quels noeuds sont les mieux connectés. 
Les arêtes du graphe On pour label le score d'alignement des 2 nœuds, Et l'épaisseur est proportionnelle au score. Si celui-ci est négatif alors l'arrêté aura une épaisseur minimale.      

"""
import Bio
import networkx as  nx
from Bio import pairwise2
import re
import matplotlib.pyplot as plt
import os

print("-"*40)
print("Ce programme permet à partir d'un fichier fasta comprenant plusieurs séquences de réaliser un graph dont les noeuds sont les séquences,les arrêtes sont proportionelles au score en largeur, la couleur et la taille des noeuds dépendent de la moyenne des scores de la séquence.")
print("-"*40)

#Saisir le chemin du fichier et vérifier sa présence
present=False
while present==False:
    file=input("Saisir le chemin du fichier fasta:  ")
    if file in os.listdir("./"):
        present=True

#Saisir les paramètres d'alignement
param1=""
param2=""
penalize_mismatch=""
penalize_gap=""
score_match=""
score_mismatch=""
score_gap_open=""
score_gap_extend=""

#Afficher la question tant que la réponse n'est pas oui ou non
#Pénalité match/mismatch
while penalize_mismatch !="oui" and penalize_mismatch !="non":
    penalize_mismatch=input("Voulez vous personaliser les pénalités match/mismatch ? Répondre oui/non:  ")
if penalize_mismatch=="oui":
    param1="m"
    while True:
        try:
            score_match=int(input("Saisir le score de match (un entier positif): "))
        except:
            continue
        if score_match>0:
           break

    while True:
        try:
            score_mismatch=int(input("Saisir le score de mismatch (un entier négatif): "))
        except:
            continue
        if score_mismatch<0:
                break
#Pénalité gap
while penalize_gap !="oui" and penalize_gap !="non":
    penalize_gap=input("Voulez vous personaliser les pénalités de gap ? Répondre oui/non: ")

if penalize_gap=="oui":
    
    param2="s"
    
    while True:
        try:
            score_gap_open=int(input("Saisir la pénalité de gap d'ouverture (doit être négative): "))
        except:
            continue
        if score_gap_open<0:
                break
    
    while True:
        try:
            score_gap_extend=int(input("Saisir la pénalité de gap d'extension (elle doit être négative et moins forte que pour les gap d'ouverture):"))
        except:
            continue
        if score_gap_extend<0 and score_gap_extend>score_gap_open:
                break
            
#Fonction provenant du cours
# Fonction pour créer un dictionnaire avec comme clés les ID et valeurs les séquences à! partir d'un fasta
def getfasta(file):
    nameHandle = open(file,'r')
    fastas = {}
    for line in nameHandle:
        line=line.strip()
        if line[0]=='>':
            header = line[1:]
            fastas[header]=''
        else:
             fastas[header]+=line
    nameHandle.close()
    return fastas


dico_prot=getfasta(file)

nb_prot=len(dico_prot)

# Création du graph
G=nx.Graph()
for k in range(nb_prot):
    key=list(dico_prot.keys())[k]
    seq=dico_prot[key]
   
    for k1 in range(k+1,nb_prot):
        key1=list(dico_prot.keys())[k1]
        # Alignement Selon les paramètres 
        if param1=="m":
            if param2=="s":
                score = pairwise2.align.globalms(dico_prot[key],dico_prot[key1],score_match,score_mismatch,score_gap_open,score_gap_extend,score_only=True)
            else:
                 score = pairwise2.align.globalmx(dico_prot[key],dico_prot[key1],score_match,score_mismatch,score_only=True)
        else:
            if param2=="s":
                 score = pairwise2.align.globalxs(dico_prot[key],dico_prot[key1],score_gap_open,score_gap_extend,score_only=True)
            else:
                 score = pairwise2.align.globalxx(dico_prot[key],dico_prot[key1],score_only=True)
        G.add_edge(key,key1,weight=score)

#Dictionnaire contenant la moyenne des scores pour chaque protéine pour la taille des nœuds 
dico_score={}
nodes=list(G.nodes())
for n in range(len(nodes)):
    s=0
    for n1 in range(len(G.nodes)):
        if n!=n1:
            s+=G[nodes[n]][nodes[n1]]["weight"]
    dico_score[nodes[n]]=int(s)    
    
#Dictionnaire avec les scores pour de largeur des arrêtes 
weights = [G[u][v]['weight']/10 for u,v in G.edges()]
list_size=[max(s*100,50000) for s in dico_score.values()]
fig = plt.figure(1, figsize=(100, 50))
pos = nx.circular_layout(G)

#Affichage du graphe 
nodes=nx.draw_networkx_nodes(G,pos,node_size=list_size,cmap=plt.cm.Reds, node_color=list(dico_score.values()))
edges=nx.draw_networkx_edges(G,pos,width=[max(w,20) for w in weights],edge_color="blue")
edge_labels = nx.get_edge_attributes(G, "weight")
el=nx.draw_networkx_edge_labels(G,pos,edge_labels,font_size=50)
nl = nx.draw_networkx_labels(G,pos,bbox=dict(facecolor="white", edgecolor='blue', boxstyle='round,pad=0.2'),font_size=50)


plt.sci(nodes)
cbar=plt.colorbar()
cbar.ax.tick_params(labelsize=80)
cbar.ax.set_title(label="Score moyen",fontsize=75)
plt.sci(edges)
plt.title("Graph représentant les scores d'alignement entre chaque gène",fontsize=100)
plt.savefig("graph.png")

print("-"*40)
print("Le graphe 'graph.png' a été enregistré dans le répertoire courant")
print("-"*40)

