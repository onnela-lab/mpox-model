import networkx as nx
import random
import collections
import numpy as np
import pandas as pd
import itertools
import sys
import pickle
import os

import mpox_utils

from mpox_utils import *




### Define Functions

# This is a simplification of the the function that runs the simulation in mpox_utils. It is created to only add and 
#     update edges to save on computational time

def simulate_edges_only(seed, N, n_initial, steps, intervention_start = 70, behavior_change = 0, behavior_change_perc=0.5):
    
    np.random.seed(seed)
    random.seed(seed)
    
    ############ Random Numbers
    mu_e = 7.6 #from Charniga 2022
    sigma_e = 4.9 #from Charniga 2022
    clip_e = (0 - mu_e)/sigma_e # calculate to clip at 0 (can't be infectious before exposed)
    mu_i = 27 # assumption, consistent with others (Spicknall 2022, Ogoina 2023)
    sigma_i = 3 # assumption, consistent with others (Spicknall 2022, Ogoina 2023)
    clip_i = (0-mu_i)/sigma_i # calculate to clip at 0 (can't recover before infectious)
    min_u = 0
    max_u = 14
    pcontact_main = 0.22
    pcontact_casual = 0.14
    pcontact_onetime = 1
    
    onetime_prob = [0.065/7, 0.087/7, 0.086/7, 0.056/7, 0.055/7, 0.055/7]
    activity_prob = [0, 0.001, 0.0054, 0.0101, 0.0315, 0.286]
    ##########
    
    ### Create graph
    G, total_degseq, main_degseq, onetime_degseq, rel_activity, activity_strat = CreateGraph(N, onetime_prob, activity_prob)
    
    G, main_relationships, casual_relationships, onetime_relationships = init_relationships(G, onetime_degseq, main_degseq, rel_activity)
    
    #### Set up lists that track things between time steps
    # First year, no main/casual relationships should be expired (so no prior relationships to avoid)
    want_main = []
    want_casual = []
    ex_list_main = [None]*N
    ex_list_casual = [None]*N

    num_edges = [len(list(G.edges()))]
    num_onetime = [len(onetime_relationships)]
    num_contacts = np.zeros(steps)    
    
    # track all nodes
    
    # Determine which edges are present at time 0
    dat_t0 = set((x, y, z, tuple(d.items())) for (x, y, z, d) in list(G.edges(keys=True, data=True)))
    edges_dat = dat_t0.copy()
    
    remove_time = []
    
    for i in range(1,steps):   
         # allow for behavior change regarding onetime relationships (50% reduction)
        if i == intervention_start:
            if behavior_change == 1: # universal behavior change
                activity_prob = list(map(lambda x: x*behavior_change_perc, activity_prob))
                
            if behavior_change == 2: # targeted behavior change (only those in the top two strata of sexual activity)
                activity_prob = activity_prob[0:4] + list(map(lambda x: x*behavior_change_perc, activity_prob[4:]))
        
        ### Relationship Updates
        G, onetime_relationships, onetime_degseq_new = update_onetimerelationships(G, onetime_relationships, onetime_prob, rel_activity, activity_prob, activity_strat)
        G, want_main, want_casual, ex_list_main, ex_list_casual = update_relationships(G, want_main, want_casual, ex_list_main, ex_list_casual)
        
        
    
        # Add any new edges to a set of cumulative edges
        new_dat = set((x, y, z, tuple(d.items())) for (x, y, z, d) in G.edges(keys=True, data=True))

        # Identify edges removed in this step
        removed = set(edges_dat) - new_dat

        # Track removed edges with their absence duration
        remove_time_new = [(edge, 0) for edge in removed]
        remove_time.extend(remove_time_new)

        # Increment absence duration for all tracked removed edges
        remove_time = [(x[0], x[1] + 1) for x in remove_time]

        # Remove edges that have been absent for more than 7 steps
        fall_off = {x[0] for x in remove_time if x[1] > 7}

        # Filter out "fallen off" edges from cumulative edge set
        edges_dat = edges_dat | new_dat  # Update cumulative edge set with new edges
        edges_dat = {x for x in edges_dat if x not in fall_off}

        # Remove "fallen off" edges from the tracking list
        remove_time = [x for x in remove_time if x[0] not in fall_off]
        
        if i == 1:
            dat_t1 = edges_dat.copy()
        
        if i == 7:
            dat_t7 = edges_dat.copy()
            
        if i == 14:
            dat_t14 = edges_dat.copy()
            
        if i == 21:
            dat_t21 = edges_dat.copy()
            
        if i == 28:
            dat_t28 = edges_dat.copy()
            
        if i == 35:
            dat_t35 = edges_dat.copy()
            
        if i == 42:
            dat_t42 = edges_dat.copy()
            
        if i == 49:
            dat_t49 = edges_dat.copy()
            
        if i == 69:
            dat_t69 = edges_dat.copy()
            
        if i == 70:
            dat_t70 = edges_dat.copy()
            
        if i == 77:
            dat_t77 = edges_dat.copy()
            
        if i == 84:
            dat_t84 = edges_dat.copy()
        
        
    return(G,dat_t0, dat_t1, dat_t7, dat_t14, dat_t21, dat_t28, dat_t35, dat_t42, dat_t49 , dat_t69, dat_t70, dat_t77, dat_t84)



def get_measures_clustering(G):
    G2 = nx.Graph(G)
    degs = [d for n, d in G.degree()]
    return(nx.transitivity(G2), nx.average_clustering(G2), len(max(nx.connected_components(G2), key=len)), 
          np.mean(degs), np.max(degs))




####### Run 
# Get same list of random seeds each time
random.seed(1)
rand_seeds = random.sample(range(1000000), k = 1000)

# Get random seed for this array number
seed_num = int(sys.argv[1])
print("SEED")
print(seed_num)
seed = rand_seeds[seed_num]

# set seeds for simulation
random.seed(seed)
np.random.seed(seed)

date = sys.argv[2]

sim_string = 'clustering_SUPPLEMENT'

# Get other arguments
N = [5000, 10000, 20000, 40000, 80000]
p_infect = 0.9
steps = 85
intervention_start = 70
behavior_change = [0,2]

nrows = len(N) * len(behavior_change) * 10

results = np.zeros((nrows, 9))
curr_row = 0

for n in range(len(N)):
    for j in range(len(behavior_change)):
            
        # set seeds for simulation

        n_initial= int(N[n]/1000)

        (G, dat_t0, dat_t1, dat_t7, dat_t14, dat_t21, dat_t28, 
         dat_t35, dat_t42, dat_t49, dat_t69, dat_t70, dat_t77, dat_t84) = simulate_edges_only(seed, N[n], n_initial, steps, intervention_start, behavior_change[j])


        G0 = nx.MultiGraph()
        G0.add_nodes_from(list(G.nodes()))
        G0.add_edges_from((node1, node2, key, dict(attributes)) for node1, node2, key, attributes in dat_t0);

        G1 = nx.MultiGraph()
        G1.add_nodes_from(list(G.nodes()))
        G1.add_edges_from((node1, node2, key, dict(attributes)) for node1, node2, key, attributes in dat_t1);

        G7 = nx.MultiGraph()
        G7.add_nodes_from(list(G.nodes()))
        G7.add_edges_from((node1, node2, key, dict(attributes)) for node1, node2, key, attributes in dat_t7);

        G14 = nx.MultiGraph()
        G14.add_nodes_from(list(G.nodes()))
        G14.add_edges_from((node1, node2, key, dict(attributes)) for node1, node2, key, attributes in dat_t14);

        G28 = nx.MultiGraph()
        G28.add_nodes_from(list(G.nodes()))
        G28.add_edges_from((node1, node2, key, dict(attributes)) for node1, node2, key, attributes in dat_t28);

        G49 = nx.MultiGraph()
        G49.add_nodes_from(list(G.nodes()))
        G49.add_edges_from((node1, node2, key, dict(attributes)) for node1, node2, key, attributes in dat_t49);

        G69 = nx.MultiGraph()
        G69.add_nodes_from(list(G.nodes()))
        G69.add_edges_from((node1, node2, key, dict(attributes)) for node1, node2, key, attributes in dat_t69);

        G70 = nx.MultiGraph()
        G70.add_nodes_from(list(G.nodes()))
        G70.add_edges_from((node1, node2, key, dict(attributes)) for node1, node2, key, attributes in dat_t70);

        G77 = nx.MultiGraph()
        G77.add_nodes_from(list(G.nodes()))
        G77.add_edges_from((node1, node2, key, dict(attributes)) for node1, node2, key, attributes in dat_t77);

        G84 = nx.MultiGraph()
        G84.add_nodes_from(list(G.nodes()))
        G84.add_edges_from((node1, node2, key, dict(attributes)) for node1, node2, key, attributes in dat_t84);

        clustering_G0 = get_measures_clustering(G0)
        clustering_G1 = get_measures_clustering(G1)
        clustering_G7 = get_measures_clustering(G7)
        clustering_G14 = get_measures_clustering(G14)
        clustering_G28 = get_measures_clustering(G28)
        clustering_G49 = get_measures_clustering(G49)
        clustering_G69 = get_measures_clustering(G69)
        clustering_G70 = get_measures_clustering(G70)
        clustering_G77 = get_measures_clustering(G77)
        clustering_G84 = get_measures_clustering(G84)

        results[0+curr_row,:] = [N[n], behavior_change[j], seed_num, 0] + list(clustering_G0)
        results[1+curr_row,:] = [N[n], behavior_change[j], seed_num, 1] + list(clustering_G1)
        results[2+curr_row,:] = [N[n], behavior_change[j], seed_num, 7] + list(clustering_G7)
        results[3+curr_row,:] = [N[n], behavior_change[j], seed_num, 14] + list(clustering_G14)
        results[4+curr_row,:] = [N[n], behavior_change[j], seed_num, 28] + list(clustering_G28)
        results[5+curr_row,:] = [N[n], behavior_change[j], seed_num, 49] + list(clustering_G49)
        results[6+curr_row,:] = [N[n], behavior_change[j], seed_num, 69] + list(clustering_G69)
        results[7+curr_row,:] = [N[n], behavior_change[j], seed_num, 70] + list(clustering_G70)
        results[8+curr_row,:] = [N[n], behavior_change[j], seed_num, 77] + list(clustering_G77)
        results[9+curr_row,:] = [N[n], behavior_change[j], seed_num, 84] + list(clustering_G84)
        
        curr_row += 10

folder_path = str('output/' + sim_string + '/')

if not os.path.exists(folder_path):
    os.makedirs(folder_path)
        
with open(str(str(folder_path) + 'clustering_results_' + str(seed_num) + '_' + str(date) +'.pkl'), 'wb') as file: 
        # A new file will be created 
        pickle.dump(results, file) 