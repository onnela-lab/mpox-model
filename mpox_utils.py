## Modeling Monkeypox (utilities)

import networkx as nx
import random
import collections
import numpy as np
import scipy.stats as sp
import pandas as pd
import itertools
from operator import add

#---------------------------------#

# Import vaccine coverage data - per data pulled from CDC websites on 7/11/2023 (see paper supplement for details)

vax_cover = pd.read_excel("mpox_vax_coverage_data.xlsx", header = None)

# Get number of vaccines available per day (repeat each 7 times because original is by week)
daily_num_FD = np.rint(np.repeat(vax_cover.iloc[10, 1:].astype(float),7).tolist()).astype(int)
daily_num_SD = np.rint(np.repeat(vax_cover.iloc[11, 1:].astype(float),7).tolist()).astype(int)


#---------------------------------#

### Purpose: Set up graph with nodes and determine each node's relationship structure

### Inputs:
# N: number of nodes
# onetime_prob: probability of a onetime encounter, based on relationship type.
# activity_prob: probability of a onetime encounter, based on sexual activity stratum

### Outputs:
# G: graph
# total_degseq: total degree sequence (main+casual+onetime) . List of length N
# main_degseq: degree sequence for main and casual partnerships. List of length N
# onetime_degseq: degree sequence for onetime partnerships. List of length N
# rel_activity: each node's relationship type. List of length N, values from 0-4
# activity_strat each node's sexual activity stratum. List of length N, values from 0-5

def CreateGraph(N, onetime_prob, activity_prob):
    while True:
        # Set up empty lists (using -100 to easily spot if a node gets a bad value)
        main_degseq = [-100]*N
        onetime_degseq = [-100]*N
        #onetime_rand = [-100]*N
        
        # assign independently
        rel_activity = np.random.choice([0,1,2,3,4,5], N, p=[0.471, 0.167, 0.074, 0.22, 0.047, 0.021])
        main_degseq = [x if x in [0,1,2] else x-2 for x in rel_activity]
        
        activity_strat = np.random.choice([0,1,2,3,4,5], N, p=[0.19,0.19,0.19,0.19,0.19,0.05])
        strat_prob = [activity_prob[x] for x in activity_strat]
        onetime_degseq = [np.random.geometric(p= (1-(strat_prob[x])))-1 for x in activity_strat]
        
    
        total_degseq = list(map(lambda x,y:x+y, main_degseq, onetime_degseq))
        
        # Get the total number of main and casual stubs to ensure they're an even number
        want_main = len([x for x in rel_activity if x in (3,4,5)])
        want_casual = len([x for x in rel_activity if x in (1,4)]) + 2*len([x for x in rel_activity if x in (2,5)])
        
        # need an even number of main, casual, and onetime relationships
        if sum(onetime_degseq)%2 == 0 and want_main%2==0 and want_casual%2==0:
            break
            
    G = nx.MultiGraph()
    G.add_nodes_from(range(N))
    return(G, total_degseq, main_degseq, onetime_degseq, rel_activity, activity_strat)
    
#---------------------------------#

### Purpose: Initialize infection status

## Inputs:
# G: graph
# N: number of nodes
# n_initial: number of nodes infected at time 0
# activity_strat: list length N of each node's sexual activity stratum (initially want to infect highly active nodes)

## Outputs:
# status: N x 5 matrix, one column for S, I, E, and R status and a column indicating if they've been diagnosed
# infection_tracker: Nx4 matrix that tracks each node's infection source, infection time, recovery time, and via what type of relationship they were infected

def init_infection(G,N,n_initial,activity_strat):
    # Infection status
    high_activity_nodes = [x for x in list(G.nodes()) if activity_strat[x] in (4,5)]
    i_nodes = random.sample(high_activity_nodes, k=n_initial)
    s_nodes = list(set(G.nodes()) - set(i_nodes))
    e_nodes = []
    r_nodes = []
    
    status = np.zeros((N,5))
    status[s_nodes,0] = 1
    status[e_nodes,1] = 1
    status[i_nodes,2] = 1
    status[r_nodes,3] = 1
    # 5th column is diagnosed
    
    # keep track of infection source, infection time, and recovery time
    infection_tracker = np.zeros((N, 4))
    infection_tracker[i_nodes, 0] = -1
    
    return(status, infection_tracker)

#---------------------------------#

### Purpose: Randomly define how long each node is in the 'exposed' state (if it gets infected)

## Inputs:
# N: number of nodes
# mu_e: mean exposure time
# sigma_e: variance of exposure time
# seed: random seed

## Outputs:
# etimes: List of length N of exposure times in days

def get_etimes(N, mu_e, sigma_e, clip_e, seed=None):
    # Set up list of exposure times
    
    # Set random state
    r = np.random.RandomState(seed)
    # Draw probabilities from a truncated normal distribution
    etimes = np.array(np.round(sp.truncnorm.rvs(loc = mu_e, scale = sigma_e, a = clip_e, b = 100, size = N)))
    return(etimes)

#---------------------------------#
#---------------------------------#
### Purpose: Randomly define how long each node is in the 'infected' state (if it gets infected)

## Inputs:
# N: number of nodes
# mu_i: mean infection time
# sigma_i: variance of infection time
# seed: random seed

## Outputs:
# itimes: List of length N of infection times in days

def get_itimes(N, mu_i, sigma_i, clip_i, seed=None):
    # Set up list of infection times
    
    # Set random state
    r = np.random.RandomState(seed)
    # Draw probabilities from a truncated normal distribution
    itimes = np.array(np.round(sp.truncnorm.rvs(loc = mu_i, scale = sigma_i, a = clip_i, b = 100, size = N)))
    
    return(itimes)

#---------------------------------#

### Purpose: Keep track of how when a node's vaccination becomes effective

## Inputs:
# N: number of nodes
# fd_efftime: time until first dose becomes effective (begins at 14)
# sd_time: countdown 28 days until second dose can be recieved after first dose
# sd_efftime: countdown 14 days until second dose becomes effective

## Outputs:
# vtimes: List of length N of vaccination timings in days

def get_vtimes(N, fd_efftime, sd_time, sd_efftime):
    # create an array to allow count down for vaccination time 
    # count down 14 days for efficacy of first dose with a 28-day time to second dose, then additional 14 days to full vax
    # after receiving second dose
    
    # Draw probabilities from a normal distribution
    vtimes = np.repeat(np.array([[fd_efftime, sd_time, sd_efftime]]), repeats = N, axis = 0)
    
    return(vtimes)

#---------------------------------#

### Purpose: Initialize edges in the network

## Inputs:
# G: the network
# onetime_degseq = degree sequence of one-time relationships (list length N)
# main_degseq: degree seqence of the main and casual relationships (list length N)
# rel_activity: list of length N, each item is that node's relationship activity (of the 6 types)

## Outputs:
# G: the network
# main_relationships: list of main relationships
# casual_relationships: list of casual relationships
# onetime_relationship: list of one-time relationships
    
def init_relationships(G, onetime_degseq, main_degseq, rel_activity):
    nodes = list(G.nodes())
   
    ### One-time relationships
    onetime_stubs = []
    for node in range(len(onetime_degseq)):
        onetime_stubs.extend(onetime_degseq[node] * [node])
    random.shuffle(onetime_stubs)
    
    for i in range(0, len(onetime_stubs), 2):
        maxnode = max(onetime_stubs[i], onetime_stubs[i+1])
        minnode = min(onetime_stubs[i], onetime_stubs[i+1])
        k = G.new_edge_key(minnode,maxnode)
        G.add_edge(minnode,maxnode, k, rel_type = "Onetime", rel_duration = 0)
        
    ### Main relationships
    want_mains = [x for x in nodes if rel_activity[x] in (3,4,5)]
    random.shuffle(want_mains)
    for j in range(0, len(want_mains), 2):
        maxnode = max(want_mains[j], want_mains[j+1])
        minnode = min(want_mains[j], want_mains[j+1])
        k = G.new_edge_key(minnode,maxnode)
        G.add_edge(minnode,maxnode, k, rel_type = "Main", rel_duration = np.random.geometric(1/407))
        
    ### Casual relationships
    casual_stubs = []
    for node in range(len(main_degseq)):
        if rel_activity[node] in (1,2):
            casual_stubs.extend(main_degseq[node] * [node])
        elif rel_activity[node] in (4,5): #remove degree number for main relationship
            casual_stubs.extend((main_degseq[node]-1) * [node])
    random.shuffle(casual_stubs)
    
    for m in range(0, len(casual_stubs), 2):
        maxnode = max(casual_stubs[m], casual_stubs[m+1])
        minnode = min(casual_stubs[m],casual_stubs[m+1])
        k = G.new_edge_key(minnode,maxnode)
        G.add_edge(minnode,maxnode, k, rel_type = "Casual", rel_duration = np.random.geometric(1/166))


    # Output for bookkeeping
    onetime_relationships = [(u,v,k) for u,v,k,e in G.edges(keys=True,data=True) if e["rel_type"]== "Onetime"]
    main_relationships = [(u,v,k) for u,v,k,e in G.edges(keys=True,data=True) if e["rel_type"]=="Main"]
    casual_relationships = [(u,v,k) for u,v,k,e in G.edges(keys=True,data=True) if e["rel_type"]=="Casual"]
    
    return(G, main_relationships, casual_relationships, onetime_relationships)

#---------------------------------#

### Purpose: Update infection status each day

## Inputs:
# status: matrix of node infection status
# itimes: how many days each node is infectious before it recovers (list of length N)
# etimes: how many days each node is exposed before it becomes infected (list of length N)
# treatment_delay: allows for treatment delay in diagnosis
# infection_tracker: Nx4 matrix that tracks each node's infection source, infection time, recovery time, and via what type of relationship they were infected

## Outputs:
# Updated status and infection_tracker

def update_status(status, itimes, etimes, treatment_delay, infection_tracker, step):
    
    # Get new infections by exposure time and new recoveries by infection time
    s_nodes = np.where(status[:,0]==1)[0].tolist()
    e_nodes = np.where(status[:,1]==1)[0].tolist()
    i_nodes = np.where(status[:,2]==1)[0].tolist()
    r_nodes = np.where(status[:,3]==1)[0].tolist()
    
    # update time remaining
    for x in e_nodes:
        etimes[x] -= 1        
    for y in i_nodes:
        itimes[y] -= 1
        treatment_delay[y] -= 1
    
    new_infect = [x for x in e_nodes if etimes[x] < 1]
    new_recover = [x for x in i_nodes if itimes[x] < 1]
    new_diagnosed = [x for x in i_nodes if treatment_delay[x] < 1]
    
    # update status tracker
    status[new_infect,1] = 0
    status[new_infect,2] = 1
    status[new_recover,2] = 0
    status[new_recover,3] = 1
    status[new_diagnosed,4] = 1
    
    infection_tracker[new_recover, 2] = step
    
    return(status, infection_tracker)

#---------------------------------#

### Purpose: Update onetime relationships each day

## Inputs:
# G: graph
# onetime_relationships: list of onetime relationships for that day
# onetime_prob: list of probabilities of onetime encounters for each relationship style
# rel_activity: list of each node's relationship style. List of length N
# activity_prob: list of probabilities of onetime encounters for each sexual activity stratum
# activity_strat: list of each node's activity stratum. List of length N

## Outputs:
# G: graph
# onetime_relationships: list of onetime relationships
# onetime_degseq: degree sequence for one-time relationships

def update_onetimerelationships(G, onetime_relationships, onetime_prob, rel_activity, activity_prob, activity_strat):
    G.remove_edges_from(onetime_relationships)
    onetime_degseq = [-100] * len(list(G.nodes()))
            
    while True:
        for n in list(G.nodes()):
            prob_index = rel_activity[n]
            prob = onetime_prob[prob_index]
            
            strat_index = activity_strat[n]
            strat_prob = activity_prob[strat_index]
            
            onetime_degseq[n] = np.random.geometric(p= (1-(strat_prob)))-1
            
    # need an even number of onetime relationship stubs
        if sum(onetime_degseq)%2 == 0 :
            break
            
     ### One-times
    onetime_stubs = []
    for node in range(len(onetime_degseq)):
        onetime_stubs.extend(onetime_degseq[node] * [node])
    random.shuffle(onetime_stubs)
    
    for i in range(0, len(onetime_stubs), 2):
        maxnode = max(onetime_stubs[i], onetime_stubs[i+1])
        minnode = min(onetime_stubs[i], onetime_stubs[i+1])
        k = G.new_edge_key(minnode,maxnode)
        G.add_edge(minnode,maxnode, k, rel_type = "Onetime", rel_duration = 0)
        
    onetime_relationships = [(u,v,k) for u,v,k,e in G.edges(keys=True,data=True) if e["rel_type"]== "Onetime"]
        
    return(G, onetime_relationships, onetime_degseq)
    
#---------------------------------#

### Purpose: Update main and casual relationships

## Inputs: 
# G: graph
# want_main: nodes that want a main partnership that currently don't have one
# want_casual: nodes that want another casual partnership
# ex_list_main: each node's most recent 'ex-main partner' (so that they aren't immediately re-paired with that node)
# ex_list_casual: same as above, but for casual partners

## Outputs:
# updated versions of the inputs

def update_relationships(G, want_main, want_casual, ex_list_main, ex_list_casual):      
    removed_main = []
    removed_casual = []
    
    # define edges here so that it doesn't change as edges are removed below
    edges = list(G.edges(keys=True))

    ### Step 1: Remove expired relationships
    for e in edges:
        # Update relationship duration
        if G.edges(keys=True)[e]["rel_type"] != "Onetime":
            G.edges(keys=True)[e]["rel_duration"] -= 1
            
            # If relationship is over, remove it from the graph
            if G.edges(keys=True)[e]["rel_duration"] <= 0:
                if G.edges(keys=True)[e]["rel_type"] == "Main":
                    removed_main.append(e)
                elif G.edges(keys=True)[e]["rel_type"] == "Casual":
                    removed_casual.append(e)
                    
                u,v,k = e
                G.remove_edge(u,v,k)
    
    # Get all the nodes that lost an edge (and thus will want a new one)
    for e in removed_main:
        u,v,k = e
        want_main.append(u)
        want_main.append(v)
        # keep track of the most recent partner they had to avoid immediate re-pairing
        ex_list_main[u] = v
        ex_list_main[v] = u
        
    for e in removed_casual:
        u,v,k = e
        want_casual.append(u)
        want_casual.append(v)
        # keep track of the most recent partner they had to avoid immediate re-pairing
        ex_list_casual[u] = v
        ex_list_casual[v] = u
                           
    
    ### Step 2: form new relationships
    new_edge = 0
                
    # Make new main partnerships
    if len(want_main) > 0:
        want_main_copy = want_main
        skipped_main = []
        while len(want_main_copy)>0:
            node = want_main_copy[0]
            # Get a list of potential partners
            poss = [n for n in want_main_copy if (n!=ex_list_main[node] and n!= node)]
            if len(poss)>0:
                chosen = random.choice(poss)
                maxnode = max(node, chosen)
                minnode = min(node, chosen)
                k = G.new_edge_key(minnode,maxnode)
                G.add_edge(minnode,maxnode, k, rel_type = "Main", rel_duration = np.random.geometric(1/407))
                want_main_copy.remove(maxnode)
                want_main_copy.remove(minnode)
                new_edge +=1
                
            else:
                # If there aren't any potential partners, skip this node for now
                want_main_copy.remove(node)
                skipped_main.append(node)
                continue
        # Keep the skipped nodes on the list of who wants a main partnership
        want_main = skipped_main        
                
    # Make new casual partnerships
    if len(want_casual) > 0:
        want_casual_copy = want_casual
        skipped = []
        while len(want_casual_copy)>0:
            node = want_casual_copy[0]
            # Get a list of potential partners
            poss = [n for n in want_casual_copy if (n!=ex_list_casual[node] and n!= node)]
            
            if len(poss)>0:
                chosen = random.choice(poss)
                maxnode = max(node, chosen)
                minnode = min(node, chosen)
                k = G.new_edge_key(minnode,maxnode)
                G.add_edge(minnode,maxnode, k, rel_type = "Casual", rel_duration = np.random.geometric(1/166))
                want_casual_copy.remove(maxnode)
                want_casual_copy.remove(minnode)
                new_edge +=1
                
            else:
                # If there aren't any potential partners, skip this node for now
                want_casual_copy.remove(node)
                skipped.append(node)
                continue
        # Keep the skipped nodes on the list of who wants a casual partnership
        want_casual = skipped        
                
    return(G, want_main, want_casual, ex_list_main, ex_list_casual)
                
    
#---------------------------------#

### Purpose: MPX Spread

## Inputs:
# G: network
# step: current day in the simulation
# p_infect: probability of infection
# status: array of infection status (N x 4)
# pcontact_main: probability of sexual contact for main partners
# pcontact_casual: probability of sexual contact for casual partners
# pcontact_onetime: probability of sexual contact for one-time partners
# isolation_complier: list of length N of whether nodes comply fully with isolation or partially
# treatment_delay: list of length N for days until node is diagnosed
# vax_eff: vaccine efficacy
# infection_tracker: Nx4 matrix that tracks each node's infection source, infection time, recovery time, and via what type of relationship they were infected

## Outputs:
# G: the network
# status: array of infection status (N x 4)
# total_contacts: count of sexual contacts that day
# infection_tracker: Nx4 matrix that tracks each node's infection source, infection time, recovery time, and via what type of relationship they were infected

def spread(G, N, step, p_infect, status, pcontact_main, pcontact_casual, pcontact_onetime, isolation_complier, 
           treatment_delay, vax_eff, infection_tracker):
    i_nodes = np.where(status[:,2]==1)[0].tolist()
    total_contacts = 0
    infection_count_list = [0]*N
    
    for i in range(len(i_nodes)):
        node = i_nodes[i]
        # calculate number of new infections per node
        infection_count = 0
        
        # Factor in diagnosis for isolation
        if (isolation_complier[node] == 1) & (status[node,4] == 1):
            # partial compliance
            pm = 0.11
            pc = 0.07
            po = 0
        elif (isolation_complier[node] == 2) & (status[node,4] == 1):
            # full compliance
            pm = 0
            pc = 0
            po = 0
        else:
            pm = pcontact_main 
            pc = pcontact_casual
            po = pcontact_onetime
            
            
        # Get susceptible neightbors
        neighbors = [x for x in list(G.neighbors(node)) if status[x,0] == 1]
        
        # Determine if there is contact and transmission for each neighbor
        if len(neighbors)>0:
            random.shuffle(neighbors)
            for n in neighbors:

                min_node = min(node,n)
                max_node = max(node,n)
                keylist = [key for w,v,key in G.edges(keys=True) if (w == min_node and v == max_node)]

                for k in range(len(keylist)):
                # Get the node order (min,max) to get the right edge

                    key = keylist[k]

                    # Determine if there is contact that day
                    if ((G.edges[min_node,max_node,key]['rel_type']=='Main' and random.random() < pm) 
                        or (G.edges[min_node,max_node,key]['rel_type']=='Casual' and random.random() < pc) 
                        or (G.edges[min_node,max_node,key]['rel_type']=='Onetime' and po == 1)):
                        contact = 1
                        total_contacts += 1

                    else:
                        contact = 0

                    # If infected, set status to P1
                    if contact ==1 and random.random() < p_infect * (1-vax_eff[n]): #multiply by vax effect of the target
                        status[n,0] = 0
                        status[n,1] = 1 # Set to exposed
                        infection_tracker[n, 0] = node # define who their infection source was
                        infection_tracker[n, 1] = step # track what day they were infected
                        treatment_delay[n] = treatment_seeking(step)

                        # track type of edge
                        if G.edges[min_node,max_node,key]['rel_type']=='Main':
                            infection_tracker[n,3] = 1
                        elif G.edges[min_node,max_node,key]['rel_type']=='Casual':
                            infection_tracker[n,3] = 2
                        else:
                            infection_tracker[n,3] = 3
                
    return(G, status, total_contacts, infection_tracker)
            
#---------------------------------#
### Purpose: keep track of the current delay between infection and treatment seeking (when a node gets 'diagnosed' and begins isolating)

## Inputs:
# step: current step of simulation

## Outputs:
# delay: current treatment delay
    
def treatment_seeking(step):
    if step == 1:
        delay = 15
    elif step > 42:
        delay = 5
    else:
        delay = 15 - (10/42)*step
        
    return(delay)

#---------------------------------#

### Purpose: allow for node vaccination

## Inputs:
# G: the network
# activity_strat: list length N of each node's sexual activity stratum
# rel_activity: , 
# daily_num_FD: number of first dose vaccines available
# daily_num_SD: number of second dose vaccines available
# fd_eff: efficacy of a first dose
# sd_eff: efficacy of a second dose
# vax_stat: list of length N indicating how many vaccine doses a node has recieved
# vtimes: N x 3 array of vaccine timings
# step: day of the simulation
# vax_delay: day vaccination starts
# vax_scenario: vaccination scenario. Scenario 1 is checked outside the function and defaults to everyone who is sexually active can be vaccinated,
#               scenario 2 only vaccinates those in the highest 2 sexual activity strata

## Outputs:
# updated versions of vax_stat, vtimes, vax_eff

def vaccinate(G, activity_strat, rel_activity, daily_num_FD, daily_num_SD, fd_eff, sd_eff,
              vax_stat, vtimes, step, vax_delay, vax_scenario):
    # allow for the delay in vaccination
    vax_step = step-vax_delay
    
    ## Update Vaccination Times
    # For anyone vaccinated, update their vaccination time scheme (will not update anything until vaccines started)
    fd_recieved = np.where(vax_stat==1)[0].tolist()
    sd_recieved = np.where(vax_stat==2)[0].tolist()
    
    # update time til efficacy/second dose as appropriate
    if len(fd_recieved) > 0:
        vtimes[fd_recieved,0] -= 1
        vtimes[fd_recieved,1] -= 1
    if len(sd_recieved) > 0:
        vtimes[sd_recieved,2] -=1
    
    ## update vaccine efficacy
    vax_eff = [fd_eff if (vtimes[x,0] <= 0 and vtimes[x,2] > 0) else sd_eff if vtimes[x,2]<=0 else 0 for x in list(G.nodes)]
    
    ## Allocate new vaccinations
    # identify nodes that have yet to be vaccinated
    fd_eligible = np.where(vax_stat==0)[0].tolist()
       
    # identify nodes for a second dose
    sd_eligible = np.where(vtimes[:,1] <= 0)[0].tolist() #check that enough time has passed since first dose
    
    # restrict to to those who are sexually active
    fd_eligible = [x for x in fd_eligible if (activity_strat[x] != 0) & (rel_activity[x] != 0)]
    sd_eligible = [x for x in sd_eligible if (activity_strat[x] != 0) & (rel_activity[x] != 0)]
    
    # allow changes based on vaccination scenario - default behavior (scenario 0) is no vaccination,
    # scenario 1 is checked outside the function and defaults to everyone who is sexually active can be vaccinated
    if vax_scenario == 2 & vax_step <=56:
        if vax_step <= 28: # restrict to top two activity strata for first 4 weeks
            fd_eligible = [x for x in fd_eligible if activity_strat[x] in (4,5)]
            sd_eligible = [x for x in sd_eligible if activity_strat[x] in (4,5)]
        else: #expand to top four activity strata for the next four weeks
            fd_eligible = [x for x in fd_eligible if activity_strat[x] in (2,3,4,5)]
            sd_eligible = [x for x in sd_eligible if activity_strat[x] in (2,3,4,5)]

    # Allocate available vaccines to those eligible
    ## first dose 
    
    if len(fd_eligible) < daily_num_FD[vax_step]:
       
        vax_fd = fd_eligible
        daily_num_FD[vax_step+1] += (daily_num_FD[vax_step] - len(fd_eligible))
        
    else:
        vax_fd = random.sample(fd_eligible, daily_num_FD[vax_step])
    
    if len(vax_fd) > 0:
        vax_stat[vax_fd] = 1
    
    ## second dose
    if len(sd_eligible) < daily_num_SD[vax_step]:
        
        vax_sd = sd_eligible
        daily_num_SD[vax_step+1] += (daily_num_SD[vax_step] - len(sd_eligible))
        
    else: 
        vax_sd = random.sample(sd_eligible, daily_num_SD[vax_step])
    
    if len(vax_sd) > 0:
        vax_stat[vax_sd] = 2
                
    return(vax_stat, vtimes, vax_eff)

#---------------------------------#

### Purpose: run the simulation

## Inputs:
# seed: set a random seed
# N: Number of nodes in the network
# n_initial: number of initial infections
# p_infect: probability of infection
# steps: number of days for the simulation to run
# intervention_start: day the interventions start
# behavior_change: which scenario of behavior change to use (1 = universal, 2 = just the top two strata of sexual activity, 3 = a proportion of everyone)
# isolation = which isolation scenario to use (1 = everyone diagnosed isolation, 2 = only some isolation)
# behavior_change_perc: percent change in probability of having a one-time partner
# vax_scenario: vaccination scenario to use (0 = no vaccination, 1 = use vaccination)
# vax_delay: days to delay vaccination start
# daily_num_FD: daily number of first vaccine doses available
# daily_num_SD: daily number of second vaccine doses available
# vax_inc = multiplier to increase the number of vaccines available

## Outputs:
# number_E: list of length(steps) of number of nodes exposed per day
# number_I: list of length(steps) of number of nodes infected per day
# number_R: list of length(steps) of number of nodes recovered per day
# infection_tracker: Nx4 matrix that tracks each node's infection source, infection time, recovery time, and via what type of relationship they were infected

def simulate(N, n_initial, p_infect, steps, intervention_start, behavior_change, 
             isolation, behavior_change_perc, vax_scenario, vax_delay, daily_num_FD, daily_num_SD, vax_inc = 1, fd_eff = 0.36, sd_eff = 0.66):
    
    ############ Random Numbers
    # numbers for exposed state
    mu_e = 7.6 #from Charniga 2022
    sigma_e = 4.9 #from Charniga 2022
    clip_e = (0 - mu_e)/sigma_e # calculate to clip at 0 (can't be infectious before exposed)
    # numbers for infectious state
    mu_i = 27 # assumption, consistent with others (Spicknall 2022, Ogoina 2023)
    sigma_i = 3 # assumption, consistent with others (Spicknall 2022, Ogoina 2023)
    clip_i = (0-mu_i)/sigma_i # calculate to clip at 0 (can't recover before infectious)
    min_u = 0
    max_u = 14
    pcontact_main = 0.22
    pcontact_casual = 0.14
    pcontact_onetime = 1
    
    # vax efficacy per Deputy 2023
    fd_efftime = 14
    sd_efftime = 14
    sd_time = 28
    
    onetime_prob = [0.065/7, 0.087/7, 0.086/7, 0.056/7, 0.055/7, 0.055/7]
    activity_prob = [0, 0.001, 0.0054, 0.0101, 0.0315, 0.286]
    ##########
    
    ### Create graph
    G, total_degseq, main_degseq, onetime_degseq, rel_activity, activity_strat = CreateGraph(N, onetime_prob, activity_prob)
    
    ever_onetime = onetime_degseq
    
    #### Initialize graph parameters    
    isolation_complier = [0]*N
    status, infection_tracker = init_infection(G, N, n_initial, rel_activity)
    
    etimes = get_etimes(N, mu_e, sigma_e, clip_e, seed=None)
    itimes = get_itimes(N, mu_i, sigma_i, clip_i, seed=None)
    vtimes = get_vtimes(N, fd_efftime, sd_time, sd_efftime)
    vax_stat = np.asarray([0]*N)
    vax_eff = [0]*N
    # Update Vax availability
    daily_num_FD = np.rint(daily_num_FD*vax_inc).astype(int)
    daily_num_SD = np.rint(daily_num_SD*vax_inc).astype(int)
    
    
    # if we start vaccinating before the disease spreads, add together the vaccines and distribute
    if vax_delay < 0:
        vax_step = abs(vax_delay)
        fd_avail = np.sum(daily_num_FD[:vax_step])
        sd_avail = np.sum(daily_num_SD[:vax_step])
        
        # make sure that I don't run out of vaccination days if vaccination starts early enough
        print("number of original vax steps: ", len(daily_num_FD))
        print("steps needed: ", abs(vax_delay)+steps)
        if (abs(vax_delay) + steps + 2) >= len(daily_num_FD):
            diff = (abs(vax_delay) + steps) - len(daily_num_FD) + 1
            daily_num_FD = np.pad(daily_num_FD, (0, diff), mode='constant', constant_values=0)
            daily_num_SD = np.pad(daily_num_SD, (0, diff), mode='constant', constant_values=0)
            print("number of updated vax steps: ", len(daily_num_FD))
        
        #determine eligibility
        fd_eligible = np.where(vax_stat==0)[0].tolist()
        fd_eligible = [x for x in fd_eligible if (activity_strat[x] != 0) & (rel_activity[x] != 0)]
        if vax_scenario == 2 :
            fd_eligible = [x for x in fd_eligible if activity_strat[x] in (4,5)]
        
        # in limited vaccination scenario, may have fewer people eligible at first than total vaccinations available
        if len(fd_eligible) < fd_avail:
            vax_fd = fd_eligible
            daily_num_FD[vax_step+1] += (fd_avail - len(fd_eligible))
        
        else:
            vax_fd = random.sample(fd_eligible,fd_avail)
            
        # distribute the days that they got vaccinated
        if len(vax_fd) > 0:
            vax_stat[vax_fd] = 1
            for i in range(vax_step):
                if len(vax_fd) < daily_num_FD[i]:
                    sel = vax_fd
                else:
                    sel = random.sample(vax_fd, daily_num_FD[i])
                    
                vtimes[sel,0] = fd_efftime - (vax_step - i)
                vtimes[sel,1] = sd_time - (vax_step - i)
                vax_fd = [x for x in vax_fd if x not in sel]
                if len(vax_fd) == 0:
                    continue
             
        # identify nodes for a second dose
        sd_eligible = np.where(vtimes[:,1] <= 0)[0].tolist() 
        sd_eligible = [x for x in sd_eligible if (activity_strat[x] != 0) & (rel_activity[x] != 0)]
        if vax_scenario == 2 :
            sd_eligible = [x for x in sd_eligible if activity_strat[x] in (4,5)]
        
        if len(sd_eligible) < sd_avail:
            vax_sd = sd_eligible
            daily_num_SD[vax_step+1] += (sd_avail - len(sd_eligible))
        
        else:
            vax_sd = random.sample(sd_eligible,sd_avail)

        if len(vax_sd) > 0:
            vax_stat[vax_sd] = 2
            for i in range(vax_step):
                if len(vax_sd) < daily_num_SD[i]:
                    sel_sd = vax_sd
                else:
                    sel_sd = random.sample(vax_sd, daily_num_SD[i])
                
                vtimes[sel_sd,2] = sd_efftime - (vax_step - i)
                vax_sd = [x for x in vax_sd if x not in sel_sd]       
        

    
    # Whether people will seek treatment
    seek_treatment = np.random.binomial(n=1,p=0.8,size = N).tolist()
    treatment_delay = [10000]*N #initialize that no one will get treatment within 10,000 days, changes once a node is infected
    
    # initialize individual isolation behavior - if they don't seek treatment, they won't isolate
    # 2 =  full compliance, 1 = partial compliance, 0 = no compliance
    if isolation == 1:
        # Everyone diagnosed complies with isolation
        isolation_complier = [x * 2 for x in seek_treatment]
    elif isolation == 2:
        qc = list(np.random.choice([0,1,2], size = N))
        isolation_complier = [x*y for x,y in zip(seek_treatment, qc)]
    else:
        isolation_complier = [0]*N
    
    G, main_relationships, casual_relationships, onetime_relationships = init_relationships(G, onetime_degseq, main_degseq, rel_activity)
    
    #### Set up lists that track things between time steps
    # First year, no main/casual relationships should be expired (so no prior relationships to avoid)
    want_main = []
    want_casual = []
    ex_list_main = [None]*N
    ex_list_casual = [None]*N
    
    ### Begin accounting
    number_S = [np.sum(status[:,0])]
    number_E = [np.sum(status[:,1])]
    number_I = [np.sum(status[:,2])]
    number_R = [np.sum(status[:,3])]
    num_edges = [len(list(G.edges()))]
    num_onetime = [len(onetime_relationships)]
    num_contacts = np.zeros(steps)    
    
    for i in range(1,steps):   
            
        # If there are no more exposed/infected nodes, the simulation can halt
        if np.sum(status[:,[1,2]]) == 0:
            #print("No more exposed/infected nodes")
            break
            
        # allow for behavior change regarding onetime relationships (50% reduction)
        if i == intervention_start:
            if behavior_change == 1: # universal behavior change
                activity_prob = list(map(lambda x: x*behavior_change_perc, activity_prob))
                
            if behavior_change == 2: # targeted behavior change (only those in the top two strata of sexual activity)
                activity_prob = activity_prob[0:4] + list(map(lambda x: x*behavior_change_perc, activity_prob[4:]))
               
            if behavior_change == 3: # some percentage of everyone changes behavior (results not presented in the paper)
                temp = np.add(np.random.choice([0,-2], N, p=[1-behavior_change_perc,behavior_change_perc]),activity_strat)
                activity_strat = np.where(temp<0, temp*0, temp)
                
            
                           
            
        ### Infection Updates
        # Update the parameters with new year - move people from E -> I and I -> R
        status, infection_tracker = update_status(status, itimes, etimes, treatment_delay, infection_tracker, i)
        
        # Vaccinate
        if vax_scenario != 0 and vax_delay <= i:
            vax_stat, vtimes, vax_eff = vaccinate(G, activity_strat, rel_activity, daily_num_FD, daily_num_SD, 
                                                  fd_eff, sd_eff, vax_stat, vtimes, i, vax_delay, vax_scenario)
        
        # Spread the infection
        G, status, contact_count, infection_tracker = spread(G, N, i, p_infect, status, pcontact_main, pcontact_casual, 
                                          pcontact_onetime, isolation_complier, treatment_delay, vax_eff, infection_tracker)  
        num_contacts[i] = contact_count
        
        
        ### Relationship Updates
        G, onetime_relationships, onetime_degseq_new = update_onetimerelationships(G, onetime_relationships, onetime_prob, rel_activity, activity_prob, activity_strat)
        G, want_main, want_casual, ex_list_main, ex_list_casual = update_relationships(G, want_main, want_casual, ex_list_main, ex_list_casual)

        
        ### Update bookkeeping
        number_E.append(np.sum(status[:,1]))
        number_I.append(np.sum(status[:,2]))
        number_R.append(np.sum(status[:,3]))
        
        
    return(number_E, number_I, number_R, infection_tracker)
    
    
