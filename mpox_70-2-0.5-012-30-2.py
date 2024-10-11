import networkx as nx
import random
import collections
import numpy as np
import pandas as pd
import itertools
import sys

import os
import mpox_utils 

from mpox_utils import *


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

# Get other arguments
N = 10000
n_initial= 10
p_infect = 0.9
steps = 250
intervention_start = 70
behavior_change = 2
behavior_change_perc = 0.5
vax_scenario = 2
vax_delay = 30

sim_string = str(intervention_start) + '-' + str(behavior_change) + '-' + str(behavior_change_perc) + \
          '-' + '012' + '-' + str(vax_delay) + '-' + str(vax_scenario)

rstar_list = [7,14,21,28,35,42,49,56,63,70,77,84,91,98,105]
isolation = [0,1,2]

final_infection = np.zeros((len(isolation),steps+1))
all_repro = np.zeros((len(isolation),len(rstar_list)+3))
all_tracker = np.zeros((N*3,4))

for d in range(len(isolation)):
    
    E_out, I_out, R_out, infection_tracker = simulate(N, n_initial, p_infect, steps, intervention_start, behavior_change, 
                 isolation[d], behavior_change_perc, vax_scenario, vax_delay, daily_num_FD, daily_num_SD)

    total_infection = [x+y+z for x,y,z in zip(E_out, I_out, R_out)]
    ti = total_infection + [total_infection[-1]]*(steps - len(total_infection))
    
    # save infection tracker
    start = d*N
    stop = start + N
    all_tracker[start:stop,:] = infection_tracker
    
    # Final infected
    row_num = d
    final_infection[row_num,:] = [isolation[d]] + ti

    repro_num = np.zeros(len(rstar_list)+1)
    # R0
    initial_infections = np.where(infection_tracker[:,0] == -1)[0].tolist()
    secondary_infections = np.where(np.isin(infection_tracker[:,0],initial_infections))[0].tolist()
    repro_num[0] = len(secondary_infections)/len(initial_infections)

    #Rstar
    for star in range(len(rstar_list)):
        infectious = np.where((infection_tracker[:,1] < rstar_list[star]) & (infection_tracker[:,2] > (rstar_list[star]-7)))[0].tolist()
        infected = np.where(np.isin(infection_tracker[:,0], infectious))[0].tolist()
        print("Rstar date: ", rstar_list[star])
        print("Number people infectious: ", len(infectious))
        print("Number of people infected by them: ", len(infected))
        if len(infectious) > 0:
            repro_num[star+1] = len(infected)/len(infectious)
        else: 
            repro_num[star+1] = 0

    # infections caused by one-time partners
    ot_infections = len(np.where(infection_tracker[:,3]==3)[0].tolist())
    
    all_repro[row_num,:] = [isolation[d]] + list(repro_num) + [ot_infections/final_infection[row_num,-1]]

# Write to file
df = pd.DataFrame(final_infection)
df2 = pd.DataFrame(all_repro)
df3 = pd.DataFrame(all_tracker)

folder_path = str('output/' + sim_string + '/')

if not os.path.exists(folder_path):
    os.makedirs(folder_path)


df.to_csv('output/' + sim_string + '/mpox_' + sim_string + '_' + str(seed_num) + '_' + str(date)+'.csv', index=False)
df2.to_csv('output/' + sim_string + '/rstar_' + sim_string + '_' + str(seed_num) + '_' + str(date)+'.csv', index=False)
df3.to_csv('output/' + sim_string + '/tracker_' + sim_string + '_' + str(seed_num) + '_' + str(date)+'.csv', index=False)