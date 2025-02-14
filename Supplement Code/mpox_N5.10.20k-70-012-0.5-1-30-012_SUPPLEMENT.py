import networkx as nx
import random
import collections
import numpy as np
import pandas as pd
import itertools
import sys
import os
import time
import pickle

import mpox_utils_supplement 

from mpox_utils_supplement import *


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


p_infect = 0.9
steps = 250
intervention_start = 70

behavior_change_perc = 0.5
vax_delay = 30

N = [5000, 10000, 20000]

behavior_change = [0,1,2]
isolation = 1

for d in range(len(behavior_change)):
    for num in range(len(N)):
        n_initial= int(N[num]/1000)
        
        vax_scenario = behavior_change[d]
        
        sim_string = 'N' + str(n_initial) + 'k-' + str(intervention_start) + '-' + str(behavior_change[d]) + '-' + str(behavior_change_perc) + \
              '-' + str(isolation) + '-' + str(vax_delay) + '-' + str(vax_scenario) + 'SUPPLEMENT'


        start_time = time.time()
        #E_out, I_out, R_out, infection_tracker, wait_time_main_list, wait_time_casual_list, num_contacts, main_degseq, onetime_degseq, rel_activity, activity_strat 
        results = simulate(N[num], n_initial, p_infect, steps, intervention_start, behavior_change[d], 
                         isolation, behavior_change_perc, vax_scenario, vax_delay, daily_num_FD, daily_num_SD)
        end_time = time.time()

        runtime = end_time - start_time

        final_results = [results, runtime]

        folder_path = str('output/' + sim_string + '/')

        if not os.path.exists(folder_path):
            os.makedirs(folder_path)

        with open(str('output/' + sim_string + '/mpox_' + sim_string + '_' + str(seed_num) + '_' + str(date) +'.pkl'), 'wb') as file: 
                # A new file will be created 
                pickle.dump(final_results, file) 