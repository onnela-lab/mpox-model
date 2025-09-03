import networkx as nx
import random
import collections
import numpy as np
import pandas as pd
import itertools
import sys
import os

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
vax_inc = 1
N = 10000
n_initial= 10
p_infect = 0.5
steps = 250
intervention_start = list(range(30,120,10))
behavior_change = 0
isolation_list = [1]
behavior_change_perc = [0.5]
vax_scenario = 0
vax_delay = [5,10,15,20,25,30]

# Best case, worst case (Miura et al 2022, Madewell et al 2022)mu_e = 7.6, sigma_e = 4.9, mu_i = 27, sigma_i = 3
mu_e_list = [9.9, 5.6]
sigma_e_list = [1.96, 0.74]
mu_i_list = [14, 28]
sigma_i_list = [3, 3]
scenario_list = ['best', 'worst']
scenario_int = [1,0]

rstar_list = [7,14,21,28,35,42,49,56,63,70,77,84,91,98,105]

final_infection = []

sim_string = '30to110' + '-' + str(behavior_change) + '-' + 'bcpRange' + \
                '-' + '12' + '-' + '5to30' + '-' + str(vax_scenario) + '_BW' +'_lowprob_' + 'SUPPLEMENT'

for a in range(len(mu_e_list)):
    mu_e = mu_e_list[a]
    sigma_e = sigma_e_list[a]
    mu_i = mu_i_list[a]
    sigma_i = sigma_i_list[a]
    scenario = scenario_list[a]
    
    print("Scenario: ", scenario)
    print("Mu e: ", mu_e, ", sigma e: ", sigma_e)
    print("Mu i: ", mu_i, ", sigma i: ", sigma_i)

    

    for c in range(len(behavior_change_perc)):
        for d in range(len(vax_delay)):
            for s in range(len(intervention_start)):
                for b in range(len(isolation_list)):
                    

                    E_out, I_out, R_out, infection_tracker, wait_time_main_list, wait_time_casual_list, num_contacts, main_degseq, onetime_degseq, rel_activity, activity_strat = simulate(N, n_initial, p_infect, steps, intervention_start[s], behavior_change, isolation_list[b], behavior_change_perc[c], vax_scenario, vax_delay[d], daily_num_FD, daily_num_SD, vax_inc, mu_e, sigma_e, mu_i, sigma_i)

                    total_infection = [x+y+z for x,y,z in zip(E_out, I_out, R_out)]
                    ti = total_infection + [total_infection[-1]]*(steps - len(total_infection))
                    final_infection_all = [intervention_start[s], vax_delay[d], isolation_list[b], scenario_int[a], behavior_change_perc[c]] + ti

                    # Final infected
                    if len(final_infection)<1:
                        final_infection = final_infection_all.copy()
                    else:
                        final_infection = np.vstack((final_infection, final_infection_all))
                    


folder_path = str('output/' + sim_string + '/')

if not os.path.exists(folder_path):
    os.makedirs(folder_path)

# Write to file
df = pd.DataFrame(final_infection)
#df2 = pd.DataFrame(all_repro)


df.to_csv('output/' + sim_string + '/mpox_' + sim_string + '_' + str(seed_num) + '_' + str(date)+'.csv', index=False)
#df2.to_csv('output/' + sim_string + '/rstar_' + sim_string + '_' + str(seed_num) + '_' + str(date)+'.csv', index=False)