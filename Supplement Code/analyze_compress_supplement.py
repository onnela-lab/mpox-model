import networkx as nx
import random
import collections
import numpy as np
import pandas as pd
import itertools
import sys
import pickle
import os
import zipfile


import mpox_utils_supplement 

from mpox_utils_supplement import *




# Get other arguments
N = [5000, 10000, 20000, 40000, 80000]

p_infect = 0.9
steps = 250
intervention_start = 70
behavior_change_perc = 0.5
vax_delay = 30

reps = 50

date = '2025-02-05'


behavior_change = [0,1,2]
isolation = 1

for d in range(len(behavior_change)):
    for num in range(len(N)):
        vax_scenario = behavior_change[d]
        n_initial= int(N[num]/1000)
        
        sim_string = 'N' + str(n_initial) + 'k-' + str(intervention_start) + '-' + str(behavior_change[d]) + '-' + str(behavior_change_perc) + '-' + str(isolation) + '-' + str(vax_delay) + '-' + str(vax_scenario) + 'SUPPLEMENT'
        
        total_infection = np.zeros((reps, 250))
        #infection_tracker = []
        wait_time_main_list = []
        wait_time_casual_list = []
        num_contacts = []
        main_degseq = np.zeros((N[num], reps))
        onetime_degseq  = np.zeros((N[num], reps))
        rel_activity  = np.zeros((N[num], reps))
        activity_strat = np.zeros((N[num], reps))
        runtime = []
        

        for i in range(1, reps+1):
        
            j = i-1
            with open(str('/n/home04/ecrenshaw/mpox/output/' + str(sim_string) + '/mpox_' + str(sim_string) + '_' + str(i) + '_' + str(date) + '.pkl'), 'rb') as file: 
                # A new file will be created 
                res = pickle.load(file) 

            (E_out1, I_out1, R_out1, infection_tracker1, wait_time_main_list1, 
             wait_time_casual_list1, num_contacts1, main_degseq1, onetime_degseq1, rel_activity1, activity_strat1) = res[0]
            runtime1 = res[1]
            
            rep_ti = [x+y+z for x,y,z in zip(E_out1, I_out1, R_out1)]
            ti = rep_ti + [rep_ti[-1]]*(steps - len(rep_ti))

            runtime.append(runtime1)              
            total_infection[j, :] = ti
            #infection_tracker.append(infection_tracker1)
            wait_time_main_list.append(wait_time_main_list1)
            wait_time_casual_list.append(wait_time_casual_list1)
            num_contacts.append(num_contacts1)
            main_degseq[:,j] = main_degseq1
            onetime_degseq[:,j] = onetime_degseq1
            rel_activity[:,j] = rel_activity1
            activity_strat[:,j] = activity_strat1
              
                
        # analyze the data so I don't have to save all of it
            
        ## wait times
        wait_time_main_avg = [np.mean(x) for x in wait_time_main_list]
        wait_time_casual_avg = [np.mean(x) for x in wait_time_casual_list]
            
        wait_time_main_perc = [sum(1 for value in row if value > 1) / len(row) if len(row) > 0 else 0 for row in wait_time_main_list]
        wait_time_casual_perc = [sum(1 for value in row if value > 1) / len(row) if len(row) > 0 else 0 for row in wait_time_casual_list]
            
        ## number of contacts
        rel_activity_sum = np.zeros((reps,6))
        activity_strat_sum = np.zeros((reps,6))

        rel_activity_num = np.zeros((reps,6))
        activity_strat_num = np.zeros((reps,6))
        rel_act_num = np.zeros((reps,6,6))
        rel_act_sum = np.zeros((reps,6,6))
        
        avg_contacts = np.sum(num_contacts, axis = 2)

        for i in range(reps):
            for j in range(N[num]):
                rel_activity_sum[i, int(rel_activity[j,i])] += avg_contacts[i,j]
                rel_activity_num[i, int(rel_activity[j,i])] += 1
                activity_strat_sum[i, int(activity_strat[j,i])] += avg_contacts[i,j]
                activity_strat_num[i, int(activity_strat[j,i])] += 1
                rel_act_num[i, int(rel_activity[j,i]), int(activity_strat[j,i])] += 1
                rel_act_sum[i, int(rel_activity[j,i]), int(activity_strat[j,i])] += avg_contacts[i,j]

        #don't divide by 0
        rel_activity_avg = np.divide(rel_activity_sum, rel_activity_num, where=rel_activity_num>0)
        activity_strat_avg = np.divide(activity_strat_sum, activity_strat_num, where=activity_strat_num>0)
        rel_act_avg = np.divide(rel_act_sum, rel_act_num, where=rel_act_num>0)

            
        save_data = [wait_time_main_avg, wait_time_casual_avg, wait_time_main_perc, wait_time_casual_perc, rel_activity_avg, activity_strat_avg, rel_act_avg, runtime, total_infection]
            
        # save new data
        with open(str('output/' + sim_string + '/mpox_' + sim_string + '_analyze_' + str(date) +'.pkl'), 'wb') as file: 
                # A new file will be created 
                pickle.dump(save_data, file) 
            
            
        
            