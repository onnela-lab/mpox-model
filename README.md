# mpox-model

This repository contains Python code for the study "Modeling the 2022 Mpox Outbreak with a Mechanistic Network Model".

For any further inquiries, please email emma_crenshaw@g.harvard.edu.

## File Overview

### Naming Conventions
All the files follow the same file naming convention following "mpox_": behavior change start time, behavior change scenario, reduction in the probability of a one-time partnership (0.25 indicates a 75% reduction in one-time partnership, i.e., $\pi_{0,k} \times 0.25$), isolation scenario, vaccination start time, and vaccination scenario. The scenarios are as follows:
* Behavior change:
  * 0 = No behavior change
  * 1 = Universal behavior change (all individuals participate)
  * 2 = Targeted behavior change (only individuals in the two highest strata of sexual activity participate)

* Isolation scenario:
  * 1 = Full compliance with isolation
  * 2 = Partial compliance with isolation

* Vaccination scenario:
  * 0 = No vaccination
  * 1 = Universal vaccination availability (all individuals can recieve vaccination)
  * 2 = Targeted vaccination (only individuals in the two highest strata of sexual activity can be vaccinated)

Therefore, the file "mpox_30to110-2-0.5-1-0to-30-2.py" runs code for simulations with behavior change that begins between day 30 and day 110, behavior change scenario 2, 50% reduction in the probability of a one-time partner, isolation scenario 1, vaccination beginning on day 30, and vaccination scenario 2.
  
### Core Functions and Input Data
* mpox_utils.py: all methods and helper functions to run the simulations
* mpox_vax_coverage_data.xlsx: file containing the number of vaccines available during each week of the simulation (see manuscript for further details)

### Generate Results
* Simulations with a range of intervention timings (behavior change scenario 2, isolation scenario 2 vaccination scenario 2)
  * [mpox_30to110-2-0.25-2-0to-30-2.py](mpox_30to110-2-0.25-2-0to-30-2.py)
  * [mpox_30to110-2-0.25-2-5to30-2.py](mpox_30to110-2-0.25-2-5to30-2.py)
  * [mpox_30to110-2-0.5-2-0to-30-2.py](mpox_30to110-2-0.5-2-0to-30-2.py)
  * [mpox_30to110-2-0.5-2-5to30-2.py](mpox_30to110-2-0.5-2-5to30-2.py)
  * [mpox_30to110-2-0.75-2-0to-30-2.py](mpox_30to110-2-0.75-2-0to-30-2.py)
  * [mpox_30to110-2-0.75-2-5to30-2.py](mpox_30to110-2-0.75-2-5to30-2.py)

* Simulations with a range of intervention timings (behavior change scenario 2, isolation scenario 1 vaccination scenario 2)
  * [mpox_30to110-2-0.25-1-0to-30-2.py](mpox_30to110-2-0.25-1-0to-30-2.py)
  * [mpox_30to110-2-0.25-1-5to30-2.py](mpox_30to110-2-0.25-1-5to30-2.py)
  * [mpox_30to110-2-0.5-1-0to-30-2.py](mpox_30to110-2-0.5-1-0to-30-2.py)
  * [mpox_30to110-2-0.5-1-5to30-2.py](mpox_30to110-2-0.5-1-5to30-2.py)
  * [mpox_30to110-2-0.75-1-0to-30-2.py](mpox_30to110-2-0.75-1-0to-30-2.py)
  * [mpox_30to110-2-0.75-1-5to30-2.py](mpox_30to110-2-0.75-1-5to30-2.py)

* Simulations with a range of interventions and isolation scenarios
  * [mpox_70-0-0.5-012-30-0.py](mpox_70-0-0.5-012-30-0.py)
  * [mpox_70-1-0.5-012-30-0.py](mpox_70-1-0.5-012-30-0.py)
  * [mpox_70-1-0.5-012-30-1.py](mpox_70-1-0.5-012-30-1.py)
  * [mpox_70-2-0.5-012-30-0.py](mpox_70-2-0.5-012-30-0.py)
  * [mpox_70-2-0.5-012-30-2.py](mpox_70-2-0.5-012-30-2.py)

### Data Processing and Visualization
* concatenate_simulations.ipynb: concatenates output from embarassingly parallelized computing cluster output into one file
* create_figures_main: creates plots comparing interventions and intervention timings
* create_figures_relationship_type.ipynb: creates plots showing simulation results by relationship tpe
* create_cumulative_edge_graph.ipynb: creates plots showing how edges accumulate over time in the graph
