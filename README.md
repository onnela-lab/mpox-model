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
* create_figures_relationship_type.ipynb: creates plots showing simulation results by relationship type
* create_figure_cumulative_edge_graph.ipynb: creates plots showing how edges accumulate over time in the graph

### Supplement Code
This is code written to address reviewer responses and primarily creates information contained in the appendix of the manuscript.
* [mpox_70-2-0.5-012-30-2_SUPPLEMENT.py](Supplement%20Code/mpox_70-2-0.5-012-30-2_SUPPLEMENT.py) - this code runs the sensitivity analyses looking at different values of the infection parameters
* [mpox_N5.10.20k-70-012-0.5-1-30-012_SUPPLEMENT.py](Supplement%20Code/mpox_N5.10.20k-70-012-0.5-1-30-012_SUPPLEMENT.py) - this code runs the population size sensitivity analyses for N = 5,000, N = 10,000, and N = 20,000
* [mpox_N40k-70-012-0.5-1-30-012_SUPPLEMENT.py](Supplement%20Code/mpox_N40k-70-012-0.5-1-30-012_SUPPLEMENT.py) - this code runs the population size sensitivity analyses for N = 40,000
* [mpox_N80k-70-012-0.5-1-30-012_SUPPLEMENT.py](Supplement%20Code/mpox_N80k-70-012-0.5-1-30-012_SUPPLEMENT.py) - this code runs the population size sensitivity analyses for N = 80,000
* [mpox_30to110-2-0.25range-1-0to-30-2_BW.py](Supplement%20Code/mpox_30to110-2-0.25range-1-0to-30-2_BW.py) - this code runs sensitivity analyses looking at best/worst case scenarios for the disease natural history parameters
* [mpox_30to110-2-0.25range-1-0to-30-2_BW_lowprob.py](Supplement%20Code/mpox_30to110-2-0.25range-1-0to-30-2_BW_lowprob.py)- this code runs sensitivity analyses looking at best/worst case scenarios for the disease natural history parameters in the scenario with low transmission probability
* [analyze_compress_supplement.py](Supplement%20Code/analyze_compress_supplement.py) - this code helps analyze the results from the population size code and saves datasets used for the figures
* [clustering_simulations.py](Supplement%20Code/clustering_simulations.py) - this code creates the information for the network structure part of the appendix
* [supplement_figures.py](Supplement%20Code/supplement_figures.py) - This code creates the new figures
* [mpox_utils_supplement.py](Supplement%20Code/mpox_utils_supplement.py) - this is a copy of the original utils file that tracks additional information during simulations required for the response to reviewers


## Parameters for Simulation Scenarios
The following table shows the parameter values used in simulation and is Table 1 in the manuscript. References numbers relate to the references in the manuscript.

### Network Model  
| Parameter         | Value            | Description |
|-------------------|------------------|------------|
| **rt_k**         |                  | Relationship type, percent of nodes (21,23) |
|                  | 47.1            | 0 main, 0 casual |
|                  | 16.7            | 0 main, 1 casual |
|                  | 7.4             | 0 main, 2 casual |
|                  | 22.0            | 1 main, 0 casual |
|                  | 4.7             | 1 main, 1 casual |
|                  | 2.1             | 1 main, 2 casual |
| **π_{o,k}**      |                  | Daily probability of one-time partnership formation (21,23) |
|                  | 0                | Stratum 1 |
|                  | 0.001            | Stratum 2 |
|                  | 0.0054           | Stratum 3 |
|                  | 0.0101           | Stratum 4 |
|                  | 0.0315           | Stratum 5 |
|                  | 0.286            | Stratum 6 |
| **n_{o}**       | Geometric(1 - π_{o,k}) | Number of one-time contacts per day (21) |
| **rd_{m,e}**    | Geometric(1/407) | Duration of main partnership, days (21,23) |
| **rd_{c,e}**    | Geometric(1/166) | Duration of casual partnership, days (21,23) |

### Epidemic Model  
| Parameter        | Value          | Description |
|------------------|---------------|------------|
| **β**            | 0.9           | Probability of transmission per sexual contact |
| **π_m**         | 0.22          | Daily probability of sexual contact (main partnership) (21,23) |
| **π_c**         | 0.14          | Daily probability of sexual contact (casual partnership) (21,23) |
| **t_{e,k}**     | Normal(7, 1)  | Time spent in exposed state, days (1,2,13,25-27) |
| **t_{i,k}**     | Normal(27, 3) | Time spent in infectious state, days (1,2,13,25-27) |

### Interventions  
| Parameter       | Value         | Description |
|-----------------|--------------|------------|
| **π_{b,k}**    | 0.5          | Decrease in probability of one-time partnership formation |
| **π_{v_1}**    | 0.14         | Probability of vaccination (one dose) (32) |
| **π_{v_2}**    | 0.227        | Probability of vaccination (two doses) (32) |
| **VE_1**       | 0.358        | Vaccine efficacy (one dose) (34) |
| **VE_2**       | 0.66         | Vaccine efficacy (two doses) (34) |

*Node-specific attributes use the subscript `X_{x,k}`, edge-specific attributes use `X_{x,e}`.*

