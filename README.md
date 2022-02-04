# Equitable access to COVID-19 vaccines makes a life-saving difference to all countries (Nature Human Behaviour)
## Environment
Python 3.7.1

Matlab 9.5.0.1049112 (R2018b) Update 3

pandas 1.0.5

numpy 1.20.3

matplotlib 3.0.2

The FEX packages COLORSPACE and cbrewer2 are required to generate figure 4:
https://www.mathworks.com/matlabcentral/fileexchange/28790-colorspace-transformations
https://www.mathworks.com/matlabcentral/fileexchange/58350-cbrewer2?s_tid=srchtitle

## Instructions
1 Instruction for data correction in model initialization can be found in read_me_data_corrections.txt

2 Processed data for specific simulation start dates are under the specific_date_info_modified folder.

3 You can get all data for generating figures in the manuscript by running final_exp.py. Create folders 'results' and 'figs' 
under the matlab_files folder to save data.

4 We use Matlab to generate figures. You can find all Matlab codes under the matlab_files folder (summary_update.m).

5 It will only take a few minutes to get all results in the main manuscript except the results in figure 4, which will take about 3 hours.

6 Visualization of the global mobility network is in GMN.pdf

7 Sample data for worldwide air-traffic data set is provided in G_air_2020.gexf. Full data are commercially available 
  from the Official Aviation Guide (https://www.oag.com/) and were used under license for the current study. Due to restrictions 
  in the licensing agreement with OAG, these data are not publicly available.

8 All simulation code is also available at https://doi.org/10.5281/zenodo.5810400.

## Cite this article
Ye, Y., Zhang, Q., Wei, X. et al. Equitable access to COVID-19 vaccines makes a life-saving difference to all countries. Nat Hum Behav (2022). https://doi.org/10.1038/s41562-022-01289-8
