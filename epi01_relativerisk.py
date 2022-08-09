#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 10:03:09 2021

@author: bpetros
"""
#%%
## Download relevant packages and set up environment
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scipy.stats as stats
import seaborn as sns
os.chdir('/Users/bpetros/OneDrive - Harvard University/cmu')

## Set default font to Montserrat
import matplotlib as mpl
import matplotlib.font_manager as fm

fe = fm.FontEntry(
    fname='/Users/bpetros/opt/anaconda3/pkgs/matplotlib-base-3.5.1-py38hfb0c5b7_0/lib/python3.8/site-packages/matplotlib/mpl-data/fonts/ttf/Montserrat-Medium.ttf',
    name='Montserrat')
fm.fontManager.ttflist.insert(0, fe) 
mpl.rcParams['font.family'] = fe.name 
del fe, fm

#%%

## Load in relevant data files
meta = pd.read_csv('data/in/cmu_metadata.csv')
meta.set_index('barcode', inplace = True)
meta['date'] = pd.to_datetime(meta.test_day)
# Remove dates pre-fall semester
meta = meta.loc[meta.date >= datetime.strptime('2020-08-20', '%Y-%m-%d'), :] #removed 25 cases
# Remove dates post-spring semester
meta = meta.loc[meta.date <= datetime.strptime('2021-04-30', '%Y-%m-%d'), :] #removed 1 case
# Remove dates between semesters
meta = meta.loc[(meta.date <= datetime.strptime('2020-11-20', '%Y-%m-%d')) | 
                (meta.date >= datetime.strptime('2021-01-18', '%Y-%m-%d')), :] #removed 188 cases

# sports teams
sports = pd.read_csv('data/in/sports_populations.csv')
sports.set_index('Sports_team', inplace = True)
sport_test = pd.read_csv('data/in/sport-tests.csv', index_col = 0)
# calc sports team cases and tests
cases = []
tests = []
for team in sports.index:
    cases.append(sum(meta.Sports_team == team) + sum(meta.Additional_Sports_team == team))
    try:
        tests.append(np.nansum(sport_test[sports.loc[team].team_name]))
    except:
        tests.append(np.nan)
sports['sports_cases'] = cases
sports['sports_tests'] = tests
del cases, sport_test, team, tests
sports.to_csv('data/out/sports_cases.csv')

# no association b/w sports team testing rates and sports team incidence rates
nansport = sports[~np.isnan(sports.sports_tests)]
corr, pv = stats.pearsonr(nansport.sports_tests/nansport.sports_population, nansport.sports_cases/nansport.sports_population)
del corr, nansport, pv

# residence halls
halls = pd.read_csv('data/in/hall_populations.csv')
hall_test = pd.read_csv('data/in/hall-tests.csv', index_col = 0)
# calc res hall cases
cases = []
tests = []
for hall in halls.hall:
    cases.append(sum(meta.Residence_hall == hall))
    tests.append(np.nansum(hall_test[hall]))
halls['hall_cases'] = cases
halls['hall_tests'] = tests
del cases, hall, hall_test, tests
halls.set_index('hall', inplace = True)
halls.to_csv('data/out/hall_cases.csv')

#%% Relative risk analyses

## relative risk function
def rr(pop, sub_pop, tot_case, sub_case):
    # pop = n of entire population
    # sub_pop = n of members of pop belonging to sub-group
    # tot = n of cases in entire population
    # sub_case = n of cases in sub-group
    
    rel_risk = (sub_case/sub_pop)/((tot_case - sub_case)/(pop - sub_pop))
    
    #log(rel_risk) approx. normally distributed --> calculate CI in log space
    se = np.sqrt(((sub_pop - sub_case)/sub_case)/sub_pop + ((pop - sub_pop - (tot_case - sub_case))/(tot_case - sub_case))/(pop - sub_pop))
    ci_log = [np.log(rel_risk) - (1.96 * se), np.log(rel_risk) + (1.96 * se)]
    ci = np.exp(ci_log)
    return rel_risk, ci[0], ci[1]

#%% 

# assumed values
cmu_pop = (8350 + 7670)/2 #from CMU, averaged fall and spring
cmu_tests = 10842+12891 #from Fathom, fall and spring
male = round(cmu_pop * 0.458) #from https://www.coloradomesa.edu/institutional-research/documents/student-profiles/eot_fall_total.pdf
campus = sum(halls.hall_population)

# sports team vs. not
s = rr(cmu_pop, sum(sports.sports_population), len(meta), sum(meta.Sports_team.notna()))

# sports testing ratio
spp = np.nansum(sports.sports_tests)/sum(sports.sports_population)
nspp = (cmu_tests - np.nansum(sports.sports_tests))/(cmu_pop - sum(sports.sports_population))
spp/nspp

# male vs. female
m = rr(cmu_pop, male, len(meta), sum(meta.Sex == 'M'))

# on- vs. off- campus
c = rr(cmu_pop, campus, len(meta), sum(meta.Lives_on_campus == 'Yes'))

# hall testing ratio
hpp = np.nansum(halls.hall_tests)/sum(halls.hall_population)
nhpp = (cmu_tests - np.nansum(halls.hall_tests))/(cmu_pop - sum(halls.hall_population))
hpp/nhpp

del spp, nspp, hpp, nhpp

#%% 

# plot relative risks
plt.rcParams['figure.dpi'] = 800
plt.rcParams['savefig.dpi'] = 800
 
l = ['Athletes vs. Non-Atheletes', 'Males vs. Females', 'On- vs. Off-Campus']
r = [s[0], m[0], c[0]]
err = [[s[1], m[1], c[1]], [s[2], m[2], c[2]]]

sns.set_style("whitegrid", {'axes.grid' : False})
sns.set(font="Montserrat")
fig, ax = plt.subplots(figsize = (4, 5))
ax.axhline(1, ls='--', linewidth=1, color='black')
sns.despine()

n = 0
for i in range(len(r)):
    x = [n,n]
    y = [err[0][i], err[1][i]]
    ax.plot(x, y, "_-", markersize = 15, markeredgewidth= 3, linewidth = 3, color="black")
    x = [n]
    y = [r[i]]
    ax.plot(x, y, "o", color="black", markersize = 10)
    n += 1
ax.set_ylim([0.75, 2.9])
ax.set_ylabel("Relative Risk", fontsize = 14, font = 'Montserrat')
ax.set(xticks=[0, 1, 2])
ax.set_xticklabels(l, rotation=45, fontsize = 12, font = 'Montserrat')
plt.savefig('data/out/relative_risk.png')
plt.savefig('data/out/relative_risk.svg')
del ax, c, err, fig, i, l, m, n, r, s, x, y


# save data as cmu_meta.spydata


#%% Determine proportion of students with given risk factor for Fig 1A

# values from CMU
fall_pop = 8350
spring_pop = 7670
sport_pop = sum(sports.sports_population)

# test positive, by semester
fall_meta = meta.loc[(meta.date < datetime.strptime('2021-01-01', '%Y-%m-%d')) & meta.Class_Year.notna(), :]
spring_meta = meta.loc[(meta.date > datetime.strptime('2021-01-01', '%Y-%m-%d')) & meta.Class_Year.notna(), :]


# fall demographics
fall_ug = sum(fall_meta.Class_Year.isin(['Year A', 'Year B', 'Year C', 'Year D']))
fall_campus = sum(fall_meta.Lives_on_campus == 'Yes')
fall_sports = sum(fall_meta.Sports_team.notna())
fall_male = sum(fall_meta.Sex == 'M')


# spring demographics
spring_ug = sum(spring_meta.Class_Year.isin(['Year A', 'Year B', 'Year C', 'Year D']))
spring_campus = sum(spring_meta.Lives_on_campus == 'Yes')
spring_sports = sum(spring_meta.Sports_team.notna())
spring_male = sum(spring_meta.Sex == 'M')

