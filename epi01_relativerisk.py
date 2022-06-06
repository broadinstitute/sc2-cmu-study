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
# calc sports team cases
cases = []
for team in sports.Sports_team:
    cases.append(sum(meta.Sports_team == team) + sum(meta.Additional_Sports_team == team))
sports['sports_cases'] = cases
del cases, team
sports.set_index('Sports_team', inplace = True)
sports.to_csv('data/out/sports_cases.csv')

# residence halls
halls = pd.read_csv('data/in/hall_populations.csv')
# calc res hall cases
cases = []
for hall in halls.hall:
    cases.append(sum(meta.Residence_hall == hall))
halls['hall_cases'] = cases
del cases, hall
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
male = round(cmu_pop * 0.458) #from https://www.coloradomesa.edu/institutional-research/documents/student-profiles/eot_fall_total.pdf
campus = sum(halls.hall_population)

# sports team vs. not
s = rr(cmu_pop, sum(sports.sports_population), len(meta), sum(meta.Sports_team.notna()))

# male vs. female
m = rr(cmu_pop, male, len(meta), sum(meta.Sex == 'M'))

# on- vs. off- campus
c = rr(cmu_pop, campus, len(meta), sum(meta.Lives_on_campus == 'Yes'))

#%% 

# plot relative risks
plt.rcParams['figure.dpi'] = 800
plt.rcParams['savefig.dpi'] = 800
 
l = ['Athletes vs. Non-Atheletes', 'Males vs. Females', 'On- vs. Off-Campus']
r = [s[0], m[0], c[0]]
err = [[s[1], m[1], c[1]], [s[2], m[2], c[2]]]

sns.set_style("whitegrid", {'axes.grid' : False})
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
ax.set_xlabel("Student Groups", fontsize = 14)
ax.set_ylabel("Relative Risk", fontsize = 14)
ax.set(xticks=[0, 1, 2])
ax.set_xticklabels(l, rotation=45, fontsize = 12)
fig.tight_layout()
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

