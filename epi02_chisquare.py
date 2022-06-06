#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 14:52:33 2021

@author: bpetros
"""

#%%

## Download relevant packages and set up environment
from adjustText import adjust_text
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
    fname='Montserrat/static/Montserrat-Medium.ttf',
    name='Montserrat')
fm.fontManager.ttflist.insert(0, fe) 
mpl.rcParams['font.family'] = fe.name 
del fe

## Load in cmu_meta.spydata

#%% Chi square analysis for each categorical variable

## class year
grade = pd.read_csv('data/in/class_years.csv')
grade.set_index('year', inplace = True)

# drop years without known population sizes
grade.drop('Year E', inplace = True)
grade.drop('Year G', inplace = True)

spring = []
fall = []

# cases per class year, for fall and for spring
for i in grade.index:
    spring.append(sum((meta.Class_Year == i) & (meta.date < datetime.strptime('2021-01-01', '%Y-%m-%d'))))
    fall.append(sum((meta.Class_Year == i) & (meta.date > datetime.strptime('2021-01-01', '%Y-%m-%d'))))

grade['spring_cases'] = spring
grade['fall_cases'] = fall
grade.to_csv('data/out/class_year_cases.csv')

del fall, i, spring

# determine expected number of cases per class and calculate chi-square statistic
spring_exp = np.multiply(np.divide(grade.spring_pop, sum(grade.spring_pop)), sum(grade.spring_cases))
s = stats.chisquare(grade.spring_cases, spring_exp)
p = ((grade.spring_cases - spring_exp)/spring_exp)

fall_exp = np.multiply(np.divide(grade.fall_pop, sum(grade.fall_pop)), sum(grade.fall_cases))
s = stats.chisquare(grade.fall_cases, fall_exp)
p = ((grade.fall_cases - fall_exp)/fall_exp)

sp = pd.DataFrame(index = grade.decoded)
sp['exp'] = spring_exp.tolist()
sp['oe'] = np.divide(np.subtract(grade.spring_cases.to_numpy(),spring_exp.to_numpy()), spring_exp.to_numpy())
sp['sem'] = 'Spring 2021'
f = pd.DataFrame(index = grade.decoded)
f['exp'] = fall_exp.tolist()
f['oe'] = np.divide(np.subtract(grade.fall_cases.to_numpy(),fall_exp.to_numpy()), fall_exp.to_numpy())
f['sem'] = 'Fall 2020'
sp = pd.concat([sp, f])

del f, s, p


# plot observed vs. expected counts
f, (ax1) = plt.subplots(1,1, figsize = (6,8), dpi = 800)
sns.despine()
col = {'Fall 2020':'#000000','Spring 2021':'#737373'}
ax1.scatter(sp.exp, sp.oe, s = 25, c = sp['sem'].map(col))
line = ax1.plot(range(260), 260*[0], linestyle = 'dashed', color = 'black', alpha = 0.5)
markers = [plt.Line2D([0,0],[0,0],color = color, marker='o', linestyle = '') for color in col.values()]
ax1.legend(markers, col.keys(), numpoints = 1, fontsize = 14, loc = 'upper right')
texts = []
for x, y, text in zip(sp.exp, sp.oe, sp.index):
    texts.append(ax1.text(x, y, text, size = 10))
adjust_text(texts, add_objects = line, force_objects=(0.05, 0.05), arrowprops=dict(arrowstyle="-",
              color='black', alpha=0.5), expand_points=(1.2,1.2))
ax1.set_title('Class Year', size = 18)
ax1.set_xlabel('Expected Cases', size = 14)
ax1.set_ylabel('Observed Cases - Expected Cases)/Expected Cases', size = 14)
f.tight_layout()
plt.savefig('data/out/class_years.png')
plt.savefig('data/out/class_years.svg')

del ax1, col, f, fall_exp, line, markers, spring_exp, sp, text, texts, x, y

#%%

# high vs. medium vs. low contact sports
contact_obs = [sum(sports.sports_cases[sports.Contact_Level == 'High']), 
              sum(sports.sports_cases[sports.Contact_Level == 'Moderate']), 
              sum(sports.sports_cases[sports.Contact_Level == 'Low'])]

contact_exp = np.multiply(np.divide([sum(sports.sports_population[sports.Contact_Level == 'High']), 
              sum(sports.sports_population[sports.Contact_Level == 'Moderate']), 
              sum(sports.sports_population[sports.Contact_Level == 'Low'])], sum(sports.sports_population)), sum(sports.sports_cases))
                          
s = stats.chisquare(contact_obs, contact_exp)  
p = ((contact_obs - contact_exp)/contact_exp)

# indoor vs. outdoor sports
loc_obs = [sum(sports.sports_cases[sports['loc'] == 'indoor']),\
           sum(sports.sports_cases[sports['loc'] == 'outdoor'])]

loc_exp = np.multiply(np.divide([sum(sports.sports_population[sports['loc'] == 'indoor']), 
              sum(sports.sports_population[sports['loc'] == 'outdoor'])], sum(sports.sports_population) - \
                  sum(sports.sports_population[sports['loc'] == 'both'])), sum(sports.sports_cases) - \
                                 sum(sports.sports_cases[sports['loc'] == 'both']))
                          
s = stats.chisquare(loc_obs, loc_exp)  
p = ((loc_obs - loc_exp)/loc_exp)

del contact_obs, contact_exp, loc_obs, loc_exp, s, p

#%%      
               
# sports teams
sport_exp = np.multiply(np.divide(sports.sports_population, sum(sports.sports_population)), sum(sports.sports_cases))
s = stats.chisquare(sports.sports_cases, sport_exp)
p = ((sports.sports_cases - sport_exp)/sport_exp)

# plot observed vs. expected counts
f, (ax1) = plt.subplots(1,1, figsize = (7,7), dpi = 300)
sns.despine()
col = {'High':"#000000", 'Moderate':"#737373", 'Low':"#bdbdbd"}
ax1.scatter(sport_exp, p, s = 25, c = sports['Contact_Level'].map(col))
line = ax1.plot(range(33), 33*[0], linestyle = 'dashed', color = 'black', alpha = 0.5)
markers = [plt.Line2D([0,0],[0,0], color = color, marker='o', linestyle = '') for color in col.values()]
ax1.legend(markers, col.keys(), numpoints = 1, loc = 'lower right')
texts = []
for x, y, text in zip(sport_exp, p, sports['short_name'].tolist()):
    texts.append(ax1.text(x, y, text, size = 8))
adjust_text(texts, add_objects = line, force_objects=(3, 3), arrowprops=dict(arrowstyle="-",
              color='black', alpha=0.5), expand_points=(1.2,1.2))
ax1.set_title('Sports Teams', size = 18)
ax1.set_xlabel('Expected Cases', size = 14)
ax1.set_ylabel('(Observed Cases - Expected Cases)/Expected Cases', size = 14)
plt.savefig('data/out/sports.png')
plt.savefig('data/out/sports.svg')

del ax1, col, f, line, markers, p, s, sport_exp, text, texts, x, y

#%%

# residence halls
hall_exp = np.multiply(np.divide(halls.hall_population, sum(halls.hall_population)), sum(halls.hall_cases))
s = stats.chisquare(halls.hall_cases, hall_exp)
p = ((halls.hall_cases - hall_exp)/hall_exp)/s[0]

# plot observed vs. expected counts
f, (ax1) = plt.subplots(1,1, figsize = (6, 7), dpi = 300)
sns.despine()
ax1.scatter(hall_exp, (halls.hall_cases-hall_exp)/hall_exp, s = 50, color = 'black')
line = ax1.plot(range(15, 60, 1), 45*[0], linestyle = 'dashed', color = 'black', alpha = 0.5)
texts = []
for x, y, text in zip(hall_exp+1, (halls.hall_cases-hall_exp)/hall_exp, halls.index):
    texts.append(ax1.text(x, y, text, size = 10))
adjust_text(texts, add_objects = line, force_objects=(10, 10), arrowprops=dict(arrowstyle="-",
              color='black', alpha=0.1), expand_points=(1.2,1.2))
ax1.set_title('Residence Halls', size = 18)
ax1.set_xlabel('Expected Cases', size = 14)
ax1.set_ylabel('(Observed Cases - Expected Cases)/Expected Cases', size = 14)
plt.savefig('data/out/halls.png')
plt.savefig('data/out/halls.svg')

del ax1, f, hall_exp, line, p, s, text, texts, x, y

## save cmu_meta_plus_grades.spydata