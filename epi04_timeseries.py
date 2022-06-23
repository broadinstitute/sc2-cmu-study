#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 14:09:18 2021

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
    fname='Montserrat/static/Montserrat-Medium.ttf',
    name='Montserrat')
fm.fontManager.ttflist.insert(0, fe) 
mpl.rcParams['font.family'] = fe.name 
del fe

## Load in cmu_meta.spydata

#%%

# Assumed values
col_pop = 5.759 * 10**6 #2019, US Census Bureau
mesa_pop = 154210 #2019, US Census Bureau

# Load in CMU-specific test data
tests = pd.read_csv('data/in/testing.csv')
tests['dt'] = pd.to_datetime(tests.date)
# Remove dates between semesters
tests = tests.loc[(tests.dt <= datetime.strptime('2020-11-20', '%Y-%m-%d')) | (tests.dt >= datetime.strptime('2021-01-18', '%Y-%m-%d')), :] 

# Load in county-level case data
mesa = pd.read_csv('data/in/us-counties.csv')
mesa = mesa[mesa.county == 'Mesa']
mesa = mesa[mesa.state == 'Colorado']

# Load in state-level case data
colorado = pd.read_csv('data/in/us-counties.csv')
colorado = colorado[colorado.state == 'Colorado']

# Create df of daily COVID incidence in Mesa vs. at CMU
incidence = pd.DataFrame()
incidence['date'] = mesa.date
incidence.set_index('date', inplace = True)
incidence.sort_values('date', axis = 0, inplace = True)

#%%

# Add Mesa county case, death incidence ratess to df
mesa_case = []
mesa_death = []
for i in range(len(incidence.index)):
    if i == 0:
        mesa_case.extend(mesa.cases[mesa.date == incidence.index[i]].tolist())
        mesa_death.extend(mesa.deaths[mesa.date == incidence.index[i]].tolist())
    else:
        mesa_case.extend(np.subtract(mesa.cases[mesa.date == incidence.index[i]].tolist(), mesa.cases[mesa.date == incidence.index[i-1]].tolist()))
        mesa_death.extend(np.subtract(mesa.deaths[mesa.date == incidence.index[i]].tolist(), mesa.deaths[mesa.date == incidence.index[i-1]].tolist()))
incidence['mesa_case'] = mesa_case
incidence['mesa_death'] = mesa_death
del i, mesa, mesa_case, mesa_death


# Add Colorado case, death incidence rates to df
col_case = []
col_death = []
for i in range(len(incidence.index)):
    if i == 0:
        col_case.extend([sum(colorado.cases[colorado.date == incidence.index[i]].tolist())])
        col_death.extend([sum(colorado.deaths[colorado.date == incidence.index[i]].tolist())])
    else:
        col_case.extend([np.subtract(sum(colorado.cases[colorado.date == incidence.index[i]].tolist()), sum(colorado.cases[colorado.date == incidence.index[i-1]].tolist()))])
        col_death.extend([np.subtract(sum(colorado.deaths[colorado.date == incidence.index[i]].tolist()), sum(colorado.deaths[colorado.date == incidence.index[i-1]].tolist()))])
incidence['col_case'] = col_case
incidence['col_death'] = col_death
del i, colorado, col_case, col_death


# Add CMU incidence and tests to df
cmu_case = []
cmu_test = []
for i in range(len(incidence.index)):
    cmu_case.append(sum(meta['date'] == incidence.index[i]))
    t = tests.tests[tests.date == incidence.index[i]].tolist()
    if t == []:
        cmu_test.extend([0])
    else: 
        cmu_test.extend(t)
incidence['cmu_case'] = cmu_case  
incidence['cmu_test'] = cmu_test
incidence['prct_pos'] = np.divide(cmu_case, cmu_test)
del i, cmu_case, cmu_test, t

#%%

# Convert to datetime objects and filter date range
dt = []
for i in range(len(incidence)):
    dt.append(datetime.strptime(incidence.index[i], '%Y-%m-%d'))
incidence['dt'] = dt
del dt, i

# Filter by semester date ranges (Aug 20 - Apr 30)
incidence = incidence[incidence.dt > datetime.strptime('2020-08-19', '%Y-%m-%d')]
incidence = incidence[incidence.dt < datetime.strptime('2021-05-01', '%Y-%m-%d')]

# Set values outside of semester range equal to nan
incidence.loc[(incidence.dt < datetime.strptime('2021-01-18', '%Y-%m-%d')) &
              (incidence.dt > datetime.strptime('2020-11-20', '%Y-%m-%d')), ['cmu_case']] = np.nan
incidence.loc[(incidence.dt < datetime.strptime('2021-01-18', '%Y-%m-%d')) &
              (incidence.dt > datetime.strptime('2020-11-20', '%Y-%m-%d')), ['cmu_test']] = np.nan
incidence.loc[(incidence.dt < datetime.strptime('2021-01-18', '%Y-%m-%d')) &
              (incidence.dt > datetime.strptime('2020-11-20', '%Y-%m-%d')), ['prct_pos']] = np.nan

# Calc semesterly positivity rates
fall = incidence.loc[(incidence.dt < datetime.strptime('2021-01-01', '%Y-%m-%d')),]
spring = incidence.loc[(incidence.dt >= datetime.strptime('2021-01-01', '%Y-%m-%d')),]
fall_pos = np.nansum(fall.cmu_case)/np.nansum(fall.cmu_test)
spring_pos = np.nansum(spring.cmu_case)/np.nansum(spring.cmu_test)
del fall, fall_pos, spring, spring_pos

#%%

# Find max correlation b/w CMU cases and Mesa case counts, death rate

corr_data = incidence[~incidence['cmu_case']. isna()]

def crosscor(x, y, max_lag):
    '''
    finds maximum correlation between x and y given lag in [1, max_lag]

    Parameters
    ----------
    x : np array or pd series
    y : np array or pd series
    max_lag : int

    Returns
    -------
    maximum correlation and corresponding lag time
    '''
    
    corr = []
    for i in range(1, max_lag + 1):
        corr.append(np.corrcoef(x[:-i], y[i:])[0][1])     
    return(corr, max(corr), np.argmax(corr)+1)
          
c_corr, c_case, l_case = crosscor(corr_data.mesa_case, corr_data.cmu_case, 14)
m_corr, m_case, ml_case = crosscor(corr_data.cmu_case, corr_data.mesa_case, 14)
d_corr, c_death, l_death = crosscor(corr_data.cmu_case, corr_data.mesa_death, 14) 

del c_case, c_death, d_corr, l_case, l_death, m_case, ml_case


#%%

plt.rcParams['figure.dpi'] = 800
plt.rcParams['savefig.dpi'] = 800

# Plot correlation as function of lag time
fig, (ax1) = plt.subplots(1, 1, figsize = (5, 3))
sns.despine()
c_corr.append(np.corrcoef(corr_data.cmu_case, corr_data.mesa_case)[0][1])
c_corr.extend(m_corr)
ax1.vlines(3, 0, 0.75, linestyles='dashed', color = 'black')
ax1.plot(list(range(-14, 15, 1)), c_corr)
ax1.set_xlabel('Lag Time (Days)')
ax1.set_ylabel('Correlation')
plt.savefig('data/out/rolling.png', dpi = 800)
plt.savefig('data/out/rolling.svg')

del ax1, c_corr, fig, m_corr

#%%

# Plot weekly counts and tests for CMU (per 100K)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (7, 8))
sns.despine()
ax1.semilogy(incidence.index, 100000 * incidence.cmu_case.rolling(7).sum()/cmu_pop) #semilog
ax1.semilogy(incidence.index, 100000 * incidence.cmu_test.rolling(7).sum()/cmu_pop) #semilog
ax1.legend(labels = ['Cases', 'Tests'], loc = 'lower right')
ax1.xaxis.set_major_locator(plt.MaxNLocator(10))
ax1.set_ylabel('Weekly COVID-19 Incidence per 100,000')
ax1.set_ylim(1, 10**5)
ax1.tick_params(labelbottom=False)

# Plot daily tests and % positivity for CMU
ax2.plot(incidence.index, incidence.prct_pos)
ax2.xaxis.set_major_locator(plt.MaxNLocator(10))
for label in ax2.get_xticklabels():
        label.set_rotation(90)
ax2.set_xlabel('Date')
ax2.set_ylabel('Percent of Tests Positive')
fig.tight_layout()
plt.savefig('data/out/test_pos.png', dpi = 800)
plt.savefig('data/out/test_pos.svg')

del ax1, ax2, fig, label

#%%

# Plot weekly counts for CMU, Mesa, and Colorado (per 100K)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (6, 8))
sns.despine()
ax1.plot(incidence.index, 100000 * incidence.mesa_case.rolling(7).sum()/mesa_pop)
ax1.plot(incidence.index, 100000 * incidence.col_case.rolling(7).sum()/col_pop)
ax1.plot(incidence.index, 100000 * incidence.cmu_case.rolling(7).sum()/cmu_pop)
ax1.legend(labels = ['Mesa County', 'Colorado State', 'CMU'])
ax1.xaxis.set_major_locator(plt.MaxNLocator(10))
ax1.set_ylabel('Weekly COVID-19 Incidence per 100,000')
ax1.tick_params(labelbottom=False)

# Plot weekly deaths for Mesa, Colorado (per 100K)
ax2.plot(incidence.index, 100000 * incidence.mesa_death.rolling(7).sum()/mesa_pop)
ax2.plot(incidence.index, 100000 * incidence.col_death.rolling(7).sum()/col_pop)
ax2.legend(labels = ['Mesa County', 'Colorado State'])
ax2.xaxis.set_major_locator(plt.MaxNLocator(10))
for label in ax2.get_xticklabels():
        label.set_rotation(90)
ax2.set_xlabel('Date')
ax2.set_ylabel('Weekly COVID-19 Deaths per 100,000')
fig.tight_layout()
plt.savefig('data/out/counts.png')
plt.savefig('data/out/counts.svg')

del ax1, ax2, fig, label