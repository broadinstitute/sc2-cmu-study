#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 10:16:34 2021

@author: bpetros
"""
#%%

## Download relevant packages and set up environment
from adjustText import adjust_text
import itertools
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scipy.stats as stats
import seaborn as sns
import statsmodels.api as sm
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

# Residence hall incidence rates as function of hall parameters
X = halls[['hall_population', 'occupancy', 'floors', 'dining',\
           'private_bath', 'RA', 'sqft', 'height', 'volume_per_person']]
y = halls['hall_cases']/halls['hall_population']


# Assess for multicollinearity (adapted from https://aegis4048.github.io/mutiple_linear_regression_and_visualization_in_python)
corr = X.corr(method='spearman')
corr.rename({'hall_population': 'Number of Students', 'occupancy': 'Percent Occupied',
             'floors': 'Floors', 'dining': 'Dining Hall', 'private_bath': 'In-Unit Bathroom',
             'RA': 'Number of RAs', 'sqft': 'Square Footage', 'height': 'Ceiling Height',
             'volume_per_person': 'Volume per Person'}, axis='index', inplace = True)
corr.rename({'hall_population': 'Number of Students', 'occupancy': 'Percent Occupied', 
             'floors': 'Floors', 'dining': 'Dining Hall', 'private_bath': 'In-Unit Bathroom',
             'RA': 'Number of RAs', 'sqft': 'Square Footage', 'height': 'Ceiling Height',
             'volume_per_person': 'Volume per Person'}, axis='columns', inplace = True)

# Mask upper triangle
mask = np.zeros_like(corr, dtype=np.bool)
mask[np.triu_indices_from(mask)] = True
np.fill_diagonal(mask, False)
# Generate custom diverging colormap
fig, ax = plt.subplots(figsize=(6, 5))
cmap = sns.color_palette("Greys")
# Draw heatmap with mask and aspect ratio
sns.heatmap(corr, mask=mask, cmap=cmap, vmin=-1, vmax=1, center=0, linewidths=.5)
ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
fig.suptitle('Correlation Matrix of Predictors', fontsize=16)
fig.tight_layout()
plt.savefig('data/out/multicollinearity.png', bbox_inches='tight', dpi = 800)
plt.savefig('data/out/multicollinearity.svg')

del ax, cmap, corr, fig, mask

#%%

# Run multiple regression on all possible combinations of input variables
combo = list(itertools.combinations(X.columns.tolist(), 1))
combo.extend(list(itertools.combinations(X.columns.tolist(), 2)))
combo.extend(list(itertools.combinations(X.columns.tolist(), 3)))
combo.extend(list(itertools.combinations(X.columns.tolist(), 4)))
combo.extend(list(itertools.combinations(X.columns.tolist(), 5)))
combo.extend(list(itertools.combinations(X.columns.tolist(), 6)))
combo.extend(list(itertools.combinations(X.columns.tolist(), 7)))
combo.extend(list(itertools.combinations(X.columns.tolist(), 8)))
combo.extend(list(itertools.combinations(X.columns.tolist(), 9)))

aic = []
bic = []

for i in range(len(combo)):
    model = sm.OLS(y, X[list(combo[i])])
    results = model.fit()
    aic.append(results.aic)
    bic.append(results.bic)
    del model, results
del i
    
# Select model that minimizes the BIC
model = sm.OLS(y, X[list(combo[bic.index(min(bic))])])
results = model.fit()
print(results.summary())
rmse = np.sqrt(np.mean(results.resid**2))

#%%

# Actual vs. predicted hall incidence rates
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), dpi = 800)
sns.despine()
ax1.scatter(results.predict(X[list(combo[bic.index(min(bic))])]), y, color = 'black')
ax1.plot(np.arange(0.5, 0.3, 0.05), np.arange(0.5, 0.3, 0.05), linestyle = 'dashed', color = 'black', alpha = 0.5)
r, p = stats.pearsonr(results.predict(X[list(combo[bic.index(min(bic))])]), y)
labs = []
for i in range(len(halls.index)): 
        labs.append(halls.index[i].split(' ')[1])
texts = []
for x, yi, text in zip(results.predict(X[list(combo[bic.index(min(bic))])]), y, labs):
    texts.append(ax1.text(x+0.002, yi+0.005, text, size = 12))
ax1.set_xlabel('Predicted Incidence', fontsize = 14)
ax1.set_ylabel('Actual Incidence', fontsize = 14)

del r, p, x, yi, text, texts

# Residual plot: residuals vs. predicted hall incidence rates
ax2.scatter(results.predict(X[list(combo[bic.index(min(bic))])]), results.resid, color = 'black')
ax2.plot([0.1, 0.225], [0, 0], linestyle = 'dashed', color = 'black')
r, p = stats.pearsonr(results.predict(X[list(combo[bic.index(min(bic))])]), results.resid)
texts = []
for x, yi, text in zip(results.predict(X[list(combo[bic.index(min(bic))])]), results.resid, labs):
    texts.append(ax2.text(x, yi, text, size = 12))
adjust_text(texts, force_points=(1, 1), arrowprops=dict(arrowstyle="-",
              color='black', alpha=0.5), expand_points=(2,2))
ax2.set_xlabel('Predicted Incidence', fontsize = 14)
ax2.set_ylabel('Residuals', fontsize = 14)
fig.tight_layout()
plt.savefig('data/out/regression.png', bbox_inches='tight', dpi = 800)
plt.savefig('data/out/regression.svg')

# Save residuals for downstream
res = results.resid
del ax1, ax2, fig, i, model, r, results, rmse, p, x, yi, text, texts

#%%

## Leave-one-out cross-validation (LOOCV)
pred = []
for hall in X.index:
    Xh = X.drop(hall)
    yh = y.drop(hall)
    Xk = X[X.index == hall]
    model = sm.OLS(yh, Xh[list(combo[bic.index(min(bic))])])
    results = model.fit()
    pred.extend(results.predict(Xk[list(combo[bic.index(min(bic))])]).tolist())

halls['pred'] = pred
del hall, pred

# LOOCV residuals vs. predicted value
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6))
sns.despine()
ax1.scatter(halls.pred, y - halls.pred, color = 'black')
line = ax1.plot(np.arange(0.1, 0.3, 0.05), [0] * 4, linestyle = 'dashed', color = 'black', alpha = 0.5)
texts = []
for x, yi, text in zip(halls.pred, y - halls.pred, labs):
    texts.append(ax1.text(x, yi+0.005, text, size = 12))
ax1.set_xlabel('Predicted Incidence', fontsize = 16)
ax1.set_ylabel('Residuals', fontsize = 16)
ax1.tick_params(axis='both', which='major', labelsize=12)
ax1.set_title('N-1 Model Fit', fontsize = 18)

del line, x, yi, text, texts

# LOOCV residuals vs. model residuals
ax2.scatter(res, y - halls.pred, color = 'black')
texts = []
for x, yi, text in zip(res, y - halls.pred, labs):
    texts.append(ax2.text(x+0.001, yi+0.001, text, size = 12))
adjust_text(texts, arrowprops=dict(arrowstyle="-", color='black', alpha=0.5), \
           force_points = (0.3, 0.3), expand_points= (0.1, 0.1)) 
ax2.set_xlabel('Residuals, N Model Fit', fontsize = 16)
ax2.set_ylabel('Residuals, N-1 Model Fit', fontsize = 16)
ax2.set_title('N vs. N-1 Residuals', fontsize = 18)
ax2.tick_params(axis='both', which='major', labelsize=12)
plt.savefig('data/out/leave_one_out.png', bbox_inches='tight', dpi = 800)
plt.savefig('data/out/leave_one_out.svg')

rmse = np.sqrt(np.mean((halls.pred - y)**2))

del aic, ax1, ax2, bic, combo, fig, labs, model, res, results, X, Xh, Xk, x, y, yh, yi


