#################################################################################
# check_network.py
# reads species/reaction network files and checks it
#################################################################################

from scipy.interpolate import interp1d
from numpy import array, arange, sin
import numpy as np
import math as math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import requests
import pandas as pd
import os
import timeit

# convert to float when string is float
def convert2float(value):
  try:
    float(value)
    return float(value)
  except:
    return value
#
# check where is the 1st float in a list
def check1rstfloat(inplist):
  #typelist = [type(value) for value in inplist]
  for i, value in enumerate(inplist):
    try:
      float(value)
      return i
    except:
      pass      
#
# create a list with the mass of each element
def mass_elem(num_elem):
  num_elem = remove_sp(num_elem)
  df_mass = pd.Series({'H' : 1., 'D' : 2., 'He' : 4., 'C' : 12., 'N' : 14., 'O' : 16., 'Si' : 28., 
                          'S' : 32., 'Fe' : 55.6, 'Na' : 23., 'Mg' : 24.3, 'Cl' : 35.4, 'P' : 31., 'F' : 19.})
  m_elem = num_elem.mul(df_mass, fill_value=1)
  mass = m_elem.sum()
  return mass
#
# remove all but element numbers in series
def remove_sp(df):
  try:
    del df['Species']
  except:
    pass
  try:
    del df['No']
  except:
    pass
  try:
    del df['+-']
  except:
    pass
  try:
    del df['Xinit']
  except:
    pass
  return df
#
# element mass "database" series that will be used to compute mass of species
db_Melem = pd.Series({'H' : 1., 'D' : 2., 'He' : 4., 'C' : 12., 'N' : 14., 'O' : 16., 'Si' : 28., 
                      'S' : 32., 'Fe' : 55.6, 'Na' : 23., 'Mg' : 24.3, 'Cl' : 35.4, 'P' : 31., 'F' : 19., 
                      '15N' : 15., '13C' : 13., '18O' : 18., '33S' : 33., '34S' : 34.})
#
# "false" species that can exist in chemical networks
false_spec = ['CR','CRP','Photon','Photons',"-"]
false_mass = [0. for i in false_spec]
df_false = pd.DataFrame({'Species' : false_spec, 'Mass' : false_mass},index=false_spec)
#
# "Standard" list of reactants and products used in df
list_re = ['React1','React2','React3','Prod1','Prod2','Prod3','Prod4']
list_re2 = list_re+['Tmin','Tmax','Type','Number']
list_checkre = ['check_'+ire for ire in list_re]

##------------------------------------------------------------------
##------------------------------------------------------------------
## read species and reaction files and create dataframes
##------------------------------------------------------------------
##------------------------------------------------------------------

start = timeit.default_timer()

filespec = 'files/sp_taquet18.in'
filereac = 'files/re_taquet18.in'
outsuf = 'taquet18_'

outsuf = 'results/'+outsuf

# check if directories exist
if os.path.isdir('results') == False:
  os.system("mkdir results")

# read species file
opsp = open(filespec,'r')
s = []
for line in opsp:
  # read column headers
  if line.find('Xinit') >= 0:
    col = line[1:].split() ; colsp = col
    elem = col[3:-1] ; Nelem = len(elem)
  # create and append dataframe
  if line[0:1] != '!':
    s2 = line.split()
    s3 = [convert2float(si) for si in s2]
    s.append(s3)
s_transp = list(map(list, zip(*s)))
df_spec = pd.DataFrame({col[i] : s_transp[:][i] for i in range(len(col))})
df_spec.index = df_spec['Species']
Nspecies = len(df_spec.index)
#
# compute mass of each species: multiply df_spec to mass series - 2 ms to run
ser_mass = pd.Series({ ie : db_Melem.loc[ie] for ie in elem}) # series with mass of each element
df_spec_sub = df_spec.loc[:,elem] # extract subset of df only with number of elements
df_spec_mass = df_spec_sub * ser_mass # df computing total mass of each element in species
df_spec['Mass'] = df_spec_mass.sum(axis=1) # new column for mass of each species
Mspec = pd.Series(df_spec['Mass'])
#
# create species series with "false" species that can exist in reaction files
ser_spec = df_spec['Species']
ser_spec = ser_spec.append(pd.Series(false_spec,index=false_spec),ignore_index=True)
df_spmass = df_spec[['Species','Mass']]
df_spmass = df_spmass.append(df_false)
#
# read reaction file
opre = open(filereac,'r')
s = [] ; il = 0
for line in opre:
  # read column headers
  if line.find('React1') >= 0 or line.find('react1') >= 0:
    col = line[1:-1].split()
    # find where each column begins in reaction file
    wherecol = [line.find(icol) for icol in col] 
    wherecol[0] += -1
    # find where products and parameters begin in reaction file
    whereA = line.find('A ') 
    whereprod = line.find('Prod1')
    if whereprod == -1:
      whereprod = line.find('prod1') 
    whereN = [i for i in range(len(col)) if "Number" in col[i]][0] - len(col)
    # compute number of parameters, reactants, and products
    Nreact = len(line[1:whereprod].split())
    Nprod = len(line[whereprod:whereA].split())
    Nmol = len(line[1:whereA].split())
    Nparam = len(line[whereA:-1].split())
  # create and append dataframe
  if line[0:1] != '!' and line[0:1] != '*':
    s2 = line.split() ; il += 1
    s2[whereN] = il
    # add blank species in list
    sreact = [line[wherecol[i]:wherecol[i+1]].replace(" ","") if line[wherecol[i]] != " " else "-" for i in range(Nreact)] 
    sprod = [line[wherecol[Nreact+i]:wherecol[Nreact+i+1]].replace(" ","") if line[wherecol[Nreact+i]] != " " else "-" for i in range(Nprod)] 
    whereA = check1rstfloat(s2) 
    sfloat = [float(s2[i+whereA]) for i in range(Nparam)] # list of parameters
    # adjust list of reactants and products
    if len(sreact) == 1:
      sreact +- ["-","-"]
    if len(sreact) == 2:
      sreact +- ["-"]
    if len(sprod) == 2:
      sprod +- ["-","-"]
    if len(sprod) == 3:
      sprod +- ["-"]
    sreact.sort(reverse=True)
    sprod.sort(reverse=True)
    #print(sreact)
    s.append(sreact+sprod+sfloat)
  if line[0:1] == '*':
    break
#
# create reaction dataframe
s_transp = list(map(list, zip(*s)))
colmodif = list_re + col[Nmol:]
if len(colmodif) != len(s_transp):
  print('Problem with number of columns in reaction file')
  exit()
dict_react = {colmodif[i] : s_transp[:][i] for i in range(len(colmodif))}
df_reac = pd.DataFrame(dict_react)
if 'Number' in df_reac.columns:
  df_reac['Number'] = df_reac['Number'].astype(int)
if 'Num' in df_reac.columns:
  df_reac['Num'] = df_reac['Num'].astype(int)
if 'Type' in df_reac.columns:
  df_reac['Type'] = df_reac['Type'].astype(int)
if 'Form' in df_reac.columns:
  df_reac['Form'] = df_reac['Form'].astype(int)
if 'Nrate' in df_reac.columns:
  df_reac['Nrate'] = df_reac['Nrate'].astype(int)
#
# analyse reaction file with some pandas functions
#print(df_reac['Num'].head())
#exit()#print(df_reac.shape)
#print(df_spec.info())
#print(df_reac["React1"].value_counts(dropna=False).head())
#print(df_reac.describe())

##------------------------------------------------------------------
##------------------------------------------------------------------
## analyse reaction dataframe
##------------------------------------------------------------------
##------------------------------------------------------------------
#
# check that there is no duplicates in species file
#
spec_count = df_spec['Species'].value_counts(dropna=False) #
if max(spec_count) > 1:
  print('Duplicated species in species file: STOP')
  df_spec.ix[spec_count > 1,colsp].to_csv(outsuf+'duplicate_species.csv')
  #print(df_spec.ix[spec_count > 1,'Species'])
  exit()
else:
  print('No duplicated species...')

for ir in list_re:
  df_reac['check_'+ir] = df_reac[ir].isin(ser_spec)
df_check = df_reac[list_checkre].all(axis=1)
#
# check that all species in reaction file exist in species file
#
for ir in list_re:
  df_reac['check_'+ir] = df_reac[ir].isin(ser_spec)
df_check = df_reac[list_checkre].all(axis=1)

# print reactions with wrong species
if all(df_check) == False:
  print('Some species the reaction file are not defined in the species file: STOP')
  df_reac.loc[~df_check, list_re2].to_csv(outsuf+'reacs_notinspecies.csv')
  #print(df_reac.loc[~df_check, list_re])
  exit()
else:
  print('All reactants and products are in species file...')
#
# check mass conservation of each reaction
#
df_remass = pd.DataFrame()
for ire in list_re:
  df_reac_2 = pd.merge(left=df_reac,right=df_spmass,left_on=[ire], right_index=True,how='left',sort=False) 
  df_remass[ire] = df_reac_2['Mass']
df_remass['dM'] = df_remass[list_re[0:3]].sum(axis=1)-df_remass[list_re[3:]].sum(axis=1)
#
# print reactions that are non conservative
if max(abs(df_remass['dM'])) != 0:
  print('Reactions are not conservative: STOP')
  #print(df_reac.ix[df_remass['dM'] != 0,list_re])
  df_reac.ix[df_remass['dM'] != 0,list_re2].to_csv(outsuf+'reacs_notconserv.csv')
  exit()
else:
  print('All reactions are conservative...')
#
# check duplicated reactions
#
# sort reactions according to reactants and then products
df_reac_ord = df_reac.sort_values(by=list_re).reset_index()
#
# select all duplicated reactions
df_reac_dupl = df_reac_ord[df_reac_ord.duplicated(subset=list_re+['Type'],keep=False)].reset_index(drop=True)
#
# select duplicated reactions with inconsistent ranges of temperatures
df_reac_dupl2 = pd.DataFrame(columns=list(df_reac_dupl.columns))
il = 0
for index, row in df_reac_dupl.iterrows():
  try:
    row2 = df_reac_dupl.loc[index+1]
    if row[list_re].equals(row2[list_re]) == True and \
      ((row['Tmax'] <= row2['Tmax'] and row['Tmin'] >= row2['Tmin']) or \
            (row2['Tmax'] <= row['Tmax'] and row2['Tmin'] >= row['Tmin'])):
      il+=1
      df_reac_dupl2 = df_reac_dupl2.append(row)
      df_reac_dupl2 = df_reac_dupl2.append(row2)
  except:
    break            
#
# re-sort dulication reactions according to initial number
df_reac_dupl2 = df_reac_dupl2.sort_values(by=['Number']).reset_index(drop=True)
#
# write duplicated reactions in csv file
df_reac_dupl2[col].to_csv(outsuf+'duplicate_reactions.csv')
print('Duplicated reactions are listed in output file...')


stop = timeit.default_timer()
print('Time: ', stop - start)  
exit()






