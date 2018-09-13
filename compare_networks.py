#################################################################################
# compare_networks.py
# compare two chemical networks, each containing a species and a reactions file
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

##------------------------------------------------------------------
##------------------------------------------------------------------
## read species and reaction files and create dataframes
##------------------------------------------------------------------
##------------------------------------------------------------------

start = timeit.default_timer()

network = ['taquet18', 'taquet18_2']
filespec = ['files/sp_'+net+'.in' for net in network]
filereac = ['files/re_'+net+'.in' for net in network]

outsuf = network[0]+'_'+network[1]
outsuf = 'results/'+outsuf

# check if directories exist
if os.path.isdir('results') == False:
  os.system("mkdir results")
#
# "Standard" list of reactants and products used in df
list_react = ['React1','React2','React3','Prod1','Prod2','Prod3','Prod4'] 
list_abc = ['A','B','C']
list_type = ['Type'] 
list_temp = ['Tmin','Tmax']
list_othparam = ['Form','Nrate']
list_allparams = list_react + list_abc + list_type + list_othparam
list_abc2 = [param+'_'+net for net in network for param in list_abc]
list_temp2 = [temp+'_'+net for net in network for temp in list_temp]
list_checkre = ['check_'+ire for ire in list_react]
#
Nspecies = [] ; Mspec = [] ; df_spec = [] ; df_spec_mass = [] ; df_spmass = []
ser_mass = [] ; ser_spec = []
Nreact = [] ; Nprod = [] ; Nmol = [] ; Nparam = []
col = [] ; whereA = [] ; whereprod = [] ; whereN = []
df_reac = [] 
#
# read species files
for ifi in range(2):

  opsp = open(filespec[ifi],'r')
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
  df_spec.append(pd.DataFrame({col[i] : s_transp[:][i] for i in range(len(col))}))
  df_spec[ifi].index = df_spec[ifi]['Species']
  Nspecies.append(len(df_spec[ifi].index))
  #
  # compute mass of each species: multiply df_spec to mass series - 2 ms to run
  ser_mass.append(pd.Series({ ie : db_Melem.loc[ie] for ie in elem})) # series with mass of each element
  df_spec_sub = df_spec[ifi].loc[:,elem] # extract subset of df only with number of elements
  df_spec_mass.append(df_spec_sub * ser_mass[ifi]) # df computing total mass of each element in species
  df_spec[ifi]['Mass'] = df_spec_mass[ifi].sum(axis=1) # new column for mass of each species
  Mspec.append(pd.Series(df_spec[ifi]['Mass']))
  #
  # create species series with "false" species that can exist in reaction files
  ser_spec.append(df_spec[ifi]['Species'])
  ser_spec[ifi] = ser_spec[ifi].append(pd.Series(false_spec,index=false_spec),ignore_index=True)
  df_spmass.append(df_spec[ifi][['Species','Mass']])
  df_spmass[ifi] = df_spmass[ifi].append(df_false)
  #
  # read reaction file
  opre = open(filereac[ifi],'r')
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
      whereprod.append(line.find('Prod1'))
      if whereprod[ifi] == -1:
        whereprod[ifi] = line.find('prod1') 
      whereN.append([i for i in range(len(col)) if "Number" in col[i]][0] - len(col))
      # compute number of parameters, reactants, and products
      Nreact.append(len(line[1:whereprod[ifi]].split()))
      Nprod.append(len(line[whereprod[ifi]:whereA].split())) 
      Nmol.append(len(line[1:whereA].split()))
      Nparam.append(len(line[whereA:-1].split()))
    # create and append dataframe
    if line[0:1] != '!' and line[0:1] != '*':
      s2 = line.split() ; il += 1
      #print(len(s2))
      s2[whereN[ifi]] = il
      # add blank species in list
      sreact = [line[wherecol[i]:wherecol[i+1]].replace(" ","") if line[wherecol[i]] != " " else "-" for i in range(Nreact[ifi])] 
      sprod = [line[wherecol[Nreact[ifi]+i]:wherecol[Nreact[ifi]+i+1]].replace(" ","") if line[wherecol[Nreact[ifi]+i]] != " " else "-" for i in range(Nprod[ifi])] 
      whereA = check1rstfloat(s2) 
      sfloat = [float(s2[i+whereA]) for i in range(Nparam[ifi])] # list of parameters
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
  s_transp = list(map(list, zip(*s))) # transpose reaction list of lists
  colmodif = list_react + col[Nmol[ifi]:] # 
  if len(colmodif) != len(s_transp):
    print('Problem with number of columns in reaction file')
    exit()
  dict_react = {colmodif[i] : s_transp[:][i] for i in range(len(colmodif))}
  df_reac.append(pd.DataFrame(dict_react)) 
  if 'Number' in df_reac[ifi].columns:
    #print(df_reac.info())
    #exit()
    df_reac[ifi]['Number'] = df_reac[ifi].index.values+1
    df_reac[ifi]['Number'] = df_reac[ifi]['Number'].astype(int)
  if 'Num' in df_reac[ifi].columns:
    df_reac[ifi]['Num'] = df_reac[ifi].index.values+1
    df_reac[ifi]['Num'] = df_reac[ifi]['Num'].astype(int)
  if 'Type' in df_reac[ifi].columns:
    df_reac[ifi]['Type'] = df_reac[ifi]['Type'].astype(int)
  if 'Form' in df_reac[ifi].columns:
    df_reac[ifi]['Form'] = df_reac[ifi]['Form'].astype(int)
  if 'Nrate' in df_reac[ifi].columns:
    df_reac[ifi]['Nrate'] = df_reac[ifi]['Nrate'].astype(int)

##------------------------------------------------------------------
##------------------------------------------------------------------
## check similar reactions
##------------------------------------------------------------------
##------------------------------------------------------------------
#
# order reactions according to reactants and products
df_reac_ord = []
df_reac_ord.append(df_reac[0].sort_values(by=list_react).reset_index())
df_reac_ord.append(df_reac[1].sort_values(by=list_react).reset_index())
#
# merge two dfs according to reactants, products, reaction types
list_commonparams = list_react+list_type+list_othparam
df_reac_merge = df_reac_ord[0].merge(df_reac_ord[1],how='outer',on=list_commonparams,suffixes=['_'+network[0],'_'+network[1]])
df_reac_merge_ord = df_reac_merge[list_commonparams+list_abc2+list_temp2+['Number_'+network[0]]].sort_values(by=['Number_'+network[0]]).reset_index()
#df_reac_merge_ord.to_csv(outsuf+'_reac_merge_ord.csv')
#
# check if A, B, C parameters are equal
for ia, abc in enumerate(list_abc): 
  check1 = abc+'_'+network[0] 
  check2 = abc+'_'+network[1]
  df_reac_merge_ord['check_'+abc] = df_reac_merge_ord[check1] == df_reac_merge_ord[check2]
df_reac_merge_ord['check_abc'] = df_reac_merge_ord[['check_A','check_B','check_C']].all(axis=1)
#
# check if temperature ranges are equal
for ia, T in enumerate(list_temp): 
  check1 = T+'_'+network[0] 
  check2 = T+'_'+network[1]
  df_reac_merge_ord['check_'+T] = df_reac_merge_ord[check1] == df_reac_merge_ord[check2]
df_reac_merge_ord['check_temp'] = df_reac_merge_ord[['check_Tmin','check_Tmax']].all(axis=1)
#
# select reactions that appear in both chemical networks
df_reac_merge_sub = df_reac_merge_ord[~pd.isnull(df_reac_merge_ord['A_'+network[0]])]
df_reac_merge_sub = df_reac_merge_sub[~pd.isnull(df_reac_merge_sub['A_'+network[1]])]
#
# select reactions that appear in one of the two chemical networks
df_reac_reac1 = df_reac_merge_ord[pd.isnull(df_reac_merge_ord['A_'+network[1]])]
df_reac_reac2 = df_reac_merge_ord[pd.isnull(df_reac_merge_ord['A_'+network[0]])]
#
# select reactions with different A, B, C parameters or different temperature ranges
df_reac_diffabc = df_reac_merge_sub[df_reac_merge_sub['check_abc'] == False]
df_reac_difftemp = df_reac_merge_sub[df_reac_merge_sub['check_temp'] == False]
#
# write results in csv files
listre_check = list_commonparams+list_abc2+list_temp2
df_reac_diffabc[listre_check].to_csv(outsuf+'_differentabc_reactions.csv')
df_reac_difftemp[listre_check].to_csv(outsuf+'_differenttemp_reactions.csv')
df_reac_reac1[listre_check].to_csv(outsuf+'_only'+network[0]+'_reactions.csv')
df_reac_reac2[listre_check].to_csv(outsuf+'_only'+network[1]+'_reactions.csv')


stop = timeit.default_timer()
print('Time: ', stop - start)  
exit()
