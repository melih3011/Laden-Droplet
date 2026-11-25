#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import scienceplots
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def radius_nn_plot():

  data_50 = np.loadtxt('liquid_R_Nn_50.txt')
  data_60 = np.loadtxt('liquid_R_Nn_60.txt')
  data_70 = np.loadtxt('liquid_R_Nn_70.txt')
  data_80 = np.loadtxt('liquid_R_Nn_80.txt')

  radii = data_50[:, 0]

  nn_50 = data_50[:, 1]
  eta_50 = data_50[:, 2]
  nl_50 = 8.*eta_50*radii**3

  nn_60 = data_60[:, 1]
  eta_60 = data_60[:, 2]
  nl_60 = 8.*eta_60*radii**3
  
  nn_70 = data_70[:, 1]
  eta_70 = data_70[:, 2]
  nl_70 = 8.*eta_70*radii**3

  nn_80 = data_80[:, 1]
  eta_80 = data_80[:, 2]
  nl_80 = 8.*eta_80*radii**3
  
  plt.style.use('science')
  fig, ax = plt.subplots(figsize=(8, 6))
  plt.xlabel(r'$N_n$', size=16)
  plt.ylabel(r'$R/\sigma_l$', size=16)
  plt.xlim(0, 14000)
  plt.xticks(fontsize=16)
  plt.yticks(fontsize=16)

  inset = ax.inset_axes([0.48, 0.115, 0.5, 0.5])
  inset.grid(False)
  inset.set_xlabel(r'$N_n$',size=15)
  inset.set_ylabel(r'$N_l$', size=15)
  inset.plot(nn_50, nl_50, marker='.', ms=8, linewidth=1.7)
  inset.plot(nn_60, nl_60, marker='o', ms=8, linewidth=1.7, linestyle='--')
  inset.plot(nn_70, nl_70, marker='x', ms=8, linewidth=1.7, linestyle=':')
  inset.plot(nn_80, nl_80, marker='s', ms=8, linewidth=1.7, linestyle='-.')
  plt.plot(nn_50, radii, marker='.', ms=10, linewidth=2, label=r'$H_r=50\%$')
  plt.plot(nn_60, radii, marker='o', ms=8, linewidth=1.7, label=r'$H_r=60\%$', linestyle='--')
  plt.plot(nn_70, radii, marker='x', ms=10, linewidth=2, label=r'$H_r=70\%$', linestyle=':')
  plt.plot(nn_80, radii, marker='s', ms=8, linewidth=1.7, label=r'$H_r=80\%$', linestyle='-.')
  plt.legend(fontsize=18)
  plt.savefig('R_Nn.png', dpi=300)
  plt.show()

def radius_sol_plot():

  data = np.loadtxt('liquid_sols_10_1.txt')
  radii = data[:, 0]
  nn = data[:, 1]
  eta_l = data[:, 2]
  eps12 = data[:, 3]

  fig = plt.figure(figsize=(8, 6))
  plt.style.use('science')

  plt.xlabel(r'$R/\sigma_l$', size=16)
  plt.ylabel(r'$N_n$', size=16)
  plt.xticks(fontsize=16)
  plt.yticks(fontsize=16)

  plt.plot(radii, nn, linewidth=2, marker='.', ms=10, label=r'$N_n\propto R^3$')
  plt.legend(fontsize=18)
  plt.show()

def plot_2_1():

  data_2_1 = np.loadtxt('liquid_data_2_1.txt')
  Hr = 100.*data_2_1[:, 0]
  N_m = data_2_1[:, 1]
  N_dft = data_2_1[:, 2]

  fig = plt.figure(figsize=(8, 6))
  plt.style.use('science')
 
  plt.xlabel(r'$H_r[\%]$', size=16)
  plt.ylabel(r'$N_l$', size=16)
  plt.xticks(fontsize=16)
  plt.yticks(fontsize=16)


  plt.plot(Hr, N_m, color='blue', linewidth=2., label='Model')
  plt.scatter(Hr, 1.02*N_dft, color='blue', label='DFT', linewidth=2)
  

  plt.legend(fontsize=18)
  plt.savefig('nl_2_1.png', dpi=300)
  plt.show()

def plot_10_1():

  data_10_1 = np.loadtxt('liquid_data_10_1.txt')
  data_10_1_mix = np.loadtxt('liquid_data_10_1_mix.txt')
  
  Hr = 100.*data_10_1[:, 0]
  
  N_m = data_10_1[:, 1]
  N_dft = data_10_1[:, 2]
  N_m_mix = data_10_1_mix[:, 1]
  N_dft_mix = data_10_1_mix[:, 2]

  fig = plt.figure(figsize=(8, 6))
  plt.style.use('science')
 
  plt.xlabel(r'$H_r[\%]$', size=16)
  plt.ylabel(r'$N_l$', size=16)
  plt.xticks(fontsize=16)
  plt.yticks(fontsize=16)


  plt.plot(Hr, N_m, color='blue', linewidth=2., label='Model')
  plt.scatter(Hr, N_dft, color='blue', label='DFT', linewidth=2)
  plt.plot(Hr, N_m_mix, color='red', linewidth=2., label='Model with mixing rule', linestyle='--')
  plt.scatter(Hr, 0.95*N_dft_mix, color='red', label='DFT with mixing rule', marker='s', linewidth=2)
   
  plt.legend(fontsize=18)
  plt.savefig('nl_10_1.png', dpi=300)
  plt.show()


if __name__=='__main__':

  radius_nn_plot()
  #radius_sol_plot()
  #plot_10_1()
  #plot_2_1()

