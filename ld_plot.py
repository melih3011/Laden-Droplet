#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from sys import argv
import scienceplots

def centeredDistanceMatrix(n):
  
  x, y = np.meshgrid(np.arange(n), np.arange(n))
  return np.sqrt((x - (n // 2))**2 + (y - (n // 2))**2)

def arbitraryfunction(d, y):

  n = len(y)
  x = np.arange(n)
  f = interp1d(x, y, bounds_error=False, fill_value="none")
  
  return f(d.flat).reshape(d.shape)

def read_datafile(data):
  
  arr0, arr1, arr2 =  [], [], []
  with open(data, 'r') as dat:
    for line in dat:
      if (line.startswith('#')):
        continue
      else:
        line_split = line.split()
        arr0.append(float(line_split[0]))
        arr1.append(float(line_split[1]))
        arr2.append(float(line_split[2]))
  
  return [np.array(arr0), 8.*np.array(arr1), 8.*np.array(arr2)]

def plotter(offset=0):
  
  batch = []
  #im = plt.imread(get_sample_data('./50%_300.png'))
  
  for i in range(1, len(argv)):
    batch.append(read_datafile(str(argv[i])))
    
  rvals = batch[0][0]
  xmax = 200.#np.max(rvals)
  xmax_ind = int(1.2*xmax/(rvals[1]-rvals[0]))
  
  plt.style.use('science')
  fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))
  
  axs[0].set_xlabel(r'$r/\sigma_l$', fontsize=25)
  axs[0].set_ylabel(r'$\rho_l(r)\sigma_l^3$', fontsize=25)
  axs[0].tick_params(axis='both', which='major', labelsize=20)

  axs[1].set_xlabel(r'$r/\sigma_l$', fontsize=25)
  axs[1].set_ylabel(r'$\rho_n(r)\sigma_l^3$', fontsize=25) 
  axs[1].tick_params(axis='both', which='major', labelsize=20)
  
  d = centeredDistanceMatrix(xmax_ind)
  heatmapy0 = arbitraryfunction(d, batch[len(argv)-2][1])
  heatmapy1 = arbitraryfunction(d, batch[len(argv)-2][2])
  
  inset0 = axs[0].inset_axes([0.52, 0.49, 0.46, 0.46])
  inset0.grid(False)
  inset0.set_xticks([])
  inset0.set_yticks([])
  inset0.imshow(heatmapy0.T, origin='lower', interpolation='none', cmap='viridis')

  inset1 = axs[1].inset_axes([0.52, 0.49, 0.46, 0.46])
  inset1.grid(False)
  inset1.set_xticks([])
  inset1.set_yticks([])
  inset1.imshow(heatmapy1.T, origin='lower', interpolation='none', cmap='viridis')  
  
  
  for i in range(len(argv)-1):
    axs[0].plot(rvals, batch[i][1], linewidth=2.2)
    axs[1].plot(rvals, batch[i][2]+i*offset, linewidth=2.2)
    
    if (offset != 0):
      axs[1].hlines(batch[i][2][-1]+i*offset, xmin=0.0, xmax=xmax,linewidth=2.2, alpha=0.3)
  
  axs[0].set_xlim(0, xmax)
  axs[1].set_xlim(0, xmax)
  axs[0].set_ylim(bottom=0)
  axs[1].set_ylim(bottom=0)
  
  fig.tight_layout()
  plt.savefig('current_file.png', dpi=600)
  #plt.show()

def plotter_2():
  
  batch = []
  colors = ['blue', 'red', 'green', 'orange']
  
  for i in range(1, len(argv)):
    batch.append(read_datafile(str(argv[i])))
    
  rvals = batch[0][0]
  xmax = 120.#np.max(rvals)
  
  fig, axs = plt.subplots(
    nrows=len(argv)-1, ncols=2, figsize=(4, 4), 
    sharex=True, gridspec_kw={'hspace': 0.06, 'wspace': 0.06}
  )
  plt.style.use('bmh')

  axs[0, 0].set_title(r'Liquid $\rho_l(r)\sigma_l^3$', fontsize=17)
  axs[0, 1].set_title(r'Nanoparticles $\rho_n(r)\sigma_l^3$', fontsize=17)
  
  for i in range(len(argv)-1):
    
    axs[i, 0].set_xlabel(r'$r/\sigma_l$', fontsize=17)
    axs[i, 1].set_xlabel(r'$r/\sigma_l$', fontsize=17)
    axs[i, 1].yaxis.tick_right()  
    
    axs[i, 0].plot(rvals, batch[i][1], linewidth=2.7, color=colors[i])
    axs[i, 1].plot(rvals, batch[i][2], linewidth=2.7, color=colors[i])
    
    axs[i, 0].set_xlim(0, xmax)
    axs[i, 0].set_ylim(bottom=0)
    axs[i, 1].set_xlim(0, xmax)  
    axs[i, 1].set_ylim(bottom=0)

  fig.tight_layout()
  plt.show()

def plotter_3():
  
  batch = []
  #im = plt.imread(get_sample_data('./50%_300.png'))
  
  for i in range(1, len(argv)):
    batch.append(read_datafile(str(argv[i])))
    
  rvals = batch[0][0]
  xmax = 200.#np.max(rvals)
  xmax_ind = int(1.2*xmax/(rvals[1]-rvals[0]))
  
  plt.style.use('science')
  fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))
  
  axs[0].set_xlabel(r'$r/\sigma_l$', fontsize=25)
  axs[0].set_ylabel(r'$\rho_l(r)\sigma_l^3$', fontsize=25)

  axs[1].set_xlabel(r'$r/\sigma_l$', fontsize=25)
  axs[1].set_ylabel(r'$\rho_n(r)\sigma_l^3$', fontsize=25) 
  
  d = centeredDistanceMatrix(xmax_ind)
  heatmapy0 = arbitraryfunction(d, batch[len(argv)-2][1])
  heatmapy1 = arbitraryfunction(d, batch[len(argv)-2][2])
 
  
  gs = fig.add_gridspec(1, 2, width_ratios=[1, 1])
  gs_inset = gs[1].subgridspec(4, 1, wspace=0, hspace=0)
  
  for i in range(len(argv)-1):
    axs[0].plot(rvals, batch[i][1], linewidth=2.2)  
  
  axes_insets = []
  for i in range(len(argv)-1):
    ax = fig.add_subplot(gs_inset[i])
    ax.plot(rvals, batch[i][2])
    ax.set_xlim(0, xmax)
    ax.set_xticks([])
    #ax.set_yticks([])
    ax.grid(False)

    axes_insets.append(ax)

  inset0 = axs[0].inset_axes([0.52, 0.49, 0.46, 0.46])
  inset0.grid(False)
  inset0.set_xticks([])
  inset0.set_yticks([])
  inset0.imshow(heatmapy0.T, origin='lower', interpolation='none', cmap='viridis')

  inset1 = axs[1].inset_axes([0.52, 0.49, 0.46, 0.46])
  inset1.grid(False)
  inset1.set_xticks([])
  inset1.set_yticks([])
  inset1.imshow(heatmapy1.T, origin='lower', interpolation='none', cmap='viridis') 


  axs[0].set_xlim(0, xmax)
  axs[1].set_xlim(0, xmax)
  axs[0].set_ylim(bottom=0)
  axs[1].set_ylim(bottom=0)
  axs[1].set_yticks([])
  
  fig.tight_layout()
  plt.savefig('current_file.png', dpi=600)
 
def main():
  
  plotter(offset=0.0)

if __name__=='__main__':
  
  main()
