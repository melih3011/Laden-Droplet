#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import scienceplots
from sys import argv


def get_radius(rho):

  R = 0
  while (abs((rho[R]-rho[-1])/rho[-1]) > 5.e-2):
    R += 1

  return R

def fr2int(rho, dr):

  integral = 0.0
  R = get_radius(rho)
  print('R = {:.3f}'.format(R*dr))

  dr *= 2.
  hr = dr*dr
  for i in range(R):
    integral += rho[i]*hr*(1.+4.*i+6.*i*i)+rho[i+1]*hr*(3.+8.*i+6.*i*i)

  return 4.*np.pi*integral*dr/12., np.pi/6.*3./(R*0.5*dr)**3*integral*dr/12.

def main():


  data = np.loadtxt(argv[1])
  
  r_vals = data[:, 0]
  dr = r_vals[1]-r_vals[0]
  
  rho_l = data[:, 1]
  rho_n = data[:, 2]

  print(fr2int(rho_l, dr))
  print(fr2int(rho_n, dr))


if __name__=='__main__':


  main()
