import numpy as np
import matplotlib.pyplot as plt
import scienceplots

data_80 = np.loadtxt('ld_1.00_400.00_-5.4902')
data_50 = np.loadtxt('ld_1.00_400.00_-5.9402')
data_30 = np.loadtxt('ld_1.00_400.00_-6.4374')


fig = plt.figure(figsize=(8, 6))
plt.style.use('science')

plt.xlabel(r'$r/\sigma_l$', size=18)
plt.ylabel(r'$\rho(r)\sigma_l^3$', size=18)
plt.xlim(0, 250)
plt.ylim(0, 0.75)
plt.xticks(size=16)
plt.yticks(size=16)

plt.plot(data_80[:, 0], data_80[:, 1]*8., linewidth=1.5, color='blue', label=r'$H_r=80\%$')
plt.plot(data_80[:, 0], data_80[:, 2]*8., linewidth=1.5, color='blue', linestyle='--')
plt.plot(data_50[:, 0], data_50[:, 1]*8., linewidth=1.5, color='green', label=r'$H_r=50\%$')
plt.plot(data_50[:, 0], data_50[:, 2]*8., linewidth=1.5, color='green', linestyle='--')
plt.plot(data_30[:, 0], data_30[:, 1]*8., linewidth=1.5, color='red', label=r'$H_r=30\%$')
plt.plot(data_30[:, 0], data_30[:, 2]*8., linewidth=1.5, color='red', linestyle='--')

plt.legend(fontsize=17)
plt.savefig('1_to_1_radial.pdf', dpi=300)
plt.show()


