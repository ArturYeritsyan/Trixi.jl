#!/usr/bin/python3
# Note: L2-Error (not L1)

import matplotlib.pyplot as plt
import matplotlib
import numpy as np

N = [4, 8, 16, 32, 64]

# Relaxed, k = 4, Random level assignment
# SSP64_2Nstar
L2Err_B1_SSP64 = [2.7947277388634555e-5, 6.794247613059925e-7, 2.0338129173678366e-8, 6.447083871620702e-10, 2.0773903147129062e-11]

# SSP104_2Nstar
L2Err_B1_SSP104 = [4.151092463630538e-5, 2.1302619331830362e-6, 1.2801513385936967e-7, 7.948576603281711e-9, 4.956373550694869e-10]

# SSP144_2Nstar
L2Err_B1_SSP144 = [4.5620036704879025e-5, 6.7672731467602716e-6, 4.213414225436999e-7, 2.637972827480107e-8, 1.6494654089996852e-9]


GoldenRatio = (1 + 5 ** 0.5) / 2

InchesX = 4.33071 # = 11 cm = 0.75 * textwidth (in current document)
InchesY = InchesX / GoldenRatio
#fig, ax = plt.subplots(figsize = (InchesX, InchesY) )
fig, ax = plt.subplots()

RWTH_Blue_RGB   = [(0, 84/256, 159/256)]
RWTH_Orange_RGB = [(246/256, 169/256, 0)]
RWTH_Green_RGB  = [(70/256, 171/256, 39/256)]
RWTH_Red_RGB    = [(204/256, 7/256, 30/256)]
RWTH_Petrol_RGB = [(0/256, 156/256, 161/256)]
RWTH_Yellow_RGB = [(225/256, 237/256, 0/256)]
RWTH_Purple_RGB = [(97/256, 33/256, 88/256)]

### ACTUAL PLOTTING: L1 errors ###

ax.scatter(N, L2Err_B1_SSP64, label = r'2N$^*$(6,4)', color = RWTH_Blue_RGB)
ax.plot(N, L2Err_B1_SSP64, color = RWTH_Blue_RGB[0], linestyle='dashed')

ax.scatter(N, L2Err_B1_SSP104, label = r'2N$^*$(10,4)', color = RWTH_Green_RGB)
ax.plot(N, L2Err_B1_SSP104, color = RWTH_Green_RGB[0], linestyle='dashed')

ax.scatter(N, L2Err_B1_SSP144, label = r'2N$^*$(14,4)', color = RWTH_Orange_RGB)
ax.plot(N, L2Err_B1_SSP144, color = RWTH_Orange_RGB[0], linestyle='dashed')


# ax.loglog(N, np.multiply(6e-2, np.power(np.array(N, dtype=float), -2) ), linestyle='dotted',
#           label = r'$\mathcal{O}\left(N^{-2}\right)$',
#           color = 'black') # Order two line fitted

ax.loglog(N, np.multiply(6e-3, np.power(np.array(N, dtype=float), -4) ), linestyle='dashdot',
          label = r'$\mathcal{O}\left(N^{-4}\right)$',
          color = 'black') # Order two line fitted

ax.loglog(N, np.multiply(7e-3, np.power(np.array(N, dtype=float), -5) ), linestyle='dotted',
            label = r'$\mathcal{O}\left(N^{-5}\right)$',
            color = 'black') # Order two line fitted

# Turn on logscale (no native support for logarithmic scatter)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel(r'$N$')

### GRID SECTION ###
ax.grid(axis ='both', which='major', alpha=0.1, linewidth = 1.5, color ='black')
ax.set_axisbelow(True)  # Hide grid behind bars

### LEGEND SECTION ###
ax.legend(loc = "lower left")

### TICKS SECTION ###

ax.set_xticks(N)

ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

ax.get_xaxis().set_tick_params(which='minor', size=0)
ax.get_xaxis().set_tick_params(which='minor', width=0) 

#ax.set_xticklabels([r"$6$", r"$12$", r"$24$", r"$48$", r"$96$", r"$192$"])

### LIMITS SECTION ###
eps_x = 0.1
ax.set_xlim([4 * (1 - eps_x), 64 * (1 + eps_x)])

# Make bounding lines thicker
'''
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2)
'''

### TITLE SECTION ###
plt.title(r"$L^2_{error}$ (Alfven Wave)")

# Scale while preserving aspect ratio
width, height = fig.get_size_inches()
factor = InchesX / width
fig.set_size_inches(width * factor, height * factor)

plt.tight_layout()  # Title, labels, ... to screen
plt.savefig('L2_error.png', dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()
