#!/usr/bin/python3
# Note: L2-Error (not L1)

import matplotlib.pyplot as plt
import matplotlib
import numpy as np

N = [4, 8, 16, 32, 64]

# Relaxed, Random level assignment
# k = 4
# SSP64_2Nstar
LInfErr_B1_SSP64 = [0.0001266430751047931, 5.013041192492018e-6, 1.6502627298020656e-7, 5.3581333814278764e-9, 1.6742418562643024e-10]

# SSP104_2Nstar
LInfErr_B1_SSP104 = [0.00012062587793193469, 5.810327280864058e-6, 2.7429463078654237e-7, 1.4168390793933838e-8, 7.936662438368103e-10]

# SSP144_2Nstar
LInfErr_B1_SSP144 = [0.00022055250733599152, 1.233172483017242e-5, 6.883720495842738e-7, 4.0217435492984066e-8, 2.4251880503811662e-9]

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

### ACTUAL PLOTTING: LInf errors ###

ax.scatter(N, LInfErr_B1_SSP64, label = r'$2N^*(6,4)$', color = RWTH_Blue_RGB)
ax.plot(N, LInfErr_B1_SSP64, color = RWTH_Blue_RGB[0], linestyle='dashed')

ax.scatter(N, LInfErr_B1_SSP104, label = r'$2N^*(10,4)$', color = RWTH_Green_RGB)
ax.plot(N, LInfErr_B1_SSP104, color = RWTH_Green_RGB[0], linestyle='dashed')

ax.scatter(N, LInfErr_B1_SSP144, label = r'$2N^*(14,4)$', color = RWTH_Orange_RGB)
ax.plot(N, LInfErr_B1_SSP144, color = RWTH_Orange_RGB[0], linestyle='dashed')


#ax.loglog(N, np.multiply(3e-1, np.power(np.array(N, dtype=float), -2) ), linestyle='dotted',
#          label = r'$\mathcal{O}\left(N^{-2}\right)$',
#          color = 'black') # Order two line fitted

ax.loglog(N, np.multiply(8e-2, np.power(np.array(N, dtype=float), -4) ), linestyle='dashdot',
          label = r'$\mathcal{O}\left(N^{-4}\right)$',
          color = 'black') # Order two line fitted

ax.loglog(N, np.multiply(6e-2, np.power(np.array(N, dtype=float), -5) ), linestyle='dotted',
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
plt.title(r"$L^\infty_{error}$ (Alfven Wave)")

# Scale while preserving aspect ratio
width, height = fig.get_size_inches()
factor = InchesX / width
fig.set_size_inches(width * factor, height * factor)

plt.tight_layout()  # Title, labels, ... to screen
plt.savefig('LInf_error.pgf', dpi=500, bbox_inches = 'tight', pad_inches = 0)
plt.show()
