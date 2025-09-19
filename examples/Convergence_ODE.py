#!/usr/bin/python3

import matplotlib.pyplot as plt
import matplotlib
import numpy as np

# Plot Absolute Errors of Bernoulli Equation (ODE)

dts = [0.1, 0.05, 0.025, 0.0125, 0.00625]

# SSP64_2Nstar
Err_SSP64 = [5.2902453107073200e-07, 2.7614879938298031e-08, 1.5391972141287624e-09, 9.0146112796674061e-11, 5.4229953860840396e-12]

# SSP104_2Nstar
Err_SSP104 = [1.9422875765506831e-06, 1.2361337575761411e-07, 7.7882689186736798e-09, 4.8858073142810099e-10, 3.0592195443546188e-11]

# SSP144_2Nstar
Err_SSP144 = [1.3332019179745913e-06, 8.2978697912849952e-08, 5.1708015647022876e-09, 3.2264568794460047e-10, 2.0198731576215323e-11]


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

### ACTUAL PLOTTING: Absolute errors ###

ax.scatter(dts, Err_SSP64, label = r'2N$^*$(6,4)', color = RWTH_Blue_RGB)
ax.plot(dts, Err_SSP64, color = RWTH_Blue_RGB[0], linestyle='dashed')

ax.scatter(dts, Err_SSP104, label = r'2N$^*$(10,4)', color = RWTH_Green_RGB)
ax.plot(dts, Err_SSP104, color = RWTH_Green_RGB[0], linestyle='dashed')

ax.scatter(dts, Err_SSP144, label = r'2N$^*$(14,4)', color = RWTH_Orange_RGB)
ax.plot(dts, Err_SSP144, color = RWTH_Orange_RGB[0], linestyle='dashed')

ax.loglog(dts, np.multiply(8e-3, np.power(np.array(dts, dtype=float), 4) ), linestyle='dashdot',
          label = r'$\mathcal{O}\left(\Delta t^{\,4}\right)$',
          color = 'black') # Order two line fitted


# Turn on logscale (no native support for logarithmic scatter)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel(r'$\Delta t$')

### LIMITS SECTION ###
# eps_x = 0.1
# ax.set_xlim([1/128 * (1 - eps_x), 1 * (1 + eps_x)])

### GRID SECTION ###
ax.grid(axis ='both', which='major', alpha=0.1, linewidth = 1.5, color ='black')
ax.set_axisbelow(True)  # Hide grid behind bars

### LEGEND SECTION ###
ax.legend(loc = "lower right")

### TICKS SECTION ###

ax.set_xticks(dts)

ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

ax.get_xaxis().set_tick_params(which='minor', size=0)
ax.get_xaxis().set_tick_params(which='minor', width=0) 

ax.set_xticklabels([r"$1/10$", r"$1/20$", r"$1/40$", r"$1/80$", r"$1/160$"])

### LIMITS SECTION ###
# eps_x = 0.1
# ax.set_xlim([4 * (1 - eps_x), 64 * (1 + eps_x)])

# Make bounding lines thicker
'''
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2)
'''

### TITLE SECTION ###
plt.title(r"Absolute Errors (Bernoulli Equation)")

# Scale while preserving aspect ratio
width, height = fig.get_size_inches()
factor = InchesX / width
fig.set_size_inches(width * factor, height * factor)

plt.tight_layout()  # Title, labels, ... to screen
plt.savefig('error_ode.png', dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.show()
