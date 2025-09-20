import numpy as np
from matplotlib import pyplot as plt

RealsScaled = [-15,10]
ImagsScaled = [-7,7]

# Plot Section

RWTH_Blue_RGB   = [(0, 84/256, 159/256)]
RWTH_Orange_RGB = [(246/256, 169/256, 0)]
RWTH_Green_RGB  = [(70/256, 171/256, 39/256)]
RWTH_Red_RGB    = [(204/256, 7/256, 30/256)]
RWTH_Petrol_RGB = [(0/256, 152/256, 161/256)]
RWTH_Yellow_RGB = [(225/256, 237/256, 0/256)]
RWTH_Purple_RGB = [(97/256, 33/256, 88/256)]

GoldenRatio = (1 + 5 ** 0.5) / 2

InchesX = 5
InchesY = InchesX / GoldenRatio
#fig, ax = plt.subplots(figsize = (InchesX, InchesY) )
fig, ax = plt.subplots()

ax.set_xlabel("Re")
ax.set_ylabel("Im")

ax.grid(axis ='y', which='both', alpha=0.1, linewidth = 1.5, color ='black')
ax.set_axisbelow(True)  # Hide grid
# Make zero line better visible
gridlines = ax.yaxis.get_gridlines()
gridlines[1].set_alpha(1)

ax.grid(axis ='x', which='both', alpha=0.1, linewidth = 1.5, color ='black')
gridlines = ax.xaxis.get_gridlines()
gridlines[-1].set_alpha(1)
#plt.show()

N = 200

xmin = 1.02 * min(RealsScaled)
xmax = 1
x = np.linspace(xmin, xmax, N)

ymin = 1.02 * min(ImagsScaled)
ymax = 1.02 * max(ImagsScaled)
y = np.linspace(ymin, ymax, N)

X, Y = np.meshgrid(x, y)

z = X + Y*1j

Degree = [14,10,6] # TODO

gamma = []
# 2N*(14,4)
gamma.append([1.0, 1.0, 0.5, 1/6, 1/24, 0.009705021572857, 0.002424593199150, 0.000578756059377, 0.000112566194836, 0.000016584859542, 0.000001788702831, 0.000000136126961, 0.000000006882578, 0.000000000206181, 0.000000000002753])
# 2N*(10,4)
gamma.append([1.0, 1.0, 0.5, 1/6, 1/24, 0.010417183051479, 0.002713139389985, 0.000552180171227, 0.000071657752512, 0.000005118309064, 0.000000150904305])
# 2N*(6,4)
gamma.append([1.0, 1.0, 0.5, 1/6, 1/24, 0.008448212703758, 0.000667275151153])
# RK4
#gamma = [1.0, 1.0, 0.5, 1/6, 1/24]


#GammaFileName = "./gamma_" + str(Degree) + ".txt"
#with open(GammaFileName, 'r') as GammaFile:
#  Lines = GammaFile.readlines()
 # for Line in Lines:
 #   Entries = Line.split()
 #   if len(Entries) > 0:  # Exclude blank lines
 #     gamma.append(float(Entries[0]))

# First fixed coefficients for 4th order accurate stability polynomial
# gamma.insert(0, 1.0/24.0)
# gamma.insert(0, 1.0/6.0)
# gamma.insert(0, 0.5)
# gamma.insert(0, 1)
# gamma.insert(0, 1)
colors = [RWTH_Orange_RGB, RWTH_Blue_RGB, RWTH_Red_RGB]
for j in range(0,len(gamma)):
  poly = gamma[j][Degree[j]] * z + gamma[j][Degree[j] - 1]

  for i in range(Degree[j] - 2, -1, -1):
    poly *= z
    poly += gamma[j][i]

  Z = abs(poly)
  plt.contourf(X, Y, Z, levels=[0,1], colors=colors[j])
  cs = plt.contour(X, Y, Z, levels=[0,1], linewidths = 0.5, colors="black")
  plt.plot([1, 2], [1, 2], color = colors[j][0], label = "2N*("+str(Degree[j])+",4)")
  
#cs = plt.contour(X, Y, Z, [1.0], colors = RWTH_Orange_RGB, linewidths = 2)
#plt.clabel(cs, inline=1, fontsize=10) # add labels to contours

# Dummy plot line for legend
#plt.plot([1, 2], [1, 2], color = RWTH_Orange_RGB[0], label = "Boundary of region of absolute stability")
#plt.plot([1, 2], [1, 2], color = RWTH_Orange_RGB[0], label = "E = 20")

plt.legend(loc = "lower left", bbox_to_anchor=(0.0, 0.075))

ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)

plt.title(r"Absolute domain of stability $|P_E(z) \leq 1 |$")

# Scale while preserving aspect ratio
width, height = fig.get_size_inches()
factor = InchesX / width
fig.set_size_inches(width * factor, 1.2*height * factor)

#plt.savefig('Boundary_StabilityRegion' + str(Degree) + '.pgf', bbox_inches = 'tight', pad_inches = 0)
plt.savefig('Boundary_StabilityRegions.pgf', bbox_inches = 'tight')

plt.show()
