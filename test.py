import matplotlib.pyplot as plt
import numpy as np
import re
import csv
import cmath
#m1 [[-2, -1, 2],[2, 1, 4],[-3, 3, -1]]
"""
desirede1 = complex(-2.79402,2.69138)
desirede2 = complex(-2.79402,-2.69138)
desirede3 = complex(3.58803)
"""
#m2 [[-1, 2, -3],[4, 5, 6],[-7, 8, -9]]
"""
desirede1 = complex(-14.0242,0)
desirede2 = complex(7.94753,0)
desirede3 = complex(1.07664,0)
"""
#m3 [[1, -5, -5],[3, 1, 2],[1, 7, 4]]

desirede1 = complex(4.14404,4.21791)
desirede2 = complex(4.14404,-4.21791)
desirede3 = complex(-2.28808,0)


iter = []
e1 = []
e2 = []
e3 = []
lines = []

with open("datafiles/m3.txt") as f:
    for l in f:
        mo = re.match(r'^([0-9]+),(-?[0-9]+\.?[0-9]+e?[+-]?[0-9]+[+-]?[0-9]*\.?[0-9]+e?[+-]?[0-9]*j?),(-?[0-9]+\.?[0-9]+e?[+-]?[0-9]+[+-]?[0-9]*\.?[0-9]+e?[+-]?[0-9]*j?),(-?[0-9]+\.?[0-9]+e?[+-]?[0-9]*[+-]?[0-9]*\.?[0-9]+e?[+-]?[0-9]*j?)?$', l)
        if mo != None:
            iter.append(int(mo.group(1)))
            e1.append(complex(mo.group(2)))
            e2.append(complex(mo.group(3)))
            e3.append(complex(mo.group(4)))
        else:
            print(l)
        lines.append(l)
        

print(len(e1))
print(len(lines))

iter = np.array(iter)

e1 = np.array(e1)
e2 = np.array(e2)
e3 = np.array(e3)

e1dif = np.abs((np.absolute(desirede1)-np.absolute(e1))/np.absolute(e1))
e2dif = np.abs((np.absolute(desirede2)-np.absolute(e2))/np.absolute(e2))
e3dif = np.abs((np.absolute(desirede3)-np.absolute(e3))/np.absolute(e3))

fig, ax = plt.subplots()
ax.scatter(iter, e1dif, label="First Eigenvalue")
ax.scatter(iter, e2dif, label="Second Eigenvalue")
ax.scatter(iter, e3dif, label="Third Eigenvalue")
#ax.scatter(iter, (e1dif+e2dif)/2, label="AVG(e1,e2)")
#ax.scatter(iter, e3dif, label="e3")
ax.legend(loc="best")
ax.grid()
ax.set_title("Convergence of eigenvalue solver")
ax.set_xscale('log')
ax.set_xlabel("Iterations")
ax.set_ylabel("Normalized difference")
ax.set_ylim(0, 1)

plt.show()