import matplotlib.pyplot as plt
import math as m
import numpy as np
import re
import glob

DOUBLEcholesky = glob.glob('DOUBLE/cholesky/IPC/*.txt')
FLOATcholesky = glob.glob('FLOAT/cholesky/*/*.txt')

DOUBLElu = glob.glob('DOUBLE/lu/*/*.txt')
FLOATlu = glob.glob('FLOAT/lu/*/*.txt')

DOUBLEmatmul = glob.glob('DOUBLE/matmul/*/*.txt')
FLOATmatmul = glob.glob('FLOAT/matmul/*/*.txt')

time = []
size = []
count1 = []
test1 = []
count2 = []
test2 = []

for n in DOUBLEcholesky:
    with open(n, 'r') as f:
        m=re.match(r"Total time: ([0-9].[0-9]+)\n\n[A-z ]+'\.\/[A-z ]*([0-9]+)[0-9\- ]*':\n\n[ ]+([0-9+,]+)[ ]+([A-z+\-_.1]+)[:uP]*[A-z \-#.% 0-9]*\n[ ]+([0-9+,]+)[ ]+([A-z+\-_.1]+)[:uP]*", f.read())
        print(f.name)
        time.append(m.group(1))
        size.append(m.group(2))
        count1.append(m.group(3))
        test1.append(m.group(4))
        count2.append(m.group(5))
        test2.append(m.group(6))

time, size, count1, test1, count2, test2 = zip(*sorted(zip(time, size, count1, test1, count2, test2)))

time = [eval(i) for i in time]
size = [eval(i) for i in size]
count1 = [eval(f.replace(',','')) for f in count1]
count2 = [eval(f.replace(',','')) for f in count2]

r = [i / j for i, j in zip(count1, count2)]

fig1, ax1 = plt.subplots()
ax1.plot(size,time)

fig2, ax2 = plt.subplots()
ax2.plot(size,count1)

fig3, ax3 = plt.subplots()
ax3.plot(size,count2)

fig4, ax4 = plt.subplots()
ax4.plot(size,r)

plt.show()