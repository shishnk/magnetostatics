import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from decimal import Decimal as dcm
from matplotlib.widgets import Button

elements = []
points = []
xL, yL = [], []
q = []

fig, ax = plt.subplots(1)

file1 = "points.txt"
path1 = os.path.join('C:\\Users\\Hukutka\\source\\repos\\Magnetostatics\\bin\\Debug\\net7.0\\output', file1)
file2 = "elements.txt"
path2 = os.path.join('C:\\Users\\Hukutka\\source\\repos\\Magnetostatics\\bin\\Debug\\net7.0\\output', file2)
file3 = "q.txt"
path3 = os.path.join('C:\\Users\\Hukutka\\source\\repos\\Magnetostatics\\bin\\Debug\\net7.0\\output', file3)

with open(path1) as file:
    for line in file:
        xUniq, yUniq = line.split()
        points.append([float(xUniq), float(yUniq)])
        xL.append(float(xUniq))
        yL.append(float(yUniq))

with open(path2) as file:
    for line in file:
        vert1, vert2, vert3, vert4, area = map(int, line.split())
        elements.append([vert1, vert2, vert3, vert4, area])

nx = elements[0][2] - 1
ny = len(elements) // nx
xUniq = []
yUniq = []

for i in range(nx + 1):
    xUniq.append(points[i][0])

for i in range(ny + 1):
    yUniq.append(points[i * (nx + 1)][1])

q = np.zeros([len(points)])

idx = 0

with open(path3) as file:
    for line in file:
        q[idx] = dcm(line)
        idx += 1

X, Y = np.meshgrid(xUniq, yUniq)
q = np.reshape(q, (len(yUniq), len(xUniq)))
colorBar = plt.contourf(X, Y, q, levels=100, cmap='jet')

class Index:
    def prev(self, event):
        for elem in elements:
            if elem[4] == 0:
                color = 'lightblue'
            elif elem[4] == 1:
                color = 'blue'
            elif elem[4] == 2:
                color = 'red'
            elif elem[4] == 3:
                color = 'lime'

            quadrangle = patches.Polygon([
                [xL[elem[0]], yL[elem[0]]],
                [xL[elem[1]], yL[elem[1]]],
                [xL[elem[3]], yL[elem[3]]],
                [xL[elem[2]], yL[elem[2]]]
            ],
                edgecolor='black', facecolor=color, linewidth=1)
            ax.add_patch(quadrangle)

        plt.plot(xL, yL, " ")
        plt.draw()

    def next(self, event):
        plt.draw()


callback = Index()
axprev = fig.add_axes([0.81, 0.01, 0.1, 0.05])
bprev = Button(axprev, 'Mesh')
bprev.on_clicked(callback.prev)

plt.colorbar(colorBar, ax=ax, format='%.0e')
plt.show()
