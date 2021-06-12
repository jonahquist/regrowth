import numpy as np
from processing_py import *
import random
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.integrate import RK45
import math

# dy/dt = f(t, y)
# y(0) = y_0

# all of the following are over pi
A_MAX = 206.0 #kg * m^-3 * y^-1  maximal assimilation
C = 73.0 # kg * m^-3 * y^-1 respiration and litter coefficient
C_2 = 73.0 # kg * m^-3 * y^-1 respiration coefficient for storage

G_1 = 365.0 # kg * m^-3 density of carbon
H_1 = 67.5 # m limit height of occurrence for photsynthesis * Beta
H_2 = 47.0 # m limit height of occurrence of meristem-sustained growth
S_1 = 0.2 # critical value of storage unit per mass
Q_R = 1000.0 # W * m^-2 characteristic PAR dependency
ALPHA_2 = 0.34 # ratio of meristematic activity
C_R = 500.0 # ppm characteristic CO2 mixing ratio dependency
T_OPT = 18.0 # deg C optimal temp for carbon assimilation
T_1 = 21

#dumb vars
S = 1
C_G = 409.8 # ppm
Q_1 = 1200 # W
T = 18
 
 #* ( (H_1 - 0.75*h) / (H_1) )


def model(t, y):
    dydt = ((A_MAX * (1 - (math.e)**(-C_G / C_R)) * (1 - (math.e)**(-Q_1 / Q_R)) * (1 - ( (T - T_OPT ) / T_1)**2)) * y - (C_2 * S * y)) / ( ((S - S_1) * G_1) / (1 + S_1))
    return dydt

# initial condition
n = 100 # number of years
t_span = np.array([0, n])
times = np.linspace(t_span[0], t_span[1], n+1)
y0 = np.array([5.0])

#solve the IVP
soln = solve_ivp(model, t_span, y0, t_eval = times)
def toHeight(t, y):
    a = 1
t = soln.t # times at which it was evaluated
Y = soln.y[0] # values at the times
print (Y[10])

# time points

class TreeGen:

    def __init__(self, max_depth):
        self.root = TreeBranch(0.707, 30, 90)
        self.max_depth = max_depth
        self.coordList = []

    # post: performs one step of tree growth; phased process
    #       called once a year
    def step(self, l):
        self.root.growCycle(l, 1)

    def treeCoords(self):
        self.coordList= []
        self.treeLoop(self.root, [0.0, 0.0])

    def treeLoop(self, c, prevCoords):
        if c:
            coords = c.getCoords(prevCoords)
            self.coordList.append(coords)
            self.treeLoop(c.left_child, coords[2:])
            self.treeLoop(c.right_child, coords[2:])
    

class TreeBranch:
    

    def __init__(self, r, a, va, b=None):
        self.ratio = r
        self.branching_angle = a
        self.parent = b
        self.left_child = None
        self.right_child = None
        self.length = 0.0
        self.radius = 0.0
        self.vector = np.array([0.0, 1.0])
        self.vector_angle = va # vertical
        self.terminal = True

    def growCycle(self, len, rad):
        
        self.length += (len/5) * random.random()
        self.radius += rad/10
        if (self.terminal):
            self.terminal = False
            # self.calcAngles()
            angle = self.vector_angle + (random.random() * self.branching_angle + 15)
            self.left_child = TreeBranch(self.ratio, self.branching_angle, angle, self)
            self.left_child.vector = np.array([np.cos(angle * np.pi/180), np.sin(angle * np.pi/180)])

            angle = self.vector_angle - (random.random() * self.branching_angle + 15)
            self.right_child = TreeBranch(self.ratio, self.branching_angle, angle, self)
            self.right_child.vector = np.array([np.cos(angle * np.pi/180), np.sin(angle * np.pi/180)])
            
        else:
            self.left_child.growCycle(len, rad * self.ratio)
            self.right_child.growCycle(len, rad * self.ratio)
            
    def calcAngles(self):
        return

    def getCoords(self, prevCoords):
        arr = [0.0, 0.0, 0.0, 0.0, 0.0]
        arr[0] = prevCoords[0]
        arr[1] = prevCoords[1]
        arr[2] = (self.vector[0] * self.length) + arr[0]
        arr[3] = (self.vector[1] * self.length) + arr[1]
        arr[4] = self.radius
        return arr

app = App(800, 800)
myTreeGen = TreeGen(0)

#
#for i in range(10):
#    print("Step " + str(i))
#    myTreeGen.step()
#    myTreeGen.treeCoords()
#    if (i % 4 == 0):    
#        app.stroke(255, 0 ,0)
#    elif (i % 3 == 0):
#        app.stroke(0, 255 ,120)
#    else:
#        app.stroke(0, 0 ,255)
#    for i in myTreeGen.coordList:
#        app.line((i[0] * 50) + 400, -1 * (i[1] * 50) + 800, (i[2] * 50) + 400, -1 * (i[3] * 50) + 800)
#        app.redraw()
#

img = app.loadImage("C:\\Users\\jonah\\Projects\\Regrowth\\leaf.png")

for i in range(8):
    myTreeGen.step(Y[i])
myTreeGen.treeCoords()

SCALE = 25

for i in myTreeGen.coordList:
    app.strokeWeight(i[4] * SCALE)
    app.line((i[0] * SCALE) + 400, -1 * (i[1] * SCALE) + 800, (i[2] * SCALE) + 400, -1 * (i[3] * SCALE) + 800)
    
    app.redraw()

app.stroke(255, 0 ,0)

plt.scatter(t, Y)
plt.xlabel("time")
plt.ylabel("y(t)")
plt.show()
