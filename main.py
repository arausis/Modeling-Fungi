import os
import sys
from scipy.integrate import odeint
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')


# A class defining a fungal strain
class fungus:
    r = None
    h = None
    m = None
    N = 0
    def __init__(self, R, H, M, n):
        self.r = R
        self.h = H
        self.m = M
        self.N = n

    def describe(self):
        print("r : " + str(self.r))
        print("h : " + str(self.h))
        print("m : " + str(self.m))
        print("N : " + str(self.N))

# A class defining the environment that contains the fungi
class environment:
    beta = None
    alpha = None
    h = None
    k = None

    def __init__(self, b,a,H,K):
        self.beta = b
        self.alpha = a
        self.h = H
        self.k = K

    def describe(self):
        print("beta : " + str(self.beta))
        print("alpha : " + str(self.alpha))
        print("h : " + str(self.h))
        print("k : " + str(self.k))


# This function takes the path to our variables folder, and returns a touple of our environment class and a list of all of our fungal strains with variables filled in
def getVars(PATH):
    try:
        files = os.listdir(PATH)
    except:
        print("Cannot access variables folder!")
        sys.exit()

    # The environment
    env = None

    # A list of all fungal strains
    fungi = []


    for file in files:
        if "variables" in file:
            with open(PATH + "/" + file, "r") as evars:
                varlist = evars.readlines()
                b = None
                a = None
                H = None
                K = None
                for v in varlist:
                    var = v.split(":")
                    if var[0] == "beta":
                        b = float(var[1].strip())
                    elif var[0] == "alpha":
                        a = float(var[1].strip())
                    elif var[0] == "h":
                        H = float(var[1].strip())
                    elif var[0] == "k":
                        K = float(var[1].strip())
                    else:
                        raise Exception("Unexpected variable in variables.txt") 

                env = environment(b, a, H, K)
        else:
            with open(PATH + "/" + file, "r") as fvars:
                varlist = fvars.readlines()
                R = None
                H = None
                M = None
                n = 0
                for v in varlist:
                    var = v.split(":")
                    if var[0] == "r":
                        R = float(var[1].strip())
                    elif var[0] == "h":
                        H = float(var[1].strip())
                    elif var[0] == "m":
                        M = float(var[1].strip()) 
                    elif var[0] == "N":
                        n = float(var[1].strip())
                    else:
                        raise Exception("Unexpected variable in " + file) 
                fungi.append(fungus(R,H,M,n))
    return(env, fungi) 

# A simpler fungi class that has pre-calculated an "r-prime"
'''
r-prime: the entire constant that we can multiply to our basic carrying capacity function 

'''
class runfungi:
    n0 = 0
    rprime = 0

    def __init__(self, f, e):
        self.n0 = f.N
        self.rprime = f.r * (1 - (e.beta * abs(e.h - f.h) / f.m) )

    def describe(self):
        print("n0 : " + str(self.n0))
        print("r-prime: " + str(self.rprime))

class model:
    fungi = []
    deltaT = 0.1 
    k = 0
    x0 = []

    def __init__(self, e, f):
        self.k = e.k
        for i in f:
            self.x0.append(i.N)
            self.fungi.append(runfungi(i, e))

    def calculate(self, x, t):
        dndt_list = []
        N = 0

        for i in x:
            N += i

        for j in range(len(self.fungi)):
            dndt = self.fungi[j].rprime * x[j] * (1- (N / self.k) )
            dndt_list.append(dndt)


        return dndt_list

'''
------------------------------------------------------------------------
The main body of the code:

'''


'''
Some things that we might want to change
'''
# The directory of your variable folder
PATH = "./vars" 
# How long to run the simulation for (maximum t-value)
tmax = 5
# How many iterations of our t-value we want (how smooth the graph is)
tstep = 1000


'''
Initializing the model
'''

# Gets all of the variables from external folders
e, f = getVars(PATH)

# Displays Information about the Environment
print("\n\n Environment: ")
e.describe()

# Displays information about each Fungus
print("\n\n Fungi: ")
for i in f:
    i.describe()
    print("\n")

# Pre-calculates an "r-prime" value (see more above), and then displays a list of r-prime values
print("\n\n Runtime Fungi: ")
mod = model(e, f)

for i in mod.fungi:
    i.describe()
    print("\n")


'''
Running the model
'''
# Runs the ODE
t = np.linspace(0,15,1000)
x = odeint(mod.calculate, mod.x0 ,t)
 

#Plots the Data (note, this plots the population of each fungal strain, NOT THE FUNGAL DECOMPOSITION RATE)
for i in range(len(x[0])):
    plt.semilogy(t, x[:,i])

plt.show(block=True)
