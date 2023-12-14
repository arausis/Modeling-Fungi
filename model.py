import os
import matplotx
import sys
from scipy.integrate import odeint
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')
import threading


'''
Some things that we might want to change
'''
# The directory of your variable folder
PATH = "./vars" 
# How long to run the simulation for (maximum t-value)
tmax = 10
# How many iterations of our t-value we want (how smooth the graph is)
tstep = 1000
# The names of each environment, in order of their repsective h value in the variables file
envnames = ["arid", "semi-arid"]

# A class defining a fungal strain
class fungus:
    r = None
    h = None
    m = None
    N = 0 
    S = False
    c = 0
    def __init__(self, R, H, M, n, s, c): 
        self.r = R
        self.h = H
        self.m = M
        self.N = n
        self.S = s
        self.c = c

    def describe(self):
        print("r : " + str(self.r))
        print("h : " + str(self.h))
        print("m : " + str(self.m))
        print("N : " + str(self.N))
        print("Sapotrophic?: " + str(self.S))
        print("c : " + str(self.c))

# A class defining the environment that contains the fungi
class environment:
    beta = None
    alpha = None
    h = None
    k = None
    name = None

    def __init__(self, b,a,H,K, n):
        self.beta = b
        self.alpha = a
        self.h = H
        self.k = K
        self.name = n

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
        if "environments" not in file:
            with open(PATH + "/" + file, "r") as fvars:
                varlist = fvars.readlines()
                R = None
                H = None
                M = None
                n = 0
                s = False
                c = 0
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
                    elif var[0] == "S":
                        if str(var[1]).strip() == "T":
                            s = True
                        else:
                            s = False
                    elif var[0] == "c":
                        c = float(var[1].strip())
                    else:
                        raise Exception("Unexpected variable in " + file) 
                fungi.append(fungus(R,H,M,n,s,c))
            
    try:
        envfiles = os.listdir(str(PATH + "/environments"))
    except:
        print("Cannot access environments folder!")
        sys.exit()

    # The environment
    env = []

    for file in envfiles:
        with open(PATH + "/environments/" + file, "r") as evars:
            varlist = evars.readlines()
            b = None
            a = None
            H = None
            K = None
            name = None
            for v in varlist:
                var = v.split(":")
                if var[0] == "beta":
                    b = float(var[1].strip())
                elif var[0] == "alpha":
                    a = float(var[1].strip())
                elif var[0] == "k":
                    K = float(var[1].strip())
                elif "h" in var[0] :
                    H = float(var[1].strip())
                elif var[0] == "name":
                    name = str(var[1].strip())
                else:
                    raise Exception("Unexpected variable in variables.txt") 
            env.append(environment(b, a, H, K, name))
    return(env, fungi) 

# A simpler fungi class that has pre-calculated an "r-prime"
'''
r-prime: the entire constant that we can multiply to our basic carrying capacity function 

'''

class model:
    fungi = []
    deltaT = 0.1 
    k = 0
    x0 = []
    name = None
    h = None

    def __init__(self, e, f):
        self.k = e.k
        self.fungi = []
        self.x0 = []
        self.name = e.name
        self.h = e.h
        for i in f:
            self.x0.append(i.N)
            self.fungi = f

    def calculate(self, x, t):
        dndt_list = []
        N = 0
        for i in x:
            N += i

        for j in range(len(self.fungi)):
            dndt = self.fungi[j].r * (self.fungi[j].c - abs(self.fungi[j].h - self.h)/self.fungi[j].m - (N / self.k)) * x[j]          
            dndt_list.append(dndt)



        return dndt_list

'''
------------------------------------------------------------------------
The main body of the code:

'''

'''
Initializing the model
'''

# Gets all of the variables from external folders
e, f = getVars(PATH)

# Displays Information about the Environment
print("\n\n Environment: ")
for environ in e:
    environ.describe()
    print("\n")

# Displays information about each Fungus
print("\n\n Fungi: ")
for i in f:
    i.describe()
    print("\n")


mod = []
for environ in e:
    print(len(f))
    temp = model(environ, f)
    #print(len(temp.fungi))
    mod.append(temp)
    #print(mod)
    #print(len(mod[len(mod)-1].fungi))


'''
Running the model
'''

for j in range(len(mod)):

    # Runs the ODE
    t = np.linspace(0,tmax, tstep)
    x = odeint(mod[j].calculate, mod[j].x0 ,t)

    #Plots the Data (note, this plots the population of each fungal strain, NOT THE FUNGAL DECOMPOSITION RATE)
    x2 = np.zeros(len(x))
    for i in range(len(x[0])):
        plt.semilogy(t, x[:,i], label="Fungus " + str(i + 1))
        if mod[j].fungi[i].S:
            for k in range(len(x)):
                x2[k] += x[k][i] 

    #matplotx.line_labels()

    plt.figure(1)
    plt.title("Population of Fungal Species in " + mod[j].name)
    plt.xlim([0, tmax])

    #Plots Fugal decomposition rate
    plt.figure(2)
    plt.semilogy(t,x2)
    plt.title("Decomposition Rate in " + mod[j].name)
    plt.xlim([0, tmax])

    plt.show(block=True)
