import os
from helper import getVars
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
            dndt = self.fungi[j].r * (1 - abs(self.fungi[j].h - self.h)/self.fungi[j].m - (N / self.k)) * x[j]          
            dndt_list.append(dndt)



        return dndt_list

if __name__ == "__main__":
    # Initializing the model
    e, f = getVars(PATH)

    # Displays Information about the Environment
    for environ in e:
        environ.describe()
        print("\n")

    # Displays information about each Fungus
    print("\n\n Fungi: ")
    for i in f:
        i.describe()
        print("\n")


    # Create some model objects
    mod = []
    for environ in e:
        temp = model(environ, f)
        mod.append(temp)

    # Running our models
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
