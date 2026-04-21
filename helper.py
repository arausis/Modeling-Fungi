import sys
import os

# A class defining a fungal strain
class fungus:
    r = None
    h = None
    m = None
    N = 0 
    S = False
    def __init__(self, R, H, M, n, s): 
        self.r = R # Hyphal extension rate
        self.h = H # Climate Value (ideal moisture)
        self.m = M # Climate Tolerance (moisture tolerance)
        self.N = n # Population
        self.S = s # Sapotrophic y/n

    def describe(self):
        print("r : " + str(self.r))
        print("h : " + str(self.h))
        print("m : " + str(self.m))
        print("N : " + str(self.N))
        print("Sapotrophic?: " + str(self.S))

# A class defining the environment that contains the fungi
class environment:
    h = None
    k = None
    name = None

    def __init__(self,H,K, n):
        self.h = H # Climate Value (moisture level)
        self.k = K # Carrying Capacity (in fungal biomass)
        self.name = n 

    def describe(self):
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
                    else:
                        raise Exception("Unexpected variable in " + file) 
                fungi.append(fungus(R,H,M,n,s))
            
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
                if var[0] == "k":
                    K = float(var[1].strip())
                elif "h" in var[0] :
                    H = float(var[1].strip())
                elif var[0] == "name":
                    name = str(var[1].strip())
                else:
                    raise Exception("Unexpected variable in variables.txt") 
            env.append(environment(H, K, name))
    return(env, fungi) 
