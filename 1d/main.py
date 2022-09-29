W = [0.25, 0.5, 0.25]

from matplotlib import pyplot as plt, colors
from numpy import log10
import numpy as np




from vector_ops import local_average



class Scan():
""" This class scans a simulation, Sim, through a parameter space of rate constants.
It outputs a set of histories A(0,t) and A(L/2,t).
"""
    def __init__(self):
        self.kRA = [1,10]
        self.kAS = [1,10]
        self.kSR = [1,10]
        self.L = 100
        self.x_sample = 0
        results = {}

    def go(self,):
        for kRA in self.kRA:
            for kAS in self.kAS:
                for kSR in self.kSR: 
                    sim = Sim(kRA = kRA, kAS = kAS, kSR = kSR, L = 100)
                    sim.debug = False

                    sim.run()
                
                    # find where in x the initial spike in A was
                    A = sim.xhistory
                    results.set({
                        "kRA": kRA, 
                        "kAS": kAS,
                        "kSR": kSR,
                        "x": self.x_sample,
                        "A": [A[t][self.x_sample] for t in range(len(sim.xhistory))]
                    })

                    #TODO: Consider making this an generator so you can plot incremental results
                    

class Sim():

    def __init__(self, kRA = 1, kAS = 1, kSR = 1, L = 100):
        """ The user of the class should explicitly set paratemeters after initializing."""
        # domain -- length and time
        self.L  = L  # assume dx = 1 
        # rate constants
        self.kRA = kRA  # activation:  Inactive -> Active 
        self.kAS = kAS  # deactivation:  Active -> Spent state
        self.kSR = kSR  # recovery:   Spent -> ready

        self.debug = True # print all logs by default

        # state will be saved every this many seconds or every step if this is smaller
        self.sampling_period_steps = 10

        self.states = ["active", "spent", "ready"]

        self.state = {}

        # save the history of the active concentration here
        self.xhistory = []


        self.t = 0

    def prepare(self):

        for k in self.states: 
            self.state[k] = [0 for i in range(0,self.L)]
        self.state["ready"] = [1 for i in range(0,self.L)]
        
        # the initial condition is a spike of active concentration at L/2
        i = self.L//2
        self.state["active"][i] = 0.5
        self.state["ready"][i] = 0.5



    def d(self):
        """Compute the reaction rates -- rate of change of the state."""
        #TODO: There is no reason we should not handle this the same way we do the state
        N = int(self.L)
        dx = {} 
        for k in self.states:
           dx[k] = list(0 for i in range(N))

      
        # compute a version of A where neighboring cells affect each other
        A = self.state["active"]
        Ad = [sum( A[(i+di)%N] * W[1+di] for di in (-1,0,1)) for i in range(N)]

      
        activation = [0]*N
        deactivation = [0]*N
        recovery = [0]*N

        for i in range(N):
            activation[i] =   self.kRA * self.state["ready"][i] * Ad[i]
            deactivation[i] = self.kAS * self.state["active"][i]
            recovery[i] =     self.kSR * self.state["spent"][i]

         
            dx["active"][i]  = activation[i] - deactivation[i]
            dx["spent"][i]   = deactivation[i] - recovery[i]
            dx["ready"][i]   = recovery[i] - activation[i]

        return dx
 

    def log(self, current_step):

        if current_step % self.sampling_period_steps == 0:
            self.xhistory.append(self.state["active"].copy())
            if self.debug:
                print("At step %d"%(current_step) )

    def run(self, num_steps = 0, dt = 0.1):
        self.prepare()

        self.log(0)

        for t in range(num_steps):
            x = self.state
            dx = self.d()
            self.t += dt
            for k in self.states:
                for j in range(len(x[k])):
                    x[k][j] += dx[k][j] * dt 
                    #TODO: control for numerical errors -- sum of all x[k] == 1

            self.log(t+1)

