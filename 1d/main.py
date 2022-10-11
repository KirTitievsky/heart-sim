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
        self.k_activ = list(np.logspace(-1,1,5))
        self.k_deact = list(np.logspace(-1,1,5))
        self.k_recov = list(np.logspace(-1,1,5))

        self.L = 100
        self.x_sample = 0
        self.results = []
        self.num_steps = 1000
        self.period = 0


    def go(self):
        for k_activ in self.k_activ:
            for k_deact in self.k_deact:
                for k_recov in self.k_recov: 
                    sim = Sim(k_activ = k_activ, k_deact = k_deact, k_recov = k_recov, L = 100)
                    sim.debug = False

                    sim.run(num_steps = self.num_steps)
                
                    # find where in x the initial spike in A was
                    A = sim.xhistory
                    self.results.append({
                        "k_activ": k_activ, 
                        "k_deact": k_deact,
                        "k_recov": k_recov,
                        "x": self.x_sample,
                        "A": [A[t][self.x_sample] for t in range(len(sim.xhistory))],
                        "t": sim.thistory
                    })

                    #TODO: Consider making this an generator so you can plot incremental results

            
class Sim():

    def __init__(self, k_activ = 1, k_deact = 1, k_recov = 1, L = 100):
        """ The user of the class should explicitly set paratemeters after initializing."""
        # domain -- length and time
        self.L  = L  # assume dx = 1 
        # rate constants
        self.k_activ = k_activ  # activation:  Inactive -> Active 
        self.k_deact = k_deact  # deactivation:  Active -> Spent state
        self.k_recov = k_recov  # recovery:   Spent -> ready

        self.debug = True # print all logs by default

        # state will be saved every this many seconds or every step if this is smaller
        self.sampling_period_steps = 10

        self.states = ["active", "spent", "ready"]

        self.state = {}

        # save the history of the active concentration here
        self.xhistory = []
        self.shistory = []
        self.thistory = []


        self.t = 0


    def pulse(self, A_peak = 1.0):
        """Creates an pulse of activated cells at x = L/2. Sets A(t, x=0) = A_peak if A(t,x=0) < A_peak)"""
        # the initial condition is a spike of active concentration at L/2
        i = self.L//2
        A,S,R = [self.state[k][i] for k in ("active", "spent", "ready")]

        if A < A_peak:
            # we activate only the ready cells
            self.state["active"][i] = min(A + R, A_peak)
            self.state["ready"][i] = 1 - S - self.state["active"][i] 

    def prepare(self):

        for k in self.states: 
            self.state[k] = [0 for i in range(0,self.L)]
        self.state["ready"] = [1 for i in range(0,self.L)]
        
        self.pulse()


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
            activation[i] =   self.k_activ * self.state["ready"][i] * Ad[i]
            deactivation[i] = self.k_deact * self.state["active"][i]
            recovery[i] =     self.k_recov * self.state["spent"][i]

         
            dx["active"][i]  = activation[i] - deactivation[i]
            dx["spent"][i]   = deactivation[i] - recovery[i]
            dx["ready"][i]   = recovery[i] - activation[i]

        return dx
 

    def log(self, current_step):

        if current_step % self.sampling_period_steps == 0:
            self.xhistory.append(self.state["active"].copy())
            self.shistory.append(self.state["spent"].copy())
            self.thistory.append(self.t)

            if self.debug:
                print("At step %d"%(current_step) )

    def run(self, num_steps = 0, dt = 0.1):
        self.prepare()

        self.log(0)

        for t in range(num_steps):
            if t % int(self.period/dt) == 0:
                self.pulse()

            x = self.state
            dx = self.d()
            self.t += dt
            for k in self.states:
                for j in range(len(x[k])):
                    x[k][j] += dx[k][j] * dt 
                    #TODO: control for numerical errors -- sum of all x[k] == 1


            self.log(t)



