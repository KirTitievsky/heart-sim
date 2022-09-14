

def get_config():
    c = dict({
        "dt": 0.1,
        "dx": 0.1,
        "L": 10,
        "T": 200,
        "k_activation": 10,
        "k_deactivation" : 3,
        "k_recovery": 0.1
       }
    )
    return c

def init_state(c, num_cells):

    state = {}
    
    state["active"] =     [0 for i in range(0,num_cells)]
    state["resting"] = [0 for i in range(0, num_cells)]
    state["ready"] =      [1 for i in range(0, num_cells)]

    return state




def set_ics(c, state):
    """Start with a spike of activated cells in the middle"""

    i = len(state["active"])//2
    state["active"][i] = 0.5
    state["ready"][i] = 0.5

W = [0.25, 0.5, 0.25]


from vector_ops import local_average

def derivative(c, state):
    N = len(state["ready"])

    def zeros(N):
        return list(0 for i in range(N))
   
    dx = {} 
    for k in state.keys():
       dx[k] = zeros(N) 

    
    for i in range(N):
        
        A = state["active"]
        a = A[i-1]*W[0] + A[i]*W[1] + A[i+1 - ((i + 1)//N) * N]*W[2]

        activation = c["k_activation"] * state["ready"][i] * a
        deactivation = c["k_deactivation"] * state["active"][i]
        recovery = c["k_recovery"] * state["resting"][i]

        dx["active"][i] =  activation - deactivation
        dx["resting"][i] = deactivation - recovery
        dx["ready"][i]   = recovery - activation

    return dx
 

c = get_config()

L = c["L"]
dx = c["dx"]


num_cells = round(L/dx)
dx_actual = 1.0*L/num_cells


x = init_state(c, num_cells)
set_ics(c, x)




from matplotlib import pyplot as plt, colors
from numpy import log10
import numpy as np




keys = ["active", "resting", "ready"]
xhistory = []
fig, axes = plt.subplots(len(keys) + 2,1)
for i in range(len(keys)): axes[i].set_title(keys[i])

num_steps = int(c["T"]/c["dt"])

for t in range(num_steps):
    dx = derivative(c, x)

    color_base = 0.5
    r = log10(t+1) / log10(num_steps) / (color_base + 1)
    color = (r,r,r)

    if t % (num_steps // 40) == 0:
        print("Time %i %f"%(t, t*c["dt"]))     
        xhistory.append(x["active"].copy())

        for k in keys:
            axes[keys.index(k)].plot(x[k], label = "%.1f"%(t*c["dt"]), color = color, linewidth=1 )
        
    for k in keys:
        for j in range(len(x[k])):
            x[k][j] += c["dt"] * dx[k][j]* c["dt"]

xhistory = np.array([np.array(i) for i in xhistory])

## show what happens at the midpoint as a summary of the trend
axes[3].set_title("Active cell fraction at point of initial stimulus.")
axes[3].plot(xhistory[:,xhistory.shape[1]//2])
axes[4].set_title("Active cell fraction over time (bottom is 0), Log transformed")
axes[4].pcolormesh(xhistory, cmap = 'plasma', norm=colors.LogNorm(vmin=0.01, vmax=100))
fig.subplots_adjust(hspace = 1)
plt.show()

