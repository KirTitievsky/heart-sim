"""This class is initialized with the output of Sim and
helps quantify and visualize it."""
import matplotlib as mpl
import matplotlib.pyplot as p
import numpy as n

def simple(sim, logscale=True):
    """Simple visualiztion that plots A(t,x=L/2) for a single simulation in a 2-d graph."""
    x = [i[int(sim.L/2)] for i in sim.xhistory]
    p.plot(sim.thistory, x, label='Active'); 
    p.legend()
    if logscale:
        p.yscale('log')
    p.ylabel('A(t, L/2)')
    p.grid(True)
    p.xlabel('Time kA')
    p.title('k_act = %g, k_deact = %g, k_recov = %g, L = %g'%(sim.k_activ, sim.k_deact, sim.k_recov, sim.L))
    p.show()
   
def simple_ready(sim, ylim = 1):
    """Simple visualiztion that plots A(t,x=L/2) for a single simulation in a 2-d graph."""
    A = [i[int(sim.L/2)] for i in sim.xhistory]
    S = [i[int(sim.L/2)] for i in sim.shistory]
    R = [1- A[i] - S[i] for i in range(len(A))]


    p.plot(sim.thistory, A, color=(0.1, 0.1, 0.1), label='Active')
    p.plot(sim.thistory, R, color=(0.1,0.7,0.1), label='Ready')
    p.legend()
    p.ylabel('R(t, L/2)')
    p.grid(True)
    p.xlabel('Time kA')
    p.ylim(0, ylim)
    p.title('k_act = %g, k_deact = %g, k_recov = %g, L = %g'%(sim.k_activ, sim.k_deact, sim.k_recov, sim.L))
    p.show()

class Vis():
    """Produces a full visualization of the active cell field for a simulation."""
    def __init__(self, sim, lognorm = False):
        xs =  n.linspace(0, sim.L-1, int(sim.L))
        ts = n.linspace(0,sim.t, len(sim.xhistory))
        X, T = n.meshgrid(xs, ts)

        self.A = n.array([n.array(row) for row in sim.xhistory])
        fig, ax = p.subplots(2,1, figsize=[5,10])
        ax[0].set_title("Fraction of activated cells over time.")
        kwargs = {'cmap': 'plasma'}
        if lognorm:
            kwargs['norm'] = mpl.colors.LogNorm(vmin=0.01, vmax=100)
        ax[0].pcolormesh(T.transpose(),X.transpose(), self.A.transpose()+0.001, **kwargs)

        ax[1].set_title("Activation at 0 and L/2")
        # find the location of the peak 
        peak       = n.argmax(self.A[0,:])
        anti_peak  = peak - int(sim.L)//2

        ax[1].plot(ts, self.A[:,peak],label = 'max')
        ax[1].plot(ts, self.A[:,anti_peak],label = 'min')
        ax[1].set_xlabel('time')

        self.fig = fig
        self.ax = ax

from collections import namedtuple
class ScanVis():
    """Visualizes a Scan of simulations."""
    def __init__(self, scan):
        make_record = namedtuple("Dimension", ["id", "values", "display_name"])
        dims = [
            make_record("k_recov", scan.k_recov, "recov"),
            make_record("k_activ", scan.k_activ, "aciv"),
            make_record("k_deact", scan.k_deact, "deac"), 
        ]

        
        grid_dims = [len(dims[0].values), len(dims[1].values)]
        figsize = list(4*i for i in grid_dims)
        t_max = max(max(r['t']) for r in scan.results)
        subplot_kw = {"frame_on":False,'xmargin':0, 'ymargin':0, 'yticks': (0,0.5,1.0), 'xticks': (0,t_max)}
        fig, axes = \
            p.subplots(*grid_dims, figsize = figsize , sharex = True, sharey = True, subplot_kw = subplot_kw, squeeze=False)
        


        for row in range(grid_dims[0]):
            for col in range(grid_dims[1]):
                title =   "%s %.2g"%(dims[0].display_name, dims[0].values[row])
                title += " %s %.2g"%(dims[1].display_name, dims[1].values[col])

                axes[row,col].text(0.8,1, title, fontsize=10, fontfamily='serif', )

        
        # set titles 
        inner_parameter_values_string = ",".join("%.2g"%i for i in dims[-1].values)
        fig.suptitle("Activation over time for %s = %s. Brighter lines = higher values."%\
            (dims[-1].display_name,inner_parameter_values_string),\
            fontsize=14, fontweight=100
        )
        # loop through the results and plot them 
        for result in scan.results:
            row = dims[0].values.index(result[dims[0].id])
            col = dims[1].values.index(result[dims[1].id])
            label = "%s=%.1f"%(dims[2].id, result[dims[2].id])
           
            c =  result[dims[2].id]/max(dims[2].values)
            c *= 0.8
            line_color = (c,c,c)
            axes[row, col].plot(result['t'], result['A'], label = label, color = line_color, linewidth=1)
        
        self.axes = axes
