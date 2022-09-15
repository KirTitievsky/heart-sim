"""This class is initialized with the output of Sim and
helps quantify and visualize it."""
import matplotlib as mpl
import matplotlib.pyplot as p
import numpy as n


class Vis():
    def __init__(self, sim):
        xs =  n.linspace(0, sim.L-1, int(sim.L))
        ts = n.linspace(0,sim.t, len(sim.xhistory))
        X, T = n.meshgrid(xs, ts)

        self.A = n.array([n.array(row) for row in sim.xhistory])
        fig, ax = p.subplots(2,1, figsize=[5,10])
        ax[0].set_title("Fraction of activated cells over time.")
        ax[0].pcolormesh(T.transpose(),X.transpose(), self.A.transpose()+0.001\
                , cmap='plasma', norm=mpl.colors.LogNorm(vmin=0.01, vmax=100))

        ax[1].set_title("Activation at 0 and L/2")
        # find the location of the peak 
        peak       = n.argmax(self.A[0,:])
        anti_peak  = peak - int(sim.L)//2

        ax[1].plot(ts, self.A[:,peak],label = 'max')
        ax[1].plot(ts, self.A[:,anti_peak],label = 'min')
        ax[1].set_xlabel('time')

        self.fig = fig
        self.ax = ax

