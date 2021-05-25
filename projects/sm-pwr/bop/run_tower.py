#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment
# https://cortix.org
"""Cortix Run File"""

import matplotlib.pyplot as plt

import unit

from cortix import Cortix
from cortix import Network

from cooling_tower import CoolingTower

def main():

    # Debugging
    make_plots = True
    make_run   = True

    # Preamble
    end_time = 10*unit.minute
    time_step = 1.5*unit.second
    show_time = (True, 5*unit.minute)

    plant = Cortix(use_mpi=False, splash=True) # System top level

    plant_net = plant.network = Network() # Network

    # Tower

    tower = CoolingTower()  # Create reactor module

    # Steady state conditions for NuSCale case
    #primary_inflow_temp = (320.9+273.15)*unit.kelvin
    #secondary_inflow_temp = (149+273.15)*unit.kelvin
    #steamer = Steamer(primary_inflow_temp, secondary_inflow_temp)  # Create reactor module

    tower.name = 'Cooling-Tower'
    tower.save = True
    tower.time_step = time_step
    tower.end_time = end_time
    tower.show_time = show_time

    plant_net.module(tower)  # Add steamer module to network

    # Balance of Plant Network Connectivity

    plant_net.draw(engine='circo', node_shape='folder')

    # Run
    if make_run:
        plant.run()  # Run network dynamics simulation

    plant.close()  # Properly shutdow Cortix

    # Plots
    if make_plots and plant.use_multiprocessing or plant.rank == 0:

        # Steamer plots
        tower = plant_net.modules[0]

        (quant, time_unit) = tower.outflow_phase.get_quantity_history('temp')

        quant.plot(x_scaling=1/unit.minute, y_shift=273.15, x_label='Time [m]',
                   y_label=quant.latex_name+' [C]')
        plt.grid()
        plt.savefig('outflow-temp.png', dpi=300)

        (quant, time_unit) = tower.outflow_phase.get_quantity_history('flowrate')

        quant.plot(x_scaling=1/unit.minute, x_label='Time [m]',
                   y_label=quant.latex_name+' ['+quant.unit+']')
        plt.grid()
        plt.savefig('outflow-mass-flowrate.png', dpi=300)

        

if __name__ == '__main__':
    main()
