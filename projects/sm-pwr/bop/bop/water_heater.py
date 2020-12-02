#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""Cortix module"""

import logging

import unit

import math
from scipy.integrate import odeint
import numpy as np

from iapws import IAPWS97 as steam_table

from cortix import Module
from cortix.support.phase_new import PhaseNew as Phase
from cortix import Quantity

class WaterHeater(Module):
    """Water heater system.

    Notes
    -----
    These are the `port` names available in this module to connect to respective
    modules: reactor, turbine.
    See instance attribute `port_names_expected`.

    """

    def __init__(self):
        """Constructor.

        Parameters
        ----------

        """

        super().__init__()

        self.port_names_expected = ['inflow', 'outflow', 'heat']

        # General attributes
        self.initial_time = 0.0*unit.second
        self.end_time = 1.0*unit.hour
        self.time_step = 10.0*unit.second

        self.show_time = (False, 10.0*unit.second)

        self.log = logging.getLogger('cortix')
        self.__logit = True # flag indicating when to log

        # Domain attributes

        # Configuration parameters

        self.volume = 15*unit.meter**3

        # Initialization

        self.inflow_temp = (20+273)*unit.kelvin
        self.inflow_pressure = 20*unit.bar
        self.inflow_mass_flowrate = 67*unit.kg/unit.second

        self.outflow_temp = (20+273.15)*unit.kelvin
        #self.outflow_temp_ss = 422 #k
        self.outflow_pressure = 34*unit.bar
        self.outflow_mass_flowrate = 67*unit.kg/unit.second

        self.heat_source_rate = 0*unit.watt

        # Outflow phase history
        quantities = list()

        temp = Quantity(name='temp',
                        formal_name='T', unit='K',
                        value=self.outflow_temp,
                        latex_name=r'$T$',
                        info='Water Heater Outflow Temperature')

        quantities.append(temp)

        press = Quantity(name='pressure',
                         formal_name='P', unit='Pa',
                         value=self.outflow_pressure,
                         latex_name=r'$P$',
                         info='Water Heater Outflow Pressure')

        quantities.append(press)

        self.outflow_phase = Phase(time_stamp=self.initial_time,
                                   time_unit='s', quantities=quantities)

    def run(self, *args):

        # Some logic for logging time stamps
        if self.initial_time + self.time_step > self.end_time:
            self.end_time = self.initial_time + self.time_step

        time = self.initial_time

        print_time = self.initial_time
        print_time_step = self.show_time[1]

        if print_time_step < self.time_step:
            print_time_step = self.time_step

        while time <= self.end_time:

            if self.show_time[0] and \
               (print_time <= time < print_time+print_time_step):

                msg = self.name+'::run():time[m]='+ str(round(time/unit.minute, 1))
                self.log.info(msg)

                self.__logit = True
                print_time += self.show_time[1]

            else:
                self.__logit = False

            # Communicate information
            #------------------------
            self.__call_ports(time)

            # Evolve one time step
            #---------------------

            time = self.__step(time)

    def __call_ports(self, time):

        # Interactions in the feed water outflow port
        #-----------------------------------------
        # one way "to" steam geneator

        # send to
        if self.get_port('outflow').connected_port:

            msg_time = self.recv('outflow')
            assert msg_time <= time

            temp = self.outflow_phase.get_value('temp', msg_time)
            pressure = self.outflow_phase.get_value('pressure', msg_time)
            outflow = dict()
            outflow['temperature'] = temp
            outflow['pressure'] = pressure
            outflow['mass_flowrate'] = self.outflow_mass_flowrate
            self.send((msg_time, outflow), 'outflow')

        # Interactions in the inflow port
        #----------------------------------------
        # one way "from" inflow

        # receive from
        if self.get_port('inflow').connected_port:

            self.send(time, 'inflow')

            (check_time, inflow) = self.recv('inflow')
            assert abs(check_time-time) <= 1e-6

            self.inflow_temp = inflow['temperature']
            self.inflow_pressure = inflow['pressure']
            self.inflow_mass_flowrate = inflow['mass_flowrate']

        # Interactions in the heat port
        #----------------------------------------
        # one way "from" heat

        # receive from
        if self.get_port('heat').connected_port:

            self.send(time, 'heat')

            (check_time, heat) = self.recv('heat')
            assert abs(check_time-time) <= 1e-6

            self.heat_source_rate = heat

    def __step(self, time=0.0):

        # Get state values
        u_0 = self.__get_state_vector(time)

        t_interval_sec = np.linspace(time, time+self.time_step, num=2)

        max_n_steps_per_time_step = 1000 # max number of nonlinear algebraic solver
                                         # iterations per time step

        (u_vec_hist, info_dict) = odeint(self.__f_vec, u_0, t_interval_sec,
                                         rtol=1e-4, atol=1e-8,
                                         mxstep=max_n_steps_per_time_step,
                                         full_output=True, tfirst=False)

        assert info_dict['message'] == 'Integration successful.', info_dict['message']

        u_vec = u_vec_hist[1, :]  # solution vector at final time step

        temp = u_vec[0] # primary outflow temp

        # Update state variables
        outflow = self.outflow_phase.get_row(time)

        time += self.time_step

        water = steam_table(T=temp, P=self.outflow_pressure/unit.mega/unit.pascal)
        assert water.phase != 'Vapour'

        self.outflow_phase.add_row(time, outflow)
        self.outflow_phase.set_value('temp', temp, time)
        self.outflow_phase.set_value('pressure', self.outflow_pressure, time)

        return time

    def __get_state_vector(self, time):
        """Return a numpy array of all unknowns ordered as shown.

        """

        u_vec = np.empty(0, dtype=np.float64)

        temp = self.outflow_phase.get_value('temp', time)
        u_vec = np.append(u_vec, temp)

        return  u_vec

    def __f_vec(self, u_vec, time):

        temp = u_vec[0]

        # initialize f_vec to zero
        f_tmp = np.zeros(1, dtype=np.float64) # vector for f_vec return

        #-----------------------
        # primary energy balance
        #-----------------------
        water = steam_table(T=temp,
                            P=self.inflow_pressure/unit.mega/unit.pascal)

        assert water.phase != 'Vapour'

        rho = water.rho
        cp = water.cp*unit.kj/unit.kg/unit.K
        vol = self.volume

        temp_in = self.inflow_temp

        tau = vol/(self.inflow_mass_flowrate/rho)

        #-----------------------
        # calculations
        #-----------------------
        heat_source = self.heat_source_rate

        f_tmp[0] = - 1/tau * (temp - temp_in) + 1./rho/cp/vol * heat_source

        return f_tmp