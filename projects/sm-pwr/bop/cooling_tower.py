#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""Cortix module"""

import logging

import unit

from iapws import IAPWS97 as steam_table

from cortix import Module
from cortix.support.phase_new import PhaseNew as Phase
from cortix import Quantity

class CoolingTower(Module):
    """Cooling Tower Module.

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

        self.port_names_expected = ['inflow', 'outflow']

        # Units

        # General attributes
        self.initial_time = 0.0*unit.second
        self.end_time = 1.0*unit.hour
        self.time_step = 10.0*unit.second

        self.show_time = (False, 10.0*unit.second)
        self.save = True

        self.log = logging.getLogger('cortix')
        self.__logit = True # flag indicating when to log

        # Domain attributes

        # Configuration parameters

        # Initialization

        self.inflow_temp = (20+273)*unit.K
        #self.inflow_pressure = 34*unit.bar
        self.inflow_pressure = 0.008066866*unit.mega*unit.pascal
        self.inflow_mass_flowrate = 67*unit.kg/unit.second

        self.outflow_temp = 50 + 273.15
        self.outflow_pressure = 34.0*unit.bar
        self.outflow_mass_flowrate = self.inflow_mass_flowrate

        # Inflow phase history
        quantities = list()

        flowrate = Quantity(name='flowrate',
                            formal_name='q_1', unit='kg/s',
                            value=self.outflow_mass_flowrate,
                            latex_name=r'$q_1$',
                            info='Condenser Inflow Mass Flowrate')

        quantities.append(flowrate)

        temp = Quantity(name='temp',
                        formal_name='T_1', unit='K',
                        value=self.outflow_temp,
                        latex_name=r'$T_1$',
                        info='Condenser Inflow Temperature')

        quantities.append(temp)

        self.inflow_phase = Phase(time_stamp=self.initial_time,
                                  time_unit='s', quantities=quantities)

        # Outflow phase history
        quantities = list()

        flowrate = Quantity(name='flowrate',
                            formal_name='q_2', unit='kg/s',
                            value=self.outflow_mass_flowrate,
                            latex_name=r'$q_2$',
                            info='Condenser Outflow Mass Flowrate')

        quantities.append(flowrate)

        temp = Quantity(name='temp',
                        formal_name='T_2', unit='K',
                        value=self.outflow_temp,
                        latex_name=r'$T_2$',
                        info='Condenser Outflow Temperature')

        quantities.append(temp)

        press = Quantity(name='pressure',
                         formal_name='P_2', unit='Pa',
                         value=self.outflow_pressure,
                         latex_name=r'$P_2$',
                         info='Condenser Outflow Pressure')

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

            # Evolve one time step
            #---------------------

            time = self.__step(time)

            # Communicate information
            #------------------------
            self.__call_ports(time)

        self.end_time = time # correct the final time if needed

    def __call_ports(self, time):


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

        # Interactions in the outflow port
        #-----------------------------------------
        # one way "to" outflow

        # send to
        if self.get_port('outflow').connected_port:

            msg_time = self.recv('outflow')
            assert msg_time <= time

            temp = self.outflow_phase.get_value('temp', msg_time)
            pressure = self.outflow_phase.get_value('pressure', msg_time)

            outflow = dict()
            outflow['temperature'] = temp
            outflow['pressure'] = pressure
            self.outflow_mass_flowrate = self.inflow_mass_flowrate
            outflow['mass_flowrate'] = self.outflow_mass_flowrate
            self.send((msg_time, outflow), 'outflow')

    def __step(self, time=0.0):


        # Update state variables
        tower_outflow = self.outflow_phase.get_row(time)
        tower_inflow = self.inflow_phase.get_row(time)

        time += self.time_step

        self.inflow_phase.add_row(time, tower_inflow)
        self.inflow_phase.set_value('temp', self.inflow_temp, time)
        self.inflow_phase.set_value('flowrate', self.inflow_mass_flowrate , time)

        self.outflow_phase.add_row(time, tower_outflow)
        self.outflow_phase.set_value('temp', self.outflow_temp, time)
        self.outflow_phase.set_value('flowrate', self.outflow_mass_flowrate , time)
        self.outflow_phase.set_value('pressure', self.outflow_pressure, time)

        return time

    def __get_state_vector(self, time):
        """Return a numpy array of all unknowns ordered as shown.
        """

        u_vec = np.empty(0, dtype=np.float64)

        temp_out = self.outflow_phase.get_value('temp', time)
        flowrate_out = self.outflow_phase.get_value('flowrate', time)
        u_vec = np.append(u_vec, temp)
        u_vec = np.append(u_vec, flowrate)
        
        return  u_vec

    def __f_vec(self, u_vec, time):

        temp = u_vec[0]
        flowrate = u_vec[0]
        # initialize f_vec to zero
        f_tmp = np.zeros(2, dtype=np.float64) # vector for f_vec return

        #-----------------------
        # primary energy balance
        #-----------------------
        water_out = steam_table(T=temp,
                            P=self.inflow_pressure/unit.mega/unit.pascal)
        water_in = steam_table(T=self.inflow_temp,
                            P=self.inflow_pressure/unit.mega/unit.pascal)
        temp_avg = (temp+self.inflow_temp)/2
        water_avg = steam_table(T=temp_avg,
                            P=self.inflow_pressure/unit.mega/unit.pascal)
        hvap = 2256.4+ (2500.9-2256.4)*(temp_avg/100)
        
        
        assert water.phase != 'Vapour'

        cp_in = water_in.Liquid.cp*unit.kj/unit.kg/unit.K
        cp_out = water_out.Liquid.cp*unit.kj/unit.kg/unit.K
        cp_avg = water_avg.Liquid.cp*unit.kj/unit.kg/unit.K
        vol = self.volume
        
        temp_in = self.inflow_temp

        tau = vol/(self.inflow_mass_flowrate/rho)

        #-----------------------
        # calculations
        #-----------------------
        heat_source_pwr = self.external_heat_source_rate + \
                          self.electric_heat_source_rate

        heat_source_pwr_dens = heat_source_pwr/vol

        f_tmp[0] = - 1/tau * (temp - temp_in) + 1./rho/cp * heat_source_pwr_dens

        return f_tmp
