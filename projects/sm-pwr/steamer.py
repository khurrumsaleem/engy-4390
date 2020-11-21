#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""Cortix module"""

import logging

import scipy.constants as unit

import iapws.iapws97 as steam_table

from cortix import Module
from cortix.support.phase_new import PhaseNew as Phase
from cortix import Quantity

class Steamer(Module):
    """Steam generator.

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

        #self.port_names_expected = ['primary-inflow', 'primary-outflow',
        #                            'secondary-inflow', 'secondary-outflow']
        self.port_names_expected = ['primary-inflow', 'primary-outflow',
                                    'secondary-outflow']
        # Units
        unit.kg = unit.kilo*unit.gram
        unit.meter = 1.0
        unit.cm = unit.centi*unit.meter
        unit.second = 1.0
        unit.pascal = 1.0
        unit.joule = 1.0
        unit.kj = unit.kilo*unit.joule
        unit.kelvin = 1.0
        unit.watt = 1.0

        # General attributes
        self.initial_time = 0.0*unit.second
        self.end_time = 1.0*unit.hour
        self.time_step = 10.0*unit.second

        self.show_time = (False, 10.0*unit.second)

        self.log = logging.getLogger('cortix')
        self.__logit = True # flag indicating when to log

        # Domain attributes

        # Configuration parameters

        # Initialization

        self.primary_inflow_pressure = 12.8
        self.primary_inflow_temp = 20 + 273.15
        self.primary_inflow_mass_flowrate = 666

        self.secondary_inflow_pressure = 3.5
        self.secondary_inflow_temp = 20 + 273.15
        self.secondary_inflow_mass_flowrate = 0.0

        self.primary_outflow_temp = 20 + 273.15
        self.primary_outflow_mass_flowrate = 666
        self.primary_outflow_pressure = 12.8

        self.secondary_outflow_mass_flowrate = 67
        self.secondary_outflow_temp = 20 + 273.15
        self.secondary_outflow_pressure = 3.4

        # Primary outflow phase history
        quantities = list()

        temp = Quantity(name='temp',
                        formal_name='T_1', unit='K',
                        value=self.primary_outflow_temp,
                        latex_name=r'$T_1$',
                        info='Primary Outflow Temperature')

        quantities.append(temp)

        self.primary_outflow_phase = Phase(time_stamp=self.initial_time,
                                           time_unit='s', quantities=quantities)

        # Secondary outflow phase history
        quantities = list()

        flowrate = Quantity(name='flowrate',
                            formal_name='q_2', unit='kg/s',
                            value=self.secondary_outflow_mass_flowrate,
                            latex_name=r'$q_2$',
                            info='Secondary Outflow Mass Flowrate')

        quantities.append(flowrate)

        temp = Quantity(name='temp',
                        formal_name='T_2', unit='K',
                        value=self.secondary_outflow_temp,
                        latex_name=r'$T_2$',
                        info='Secondary Outflow Temperature')

        quantities.append(temp)

        press = Quantity(name='pressure',
                         formal_name='P_2', unit='Pa',
                         value=self.secondary_outflow_pressure,
                         latex_name=r'$P_2$',
                         info='Secondary Outflow Pressure')

        quantities.append(press)

        self.secondary_outflow_phase = Phase(time_stamp=self.initial_time,
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

        # Interactions in the primary-outflow port
        #-----------------------------------------
        # one way "to" primary-outflow

        # send to
        if self.get_port('primary-outflow').connected_port:

            msg_time = self.recv('primary-outflow')

            temp = self.primary_outflow_phase.get_value('temp', msg_time)
            primary_outflow = dict()
            primary_outflow['temperature'] = temp
            primary_outflow['pressure'] = self.primary_outflow_pressure
            primary_outflow['flowrate'] = self.primary_outflow_mass_flowrate
            self.send((msg_time, primary_outflow), 'primary-outflow')

        # Interactions in the secondary-outflow port
        #-----------------------------------------
        # one way "to" secondary-outflow

        # send to
        if self.get_port('secondary-outflow').connected_port:

            msg_time = self.recv('secondary-outflow')

            temp = self.secondary_outflow_phase.get_value('temp', msg_time)
            press = self.secondary_outflow_phase.get_value('pressure', msg_time)
            flowrate = self.secondary_outflow_phase.get_value('flowrate', msg_time)
            secondary_outflow = dict()
            secondary_outflow['temperature'] = temp
            secondary_outflow['pressure'] = press
            secondary_outflow['mass_flowrate'] = flowrate

            self.send((msg_time, secondary_outflow), 'secondary-outflow')

        # Interactions in the primary-inflow port
        #----------------------------------------
        # one way "from" primary-inflow

        # receive from
#        if self.get_port('primary-inflow').connected_port:
#
#            self.send(time, 'primary-inflow')
#
#            (check_time, primary_inflow) = self.recv('primary-inflow')
#            assert abs(check_time-time) <= 1e-6
#
#            self.primary_inflow_temp = primary_inflow['temperature']
#            self.primary_inflow_pressure = primary_inflow['pressure']
#            self.primary_inflow_mass_flowrate = primary_inflow['mass_flowrate']

        # Interactions in the secondary-inflow port
        #----------------------------------------
        # one way "from" secondary-inflow

        # receive from
        if self.get_port('secondary-inflow').connected_port:

            self.send(time, 'secondary-inflow')

            (check_time, secondary_inflow) = self.recv('secondary-inflow')
            assert abs(check_time-time) <= 1e-6

            self.secondary_inflow_temp = secondary_inflow['temperature']
            self.secondary_inflow_pressure = secondary_inflow['pressure']
            self.secondary_inflow_mass_flowrate = secondary_inflow['temperature']

    def __step(self, time=0.0):

        '''
        p_out_sg = params['pressure_sg_exit']
        t_in_sg = params['t_sg_entrance']
        m_dot_1 = params['coolant_mass_flowrate']
        m_dot_2 = params['mdot_secondary']
        eta_sg = params['sg_efficiency']
        p_primary = params['primary_pressure']
        u_overall = params['heat_transfer_coefficient_garbage']

        #If steady state
        if primary_inflow_temp > 321 + 273.15:
            primary_outflow_temp = 265 + 273.15
            t_exit_sg = 307 + 273.15
            return primary_outflow_temp, t_exit_sg, p_out_sg

        q_pls = u_overall * (primary_inflow_temp - t_in_sg
                             )  #Garbage overall U approximation

        cp_primary = steam_table._Region1((primary_inflow_temp),
                                          p_primary)["cp"]  #Approximated Cp

        primary_outflow_temp = primary_inflow_temp - (
            q_pls /
            (m_dot_1 * cp_primary))  #Temp from Garbage equ. and approxiamted Cp

        cp_primary = steam_table._Region1(
            (primary_outflow_temp + primary_inflow_temp) / 2,
            p_primary)["cp"]  #better Cp

        q_primary = m_dot_1 * (primary_inflow_temp -
                               primary_outflow_temp) * cp_primary  #Better Q

        h_1_sg = steam_table._Region1(148.89 + 273.15, 3.5190841)['h']
        h_2_sg = (eta_sg * q_primary / m_dot_2) + h_1_sg  #Enthalpy balance

        ### Calculations pertaining to SG exiting conditions ###

        bubl = steam_table._Region4(p_out_sg, 0)
        dew = steam_table._Region4(p_out_sg, 1)
        bubl_entropy = bubl['s']
        dew_entropy = dew['s']
        bubl_enthalpy = bubl['h']
        dew_enthalpy = dew['h']

        #if fluid exiting sg is subcooled
        if h_2_sg < bubl_enthalpy:

            t_exit_sg = steam_table._Backward1_T_Ph(p_out_sg, h_2_sg)
            return primary_outflow_temp, t_exit_sg, p_out_sg

        # if fluid exiting sg is superheated
        if h_2_sg > dew_enthalpy:

            t_exit_sg = steam_table._Backward2_T_Ph(p_out_sg, h_2_sg)

            return primary_outflow_temp, t_exit_sg, p_out_sg
        # else fluid is two phase
        t_exit_sg = steam_table._TSat_P(p_out_sg)

        #return primary_outflow_temp, t_exit_sg, p_out_sg
        '''

        # temporary to get ports tested

        primary_outflow = self.primary_outflow_phase.get_row(time)

        time += self.time_step

        self.primary_outflow_phase.add_row(time, primary_outflow)


        return time
