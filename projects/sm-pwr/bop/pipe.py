#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""Cortix module"""

import logging
import math
import unit
import numpy as np

from iapws import IAPWS97 as WaterProps
from cortix import Module
from cortix.support.phase_new import PhaseNew as Phase
from cortix import Quantity
from scipy.integrate import odeint

class Pipe(Module):
    """Pipe Module.

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
        self.end_time = 6000*unit.second
        self.time_step = 1.0*unit.second

        self.show_time = (False, 10.0*unit.second)
        self.save = True

        self.log = logging.getLogger('cortix')
        self.__logit = True # flag indicating when to log

        # Domain attributes
        self.f_mood = 0
        self.pressure_mood = 34*unit.bar
        
        self.f_sj = 0
        self.pressure_sj = 34*unit.bar
        
        self.f_chen = 0
        self.pressure_chen = 34*unit.bar
        
        self.f_man = 0
        self.pressure_man = 34*unit.bar
        
        self.f_all = 0
        self.pressure_all = 34*unit.bar
        
        # Configuration parameters
        
        self.diameter = 0.2*unit.meter
        self.thickness = 10*unit.milli*unit.meter
        self.length = 10*unit.meter
        self.pipe_shell_thermal_conductivity = 45*unit.watt/unit.k/unit.meter
        self.pipe_roughness = .061*unit.milli*unit.meter
        
        self.outside_air_temperature = (25+273)*unit.k
        self.outside_forced_airflow = False
        self.outside_airflow_velocity = 0*unit.meter/unit.second
        self.outside_air_pressure = 1*unit.atm
        

        # Initialization

        self.inflow_temp = (100+273)*unit.K
        self.inflow_pressure = 34*unit.bar
        self.inflow_mass_flowrate = 0*unit.kg/unit.second

        self.outflow_temp = 50 + 273.15
        self.outflow_pressure = 34.0*unit.bar
        self.outflow_mass_flowrate = self.inflow_mass_flowrate

        # Inflow phase history
        quantities = list()

        flowrate = Quantity(name='flowrate',
                            formal_name='q_1', unit='kg/s',
                            value=self.inflow_mass_flowrate,
                            latex_name=r'$q_1$',
                            info='Pipe Inflow Mass Flowrate')

        quantities.append(flowrate)

        temp = Quantity(name='temp',
                        formal_name='T_1', unit='K',
                        value=self.inflow_temp,
                        latex_name=r'$T_1$',
                        info='Pipe Inflow Temperature')

        quantities.append(temp)

        press = Quantity(name='pressure',
                         formal_name='P_2', unit='Pa',
                         value=self.inflow_pressure,
                         latex_name=r'$P_2$',
                         info='Pipe inflow Pressure')
        quantities.append(press)
        self.inflow_phase = Phase(time_stamp=self.initial_time,
                                  time_unit='s', quantities=quantities)

        # Outflow phase history
        quantities = list()

        flowrate = Quantity(name='flowrate',
                            formal_name='q_2', unit='kg/s',
                            value=self.outflow_mass_flowrate,
                            latex_name=r'$q_2$',
                            info='Pipe Outflow Mass Flowrate')

        quantities.append(flowrate)

        temp = Quantity(name='temp',
                        formal_name='T_2', unit='K',
                        value=self.outflow_temp,
                        latex_name=r'$T_2$',
                        info='Pipe Outflow Temperature')

        quantities.append(temp)

        press = Quantity(name='pressure',
                         formal_name='P_2', unit='Pa',
                         value=self.outflow_pressure,
                         latex_name=r'$P_2$',
                         info='Pipe Outflow Pressure')

        quantities.append(press)

        self.outflow_phase = Phase(time_stamp=self.initial_time,
                                   time_unit='s', quantities=quantities)
        
        # friction factors and pressure drops
        quantities = list()
        
        friction = Quantity(name='moodf',
                            formal_name='f_mood', unit='',
                            value=self.f_mood,
                            latex_name=r'$f_mood$',
                            info='moody friction factor')
        quantities.append(friction)          
        press = Quantity(name='moodp',
                         formal_name='P_mood', unit='Pa',
                         value=self.pressure_mood,
                         latex_name=r'$P_mood$',
                         info='Moody Pressure')
        quantities.append(press)
        
        friction = Quantity(name='sjf',
                            formal_name='f_sj', unit='',
                            value=self.f_sj,
                            latex_name=r'$f_sj$',
                            info='sj friction factor')
        quantities.append(friction)  
        press = Quantity(name='sjp',
                         formal_name='P_sj', unit='Pa',
                         value=self.pressure_sj,
                         latex_name=r'$P_sj$',
                         info='sj Pressure')
        quantities.append(press) 
        
        friction = Quantity(name='chenf',
                            formal_name='f_chen', unit='',
                            value=self.f_chen,
                            latex_name=r'$f_chen$',
                            info='chen friction factor')
        quantities.append(friction)  
        press = Quantity(name='chenp',
                         formal_name='P_chen', unit='Pa',
                         value=self.pressure_chen,
                         latex_name=r'$P_chen$',
                         info='chen Pressure')
        quantities.append(press)   

        friction = Quantity(name='manf',
                            formal_name='f_man', unit='',
                            value=self.f_man,
                            latex_name=r'$f_man$',
                            info='man friction factor')
        quantities.append(friction)  
        press = Quantity(name='manp',
                         formal_name='P_man', unit='Pa',
                         value=self.pressure_man,
                         latex_name=r'$P_man$',
                         info='man Pressure')
        quantities.append(press)
        
        friction = Quantity(name='allf',
                            formal_name='f_all', unit='',
                            value=self.f_all,
                            latex_name=r'$f_all$',
                            info='all friction factor')
        quantities.append(friction)  
        press = Quantity(name='allp',
                         formal_name='P_all', unit='Pa',
                         value=self.pressure_all,
                         latex_name=r'$P_all$',
                         info='all Pressure')
        quantities.append(press)
        
        self.state_phase = Phase(time_stamp=self.initial_time,
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

        temp = self.inflow_temp
        pressure = self.inflow_pressure
        waterP = WaterProps(T=400,P=12.8)
         
        print(time)
        flowrate = time/self.end_time * 1000*unit.kg/unit.second
        rho = waterP.rho #units kg/m3
        vel = flowrate/(math.pi/4*self.diameter**2 * rho) #m/s       
        frictions = self.__friction(waterP,flowrate)
        
        p_mood = rho/2*vel**2/self.diameter*frictions[0]
        p_sj = rho/2*vel**2/self.diameter*frictions[1]
        p_chen = rho/2*vel**2/self.diameter*frictions[2]
        p_man = rho/2*vel**2/self.diameter*frictions[3]
        p_all = rho/2*vel**2/self.diameter*frictions[4]
        # Update state variables
        state = self.state_phase.get_row(time)
        inflow = self.inflow_phase.get_row(time)
        outflow = self.outflow_phase.get_row(time)

        time += self.time_step

        self.state_phase.add_row(time, state)
        self.state_phase.set_value('moodf', frictions[0], time)
        self.state_phase.set_value('moodp', p_mood, time)
        self.state_phase.set_value('sjf', frictions[1], time)
        self.state_phase.set_value('sjp', p_sj, time)
        self.state_phase.set_value('chenf', frictions[2], time)
        self.state_phase.set_value('chenp', p_chen, time)
        self.state_phase.set_value('manf', frictions[3], time)
        self.state_phase.set_value('manp', p_man, time)
        self.state_phase.set_value('allf', frictions[4], time)
        self.state_phase.set_value('allp', p_all, time)

        self.inflow_phase.add_row(time, inflow)
        self.inflow_phase.set_value('temp', self.inflow_temp, time)
        self.inflow_phase.set_value('flowrate', flowrate , time)
        self.inflow_phase.set_value('pressure', self.inflow_pressure, time)
        
        self.outflow_phase.add_row(time, outflow)
        self.outflow_phase.set_value('temp', self.outflow_temp, time)
        self.outflow_phase.set_value('flowrate', flowrate , time)
        self.outflow_phase.set_value('pressure', self.outflow_pressure, time)
        
        return time
        
    def __friction(self,waterP,flowrate):
        #considering only one phase flow
        d = self.diameter 
        mu = waterP.mu  #dyamic viscosity (pa*s)
        Re = 4*flowrate/(1*math.pi*d*mu) #mew is the viscosity
        e = self.pipe_roughness
        
        #for Re = 4e3-1e8 and e/D = 0-1e-2
        if Re > 4000 and e/d < .01:
            #Moody 1947
            fmood = .0055*(1+(2*10**4*e/d+10**6/Re)**(1/3))
        else:
            fmood = 0

        #for Re = 5e3-1e8 and e/D = 1e-6-5e-2
        if Re > 5000 and e/d < .05:
            #Swamee and Jain 1976
            fsj = .25/(math.log10(e/d/3.7 + 5.74/(Re**0.9)))**2
        else:
            fsj = 0
            
        #for Re = 4e3-4e8
        if Re > 4000:
            #Chen 1979
            one_over_sqrt_f = -2*math.log10(e/d/3.7065 - 5.0452/Re*math.log10(1/2.8257*(e/d)**1.1098 + 5.8506/Re**0.8981))
            fchen = (1/one_over_sqrt_f)**2
        else:
            fchen = 0
        
        #for Re = 4e3-1e8 and e/D = 0-5e-2
        if Re > 4000 and e/d < .05:
            #Manadilli 1997
            one_over_sqrt_f = -2*math.log10(e/d/3.7 + 95/Re**0.983 - 96.82/Re)
            fman = (1/one_over_sqrt_f)**2
        else:
            fman = 0
            
        
        #for all flow regimes
        if Re > 4000 and e/d < .05:
            #Cheng 2008
            a = 1/(1+(Re/2720)**9)
            b = 1/(1+(Re/(160*d/e))**2)
            one_over_f = (Re/64)**a * (1.8*math.log10(Re/6.8))**(2*(1-a)*(b)) * (2.0*math.log10(3.7*d/e))**(2*(1-a)*(1-b))
            fall = 1/one_over_f
        else:
            fall = 0
            
        return fmood, fsj, fchen, fman, fall
            
            
            
            
            
        