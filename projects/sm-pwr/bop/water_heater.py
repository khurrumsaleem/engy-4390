#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file is part of the Cortix toolkit environment.
# https://cortix.org
"""Cortix module"""

import logging

import unit

import iapws.iapws97 as steam_table

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

        self.port_names_expected = ['inflow', 'outflow']

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

      
        self.heating_inflow_temp = 190 + 273.15
        self.heating_inflow_pressure = 0.1 #MPa
        
        self.heating_outflow_temp = 180 + 273.15
        
        self.inflow_pressure = 3.4
        self.inflow_temp = 20+273
        self.inflow_mass_flowrate = 67
        
        self.tau_heating = 1.5
        self.tau_secondary = 2
        
        

        self.outflow_temp = 20 + 273.15
        self.outflow_temp_ss = 422 #k
        self.outflow_mass_flowrate = 67
        self.outflow_pressure = 3.4
        
        # Heating outflow phase history
        quantities = list()

        temp = Quantity(name='temp',
                        formal_name='T_1', unit='K',
                        value=self.heating_outflow_temp,
                        latex_name=r'$T_1$',
                        info='Steamer Primary Outflow Temperature')

        quantities.append(temp)

        self.heating_outflow_phase = Phase(time_stamp=self.initial_time,
                                           time_unit='s', quantities=quantities)
        
        # Outflow phase history
        quantities = list()

        flowrate = Quantity(name='flowrate',
                            formal_name='q_2', unit='kg/s',
                            value=self.outflow_mass_flowrate,
                            latex_name=r'$q_2$',
                            info='Outflow Mass Flowrate')

        quantities.append(flowrate)

        temp = Quantity(name='temp',
                        formal_name='T_2', unit='K',
                        value=self.outflow_temp,
                        latex_name=r'$T_2$',
                        info='Outflow Temperature')

        quantities.append(temp)

        press = Quantity(name='pressure',
                         formal_name='P_2', unit='MPa',
                         value=self.outflow_pressure,
                         latex_name=r'$P_2$',
                         info='Outflow Pressure')

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
            outflow['flowrate'] = self.outflow_mass_flowrate
            self.send((msg_time, outflow), 'outflow')

        # Interactions in the inflow port
        #----------------------------------------
        # one way "from" condenser

        # receive from
        if self.get_port('inflow').connected_port:

            self.send(time, 'inflow')

            (check_time, inflow) = self.recv('inflow')
            assert abs(check_time-time) <= 1e-6

            self.inflow_temp = inflow['temperature']
            self.inflow_pressure = inflow['pressure']
            self.inflow_mass_flowrate = inflow['mass_flowrate']
            
  

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

        temp_p = u_vec[0] # primary outflow temp
        temp_s = u_vec[1] # secondary outflow temp

        # Get state values
        
        flow_rate = self.inflow_mass_flowrate
        t_exit = self.outflow_temp_ss
        p_out = self.outflow_pressure_ss
        #flow_out = flowturb+flowcond
        flow_out = flowcond
        h_exit = steam_table._Region1(t_out,p_out)['h']
        q_removed = flow_out*h_exit-flowcond*hcondin
        #update state variables
        condenser_outflow = self.outflow_phase.get_row(time)

        time += self.time_step

        self.outflow_phase.add_row(time, flow_out)
        self.outflow_phase.set_value('temp', t_exit, time)
        self.outflow_phase.set_value('flowrate', flow_out, time)
        self.outflow_phase.set_value('pressure', p_out, time)

        return time
    
    def __get_state_vector(self, time):
        """Return a numpy array of all unknowns ordered as shown.
           
        """

        u_vec = np.empty(0, dtype=np.float64)

        temp_p = self.heating_outflow_phase.get_value('temp', time)
        u_vec = np.append(u_vec, temp_p)

        temp_s = self.outflow_phase.get_value('temp', time)
        u_vec = np.append(u_vec, temp_s)
        
    
    def __f_vec(self, u_vec, time):

        temp_p = u_vec[0] # get temperature of primary outflow

        temp_s = u_vec[1] # get temperature of secondary outflow

        # initialize f_vec to zero
        f_tmp = np.zeros(2, dtype=np.float64) # vector for f_vec return

        #-----------------------
        # primary energy balance
        #-----------------------
        rho_p = 1/(steam_table._Region1(self.heating_inflow_temp,self.heating_inflow_pressure/unit.mega)["v"])
        cp_p = steam_table._Region2(self.heating_inflow_temp,self.heating_inflow_pressure/unit.mega)["cp"]
        vol_p = self.primary_volume
        
        
        temp_p_in = self.heating_inflow_temp

        tau_p = self.tau_heating

        #-----------------------
        # secondary energy balance
        #-----------------------
        rho_s = 1/(steam_table._Region1((self.secondary_inflow_temp),self.secondary_inflow_pressure/unit.mega)["v"])
        cp_s = steam_table._Region1((self.secondary_inflow_temp),self.secondary_inflow_pressure/unit.mega)["cp"]
        vol_s = self.secondary_volume

        temp_s_in = self.secondary_inflow_temp

        tau_s = self.tau_secondary
        
        #-----------------------
        # calculations
        #-----------------------
        heat_sink = self.__heat_sink_rate(temp_p_in, temp_s_in, cp_p, cp_s)

        f_tmp[0] = - 1/tau_p * (temp_p - temp_p_in) + 1./rho_p/cp_p/vol_p * heat_sink
        
        
        


        heat_source = - heat_sink
        
        temp_s_out = - 1/tau_s * (temp_s - temp_s_in) + 1./rho_s/cp_s/vol_s * heat_source

        q_total = (temp_s_out-temp_s_in)*cp_s
        q_heat = (self.t_sat - temp_s_in)*cp_s
       
        if q_total < q_heat:
            f_tmp[1] = temp_s_out
            
        else:
            h_v = steam_table._Region4(self.secondary_inflow_pressure,1)['h']
            h_l = steam_table._Region4(self.secondary_inflow_pressure,0)['h']
            h_vap = h_v-h_l
            
            quality = (h_vap - (q_total - q_heat))/h_vap
            f_tmp[1] = steam_table._Region4(self.secondary_inflow_pressure,quality)['T']
        
        
        #print(temp_p_in,f_tmp[0], f_tmp[1])
        #f_tmp[1] = - 1/tau_s * (temp_s - temp_s_in) + 1./rho_s/cp_s/vol_s * heat_source

        return f_tmp

    def __heat_sink_rate(self, temp_p, temp_s, cp_p, cp_s):
        """Cooling rate of primary."""

        c_min = cp_s*self.secondary_inflow_mass_flowrate
        
        #ntu = self.ht_coeff/(c_min)
        ntu = 4.5
        
        c_r = (c_min)/(cp_p*self.heating_inflow_mass_flowrate)
        
        
        
        eta = (1-np.exp(-ntu*(1-c_r)))/(1-c_r*np.exp(-ntu*(1-c_r)))
        
       
        q_p = - eta*c_min*(temp_p - temp_s)
        
        


        return q_p

