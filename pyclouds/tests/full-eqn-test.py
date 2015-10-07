
# coding: utf-8

# In[1]:

import numpy as np
import matplotlib.pyplot as plot


# # Ambient state

# In[2]:

from pycfd.reference.atmospheric_flow import stratification_profiles, saturation_calculation

ambient_stratification = stratification_profiles.Soong1973Dry()
T_e = lambda z: ambient_stratification.temp(z)
p_e = ambient_stratification.p


# # Model integration

# In[3]:

import odespy, numpy
from pyclouds import cloud_equations, plotting, cloud_microphysics


# In[4]:

g = ambient_stratification.g
R_d = ambient_stratification.gas_properties.R()
cp_d = ambient_stratification.gas_properties.cp()

environment = {'T_e': T_e, 'p_e': p_e, 'rho': ambient_stratification.rho, }
constants = {'R_d': R_d, 'cp_d': cp_d, 'g': g }


# In[5]:

z_points = numpy.linspace(100., 4e3, 500)

w0 = 0.1
T0 = T_e(0.0) + 0.2
# p, r, w, T, q_v, q_r, q_l, q_i
initial_condition = [101325.0, 500.0, w0, T0, 0.00, 0., 0., 0.]


# In[6]:

stopping_criterion = lambda F, z, k: k > 4

microphysics = cloud_microphysics.DummyMicrophysics()

cloud_model = cloud_equations.DryAirOnly(environment=environment, constants=constants, gamma=0.0, D=0.0, beta=0.0, microphysics=microphysics)
F, z = cloud_model.integrate(initial_condition, z_points, SolverClass=odespy.Euler, stopping_criterion=stopping_criterion)

print 

cloud_model = cloud_equations.FullThermodynamicsCloudEquations(environment=environment, constants=constants, gamma=0.0, D=0.0, beta=0.0, microphysics=microphysics)
F, z = cloud_model.integrate(initial_condition, z_points, SolverClass=odespy.Euler, stopping_criterion=stopping_criterion)
plotting.profile_plot(F, z, Te=T_e)
