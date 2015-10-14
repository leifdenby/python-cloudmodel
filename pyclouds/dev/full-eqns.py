
# coding: utf-8

# In[1]:

import numpy as np
import matplotlib
matplotlib.use('Agg')


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

environment = ambient_stratification

# In[5]:

z_max = 10e3
z_points = numpy.linspace(100., z_max, 100)

w0 = 0.1
T0 = T_e(0.0) + 0.2
# p, r, w, T, q_v, q_r, q_l, q_i
initial_condition = [500.0, w0, T0, 0.01, 0., 0., 0.]


# In[6]:

stopping_criterion = None

beta = 0.0
d_beta = 0.02
D = 0.0
dD = 0.0
gamma = 1.0

microphysics = cloud_microphysics.MoistAdjustmentMicrophysics()

default_kwargs = {
    'environment': environment,
    'microphysics': microphysics, 
}

integration_kwargs = {
    'initial_condition': initial_condition,
    'z': z_points,
}

profiles = []

cloud_model = cloud_equations.DryAirOnly(gamma=gamma, D=D, beta=beta, **default_kwargs)
profiles.append(cloud_model.integrate(**integration_kwargs))

cloud_model = cloud_equations.DryAirOnly(gamma=gamma, D=D+dD, beta=beta+d_beta, **default_kwargs)
profiles.append(cloud_model.integrate(**integration_kwargs))

# cloud_model = cloud_equations.DryAirOnly(gamma=gamma, D=D+dD, beta=beta+2*d_beta, **default_kwargs)
# profiles.append(cloud_model.integrate(**integration_kwargs))

cloud_model = cloud_equations.FullThermodynamicsCloudEquations(gamma=gamma, D=D, beta=beta, **default_kwargs)
profiles.append(cloud_model.integrate(**integration_kwargs))

cloud_model = cloud_equations.FullThermodynamicsCloudEquations(gamma=gamma, D=D+dD, beta=beta+d_beta, **default_kwargs)
profiles.append(cloud_model.integrate(initial_condition, z_points))

fig = plotting.plot_profiles(profiles, variables=['r', 'w', 'T', 'q_v', 'mse', 'T__tephigram',])
fig.savefig('/tmp/cloud-profile-comparison.png')
