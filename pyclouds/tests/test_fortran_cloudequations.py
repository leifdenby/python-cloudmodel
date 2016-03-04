import numpy as np
import inspect
import pytest

import pyclouds.cloud_equations
import pyclouds.f_cloudequations

from pycfd.reference.atmospheric_flow import stratification_profiles

@pytest.fixture
def f_ver():
    return pyclouds.f_cloudequations.spec_conc_eqns


@pytest.fixture
def py_ver(f_ver):
    gamma, beta, D = f_ver.m__gamma, f_ver.m__beta, f_ver.m__d

    constant_vars = ['cp_d', 'cp_v', 'L_v', 'L_s', 'rho_l', 'rho_i', 'g', 'R_d', 'R_v']

    constants = dict([(v, getattr(f_ver, v.lower())) for v in constant_vars ])

    ModelClass_py = pyclouds.cloud_equations.FullThermodynamicsCloudEquations

    return ModelClass_py(gamma=gamma, D=D, beta=beta, environment=None, constants=constants)


def test_derivatives(py_ver, f_ver):

    func_names = ['dw_dz', 'dT_dz', 'dr_dz']

    for func_name in func_names:
        py_func = getattr(py_ver, func_name)
        f_func = getattr(f_ver, func_name.lower())
        
        py_args = [a for a in inspect.getargspec(py_func).args if a != 'self']
        # Fortran in case insensitive
        f_args = [s.lower() for s in py_args]

        n_vars = len(py_args)

        # make a random state
        F = np.random.random((n_vars,))

        res_py = py_func(**dict(zip(py_args, F)))
        res_f = f_func(**dict(zip(f_args, F)))

        assert res_py == res_f
