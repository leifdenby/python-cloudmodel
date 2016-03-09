all: spec_conc_eqns.f90
	f2py -c spec_conc_eqns.f90 -m f_cloudequations

tests: all
	py.test
