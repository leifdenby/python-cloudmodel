# Python 1D cloud-model

The `pyclouds` package is a tool to implement and investigate new 1D
cloud-models, with particular emphasis on improving the cloud model in
CCFM ([Convective Cloud Field Model - Nober & Graf
2005](http://www.atmos-chem-phys.net/5/2749/2005/acp-5-2749-2005.html)).
Derivation of the cloud equations can be found in [L Denby
2017](https://www.repository.cam.ac.uk/handle/1810/269850).

Each model is based around a set of conservation equations describing the
vertical change in in-cloud variables in a convective cloud in the mature state.
Once specified each cloud-model can then be integrated, which is done
using one of [scipy's ODE integrators](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp)
for initial value problems. By default Runge-Kutta-Fehlberg is used due to
its adaptive error-correcting timestepping.

An example of an integrated cloud of a given radius follows below:

![Cloud vertical radius, temperature and vertical velocity](doc/cloud_profile.png)

![Cloud vertical profile of specific concentration of water vapor, liquid water
and rain water](doc/cloud_hydrometeors.png)

## Getting started

### Installing `pyclouds`

Depending on whether you wish to modify the cloud-models included with
`pyclouds` or whether you simply wish to use the cloud-equations as they
are you will either want to install `pyclouds` from a local copy or
directly from pipy (or github).

In either case it is advised that you create an isolated python
environment into which you install `pyclouds`. This can be achieved with
a tool like [conda](https://docs.conda.io/en/latest/miniconda.html#installing).

#### Installing with aim to modify `pyclouds`

To be able to modify the `pyclouds` source-code and interact with it you will
need to clone the `pyclouds` repository on github (or your fork on github)
to your local computer and install from that. That will ensure that any
changes you make locally will be available when you import the `pyclouds`
package:

```
git clone https://github.com/leifdenby/python-cloudmodel
cd python-cloudmodel
python -m pip install -e .
```

The `-e` flag instructs `pip` to install `pyclouds` so that any changes
you make to the source directory are available when you import `pyclouds`.
You'll also see the version of `pyclouds` (i.e. `pyclouds.__version__`)
includes the word "dirty" once you've made changes.

#### Installing with aim to use `pyclouds` as-is

You can simply install the most recent tagged version of `pyclouds`
directly from [pipy](https://pypi.org/) with `pip`:

```bash
python -m pip install pyclouds
```

Or install the `master` branch (note this may be ahead of the most
recently tagged version) on github

```bash
python -m pip install git+https://github.com/leifdenby/python-cloudmodel#egg=pyclouds
```

### Using `pyclouds`

The best way to get an idea of how to use `pyclouds` is to look at the
included Jupyter notebooks in [notebooks/](notebooks/) or the tests
(filenames starting with `test_`) [tests/](tests/).

## Contributing to `pyclouds`

If you spot any mistakes or have made any additions please make a
pull-request. Thank you!
