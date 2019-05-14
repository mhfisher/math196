from bokeh.io import export_png
from bokeh.plotting import figure, output_file, show

from ode import *

parameters = {
    'a': 1.0,
    'h': 1.0,
    'rho': 10.0,
    'R': 50.0,
    'g': 9.8,
    'l': 3.0
}

x_range = np.arange(0.0, parameters['a'], 0.01)
catenary_parameters = {key: value for key,value in parameters.iteritems()}
catenary_parameters['T0'] = T0_catenary(**parameters)
roadway_parameters = {key: value for key,value in parameters.iteritems()}
roadway_parameters['T0'] = T0_roadway(**parameters)
bridge_parameters = {key: value for key,value in parameters.iteritems()}
bridge_parameters['T0'] = T0_bridge(**parameters)

curves = [
    ('Catenary closed-form',) + catenary(**catenary_parameters),
    ('Negligible chain mass closed-form',) + roadway(**roadway_parameters),
    # ('Catenary ODE',) + solver(x_range, case='catenary', **catenary_parameters),
    # ('Negligible chain mass ODE',) + solver(x_range, case='roadway', **roadway_parameters),
    ('Suspension Bridge ODE',) + solver(x_range, case='bridge', **bridge_parameters),
]


plot(curves, match_first_h=True)
