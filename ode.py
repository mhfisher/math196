""" Differential equations for various cases """
from bokeh.plotting import figure, output_file, show
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

# USE THE SOLVED FORM OF T0 (and maybe weighted average for bridge?)

""" Differential Equations """
def bridge_derivative(Y, x, rho=1, R=1, g=9.8, T0=1):
    """
    Returns Y' = [y2, y2'] where Y = [y1, y2], y1 = u(x), y2 = u'(x)
    Both the mass of the chain and roadway are considered.
    """
    y1, y2 = Y
    y2_prime = (rho*np.sqrt(1 + y2**2) + R) * (g / T0)
    return [y2, y2_prime]


def catenary_derivative(Y, x, rho, R, g, T0):
    """
    Returns Y' = [y2, y2'] where Y = [y1, y2], y1 = u(x), y2 = u'(x)
    Only rho, the mass of the chain, is considered.
    """
    y1, y2 = Y
    y2_prime = ((g * rho) / T0) * np.sqrt(1 + y2**2)
    return [y2, y2_prime]


def roadway_derivative(Y, x, rho, R, g, T0):
    """
    Returns Y' = [y2, y2'] where Y = [y1, y2], y1 = u(x), y2 = u'(x)
    Only R, the mass of the roadway, is considered.
    """
    y1, y2 = Y
    y2_prime =  R * (g / T0)
    return [y2, y2_prime]


def solver(a=1.0, h=1.0, rho=1.0, R=1.0, g=9.8, T0=9.8, case=None):
    """
    Graph the right side of the curve described by the model.

    case must be one of {'catenary', 'roadway', 'bridge'}
    """
    derivatives = {
        'catenary': catenary_derivative,
        'roadway': roadway_derivative,
        'bridge': bridge_derivative
    }
    assert case in derivatives.keys(), "Please provide a valid case"

    x_range = np.arange(0, a, 0.01)
    U_x = lambda Y, x: derivatives[case](Y, x, rho, R, g, T0)
    solution = integrate.odeint(U_x, [0, 0], x_range)
    return (x_range, solution[:, 0])

""" Closed-form equations """
def catenary(a=1.0, h=1.0, rho=1.0, R=1.0, g=9.8, T0=9.8):
    """
    Return the closed form of the catenary (hyperbolic cosine)
    """
    C = (g * rho) / T0
    x_range = np.arange(0, a, 0.01)

    return (x_range, [(1/C) * (np.cosh(C*x) - np.cosh(C*a)) + h for x in x_range])


def roadway(a=1.0, h=1.0, rho=1.0, R=1.0, g=9.8, T0=9.8):
    """
    Return the closed form solution of the negligible chain mass (x^2)
    """
    x_range = np.arange(0, a, 0.01)
    return (x_range, [((g*R)/(2*T0))*(x**2 - a**2) + h for x in x_range])


""" Plotting utilities"""
def plot(curves, match_first_h=False):
    """
    Use Bokeh to plot the provided curves. If match_first_h is set, the value
    of the curves at a will be rectified to be the first curve's valeu at a.
    curves should be a list where each element is a tuple (name, x_range, u)
    """
    assert [len(curve) == 3 for curve in curves], \
        "Incorrect curves argument: {0}".format(curves)

    if match_first_h:
        for index,curve_tuple in enumerate(curves):
            name, x_range, u = curve_tuple
            first_h = curves[0][2][-1]
            h_diff = first_h - u[-1]
            curves[index] = (name, x_range, [y + h_diff for y in u])

    p = figure(plot_width=800, plot_height=400, x_axis_type="linear")
    for name, x_range, u in curves:
        r, g, b = np.random.randint(0, high=255, size=3)
        if 'diff eq' in name:
            p.circle(x_range, u, color=(r, g, b), legend=name)
        else:
            p.line(x_range, u, color=(r, g, b), legend=name, line_width=3)
    p.legend.location = "bottom_right"
    show(p)

# def l(T0=1, g=1, rho=1, a=1):
#     return (T0 / (g*rho))*np.sinh((g*rho*a)/T0)

# T0 = inversefunc(l, )
parameters = {
    'a': 1.0,
    'h': 1.0,
    'rho': 1.0,
    'R': 1.0,
    'g': 9.8,
    'T0': 10.0,
}

x_range, u_c = catenary(**parameters)
x_range, u_r = roadway(**parameters)
lol = [u_c[i] + u_r[i] for i in range(len(u_c))]
curves = [
    ('catenary diff eq',) + solver(case='catenary', **parameters),
    ('roadway diff eq',) +  solver(case='roadway', **parameters),
    ('suspension bridge',) + solver(case='bridge', **parameters),
    ('suspension bridge lol',) + (x_range, lol),
    ('catenary closed form',) + catenary(**parameters),
    ('roadway closed form',) + roadway(**parameters)
]

plot(curves, match_first_h=True)
# h_diff = solution_dict['bridge'][-1] - solution_dict['catenary'][-1]
# solution_dict['catenary'] = [y + h_diff for y in solution_dict['catenary']]