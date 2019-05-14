""" Differential equations for various cases """
from bokeh.io import export_png
from bokeh.plotting import figure, output_file, show
import math
import numpy as np
from pynverse import inversefunc
from scipy import integrate

""" Differential Equations """
def bridge_derivative(Y, x, rho=1, R=1, g=9.8, T0=1, **kwargs):
    """
    Returns Y' = [y2, y2'] where Y = [y1, y2], y1 = u(x), y2 = u'(x)
    Both the mass of the chain and roadway are considered.
    """
    y1, y2 = Y
    y2_prime = (rho*np.sqrt(1 + y2**2) + R) * (g / T0)
    return [y2, y2_prime]

def catenary_derivative(Y, x, rho, R, g, T0, **kwargs):
    """
    Returns Y' = [y2, y2'] where Y = [y1, y2], y1 = u(x), y2 = u'(x)
    Only rho, the mass of the chain, is considered.
    """
    y1, y2 = Y
    y2_prime = ((g * rho) / T0) * np.sqrt(1 + y2**2)
    return [y2, y2_prime]


def roadway_derivative(Y, x, rho, R, g, T0, **kwargs):
    """
    Returns Y' = [y2, y2'] where Y = [y1, y2], y1 = u(x), y2 = u'(x)
    Only R, the mass of the roadway, is considered.
    """
    y1, y2 = Y
    y2_prime =  R * (g / T0)
    return [y2, y2_prime]


def solver(x_range, a=1.0, h=1.0, rho=1.0, R=1.0, g=9.8, T0=9.8, case=None, n=1000.0, l=0.0):
    """
    Solve for the first derivative of u for the given case.

    case must be one of {'catenary', 'roadway', 'bridge'}
    """
    derivatives = {
        'catenary': catenary_derivative,
        'roadway': roadway_derivative,
        'bridge': bridge_derivative,
    }
    assert case in derivatives.keys(), "Please provide a valid case"

    # if not x_range:
    #     x_range = np.arange(0, a, 0.01)
    U_x = lambda Y, x: derivatives[case](Y, x, rho, R, g, T0, a=a, n=n)
    solution = integrate.odeint(U_x, [0, 0], x_range)
    return (x_range, solution[:, 0])


""" Closed-form equations """
def catenary(a=1.0, h=1.0, rho=1.0, R=1.0, g=9.8, T0=9.8, l=2.0):
    """
    Return the closed form of the catenary (hyperbolic cosine)
    """
    print('catenary T0 ' + str(T0))
    C = (g * rho) / T0
    x_range = np.arange(0, a, 0.01)

    result_tuple = (x_range, [(1/C) * (np.cosh(C*x) - np.cosh(C*a)) + h for x in x_range])

    return result_tuple


def roadway(a=1.0, h=1.0, rho=1.0, R=1.0, g=9.8, T0=9.8, l=2.0):
    """
    Return the closed form solution of the negligible chain mass (x^2)
    """
    print('roadway T0 ' + str(T0))
    x_range = np.arange(0, a, 0.01)
    return (x_range, [((g*R)/(2*T0))*(x**2 - a**2) + h for x in x_range])


""" Plotting utilities"""
def plot(curves, match_first_h=False, filename=None):
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
            # h_diff is the difference between the values of the curves at a
            # (the last point)
            first_h = curves[0][2][-1]
            h_diff = first_h - u[-1]
            curves[index] = (name, x_range, [y + h_diff for y in u])

    p = figure(plot_width=800, plot_height=400, x_axis_type="linear")
    for name, x_range, u in curves:
        r, g, b = np.random.randint(0, high=255, size=3)
        if 'Bridge' in name:
            p.diamond(x_range, u, color='black', legend=name)
        elif 'ODE' in name:
            p.circle(x_range, u, color='black', legend=name)
        elif 'Catenary' in name:
            p.line(x_range, u, color='red', legend=name, line_width=3)
        elif 'Negligible' in name:
            p.line(x_range, u, color='blue', legend=name, line_width=3)
    p.legend.location = "bottom_right"
    show(p)

""" Functions to calculate the horizontal tension T0 """
def T0_catenary(l=2.0, g=9.8, rho=1.0, a=1.0, **kwargs):
    l_fn = lambda T0: (T0 / (g*rho)) * np.sinh(((g*rho) / T0)*a)
    result = inversefunc(l_fn, y_values=l, domain=[0, 10**9], image=[0, 10**9],
                       open_domain=True)
    print(result)
    return result


def T0_roadway(l=2.0, g=9.8, R=1.0, a=1.0, **kwargs):
    integrand = lambda x,T0: np.sqrt(1.0 + ((g*R*x)/T0)**2)
    l_fn = lambda T0: integrate.quad(integrand, 0, a, args=(T0,))[0]
    result = inversefunc(l_fn, y_values=l, domain=[0, 10**9], image=[0, 10**9],
                         open_domain=True)
    print(result)
    return result


def B(v, T0, g, R, rho):
    r = (R / rho)
    result = (r/np.sqrt(r**2 - 1))*(-1.0*np.arctanh((r*v)/(np.sqrt(r**2 -1)*np.sqrt(v**2 + 1))))
    result += (r/np.sqrt(r**2 - 1))*np.arctanh(v / np.sqrt(r**2 - 1))
    result +=  np.arcsinh(v)
    return (T0/(g*rho))*result


def T0_bridge(l=2.0, g=9.8, rho=1.0, R=1.0, a=1.0, **kwargs):
    integrand = lambda x,T0: np.sqrt(1.0 + inversefunc(B, args=(T0, g, R, rho),
                                     y_values=x, domain=[0, 5], image=[0, 10],
                                     open_domain=False)**2)
    l_fn = lambda T0: integrate.quad(integrand, 0, a, args=(T0,))[0]
    result = inversefunc(l_fn, y_values=l, domain=[0.0, 10**3], image=[0, 10**3],
                   open_domain=True)
    print('T0 bridge: ' + str(result))
    return result

