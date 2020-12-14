
# https://www.kdnuggets.com/2016/11/linear-regression-least-squares-matrix-multiplication-concise-technical-overview.html
# https://stackoverflow.com/questions/53698635/how-to-define-a-plane-with-3-points-and-plot-it-in-3d
# https://math.stackexchange.com/questions/99299/best-fitting-plane-given-a-set-of-points

import numpy as np
import matplotlib
# matplotlib.use('TkAgg')
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


def calc_plane(X_val, Y_val, Z_val):
    # Setup matrices
    Npoints = np.shape(X_val)[0]

    # The model is: Z_val = p[0] + p[1]*X_val + p[2]*Y_val
    Npar = 3
    X = np.zeros((Npoints, Npar))
    y = np.zeros((Npoints, 1))
    X[:, 0] = np.ones(Npoints)
    X[:, 1] = X_val
    X[:, 2] = Y_val
    y[:, 0] = Z_val

    # Solve for projection matrix
    par = np.linalg.inv(X.T.dot(X)).dot(X.T).dot(y)
    residuals = y - X.dot(par)
    SSE = np.linalg.norm(residuals)

    print("solution: z = %f + %f x + %f y" % (par[0], par[1], par[2]))
    print("residuals:")
    print(residuals)
    print("SSE: %f" % SSE)

    return par


def plot_plane(x_val, y_val, z_val, par, zstring, zmin, zmax):
    # plot raw data
    fig = plt.figure(figsize=(8, 6))
    # ax = fig.add_subplot(111, projection="3d")
    ax = Axes3D(fig)
    angle = 45
    ax.view_init(45, angle)


    # ax.scatter(x_val, y_val, z_val, s=100, color='r')
    ax.plot3D(x_val, y_val, z_val, 'o-', markersize=10, color='r', zorder=2)
    mmx = get_fix_mins_maxs(-10, 10)
    ax.set_xlim3d(mmx)
    mmy = get_fix_mins_maxs(20, 60)
    ax.set_ylim3d(mmy)
    mmz = get_fix_mins_maxs(zmin, zmax)
    ax.set_zlim3d(mmz)

    plt.title('Bilinear fit to {0} according to NTA8800 Appendix Q'.format(zstring))
    ax.set_xlabel('Tin')
    ax.set_ylabel('Tout')
    ax.set_zlabel(zstring)

    # plot plane
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # Find regression line
    xx, yy = np.meshgrid(np.arange(xlim[0], xlim[1], 2.0),
                         np.arange(ylim[0], ylim[1], 2.0))
    zz = np.array(par[0] + par[1] * xx + par[2] * yy)

    # zz = np.zeros(xx.shape)
    # for row in range(xx.shape[0]):
      #  for col in range(xx.shape[1]):
       #     zz[row, col] = par[0] + par[1] * xx[row, col] + par[2] * yy[row, col]

    # plot the surface
    # ax.plot_wireframe(xx, yy, zz, color='lime', alpha=0.5)
    ax.plot_surface(xx, yy, zz, color='lime', alpha=0.5)
    CS1 = ax.contour(xx, yy, zz, zdir = 'z', cmap=cm.winter)

    plt.show()


def get_fix_mins_maxs(mins, maxs):
    deltas = (maxs - mins) / 12.
    mins = mins + deltas / 4.
    maxs = maxs - deltas / 4.

    return [mins, maxs]


if __name__ == "__main__":
    Tin = np.array([7, 7, -7])
    Tout = np.array([35, 55, 35])
    COP = np.array([4.0, 3.0, 2.5])
    coef = calc_plane(Tin, Tout, COP)
    plot_plane(Tin, Tout, COP, coef, 'COP', 1.0, 5.0)

    P_max = np.array([6.0, 2.0, 3.0])
    P_coef = calc_plane(Tin, Tout, P_max)
    plot_plane(Tin, Tout, P_max, P_coef, 'Power', 0.0, 10.0)


