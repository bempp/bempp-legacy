def plot_slice(plot_f, dolfin_function=None, x=None, y=None, z=None, n=(151,151), extent=(-5,5,-5,5),
                 cmap='coolwarm', filename=None, title=None, show=None):
    """Plots a slice of the function plot_f at a constant x, y or z.

    Parameters
    ----------
    plot_f : function
        A Python function which takes one input, points, and returns
        the value at these points.

    dolfin_function : function
        A Python function which takes one input, point, and returns the
        value at that point, or None if the point is outside the FEM
        domain.

    x : float
        The x-coordinate to plot a slice at.

    y : float
        The y-coordinate to plot a slice at.

    z : float
        The z-coordinate to plot a slice at.

    n : tuple
        The number of points to use on each axis.

    extent : tuple
        The bounds of the two axes.

    cmap : string
        The colour map to use on the plot.

    filename : string
        The file name to save the plot to.

    title : string
        The title of the plot.

    show : bool
        If True, the plot will be displayed. If False, it wil not.

    """
    import numpy as np
    if show is None:
        show = filename==None

    if x is not None and y is None and z is None:
        x_p, y_p, z_p = np.mgrid[x:x:1j, extent[0]:extent[1]:n[0]*1j, extent[2]:extent[3]:n[1] * 1j]
    elif x is None and y is not None and z is None:
        x_p, y_p, z_p = np.mgrid[extent[0]:extent[1]:n[0]*1j, y:y:1j, extent[2]:extent[3]:n[1] * 1j]
    elif x is None and y is None and z is not None:
        x_p, y_p, z_p = np.mgrid[extent[0]:extent[1]:n[0]*1j, extent[2]:extent[3]:n[1] * 1j, z:z:1j]
    else:
        raise TypeError("Exactly one of x, y and z must be set.")
    points = np.vstack((x_p.ravel(), y_p.ravel(), z_p.ravel()))


    if dolfin_function is None:
        plot_me = np.real(plot_f(points))
    else:
        plot_me = np.zeros(points.shape[1], dtype=np.float64)
        bem_x = np.full(points.shape[1], False, dtype=np.bool)
        fem_plot = []
        for i,p in enumerate(points.T):
            value = dolfin_function(np.array(p))
            if value is None:
                bem_x[i] = True
            else:
                fem_plot.append(value)
        plot_me[np.logical_not(bem_x)] = np.real(np.array(fem_plot))
        plot_me[bem_x] = np.real(plot_f(points[:,bem_x]))

    import matplotlib
    from matplotlib import pyplot as plt

    plt.imshow(plot_me.reshape(n).T,
               cmap=cmap, origin='lower',
               extent=extent)
    plt.colorbar()
    if title is None:
        if x is not None:
            title = "Plot at x="+str(x)
        if y is not None:
            title = "Plot at y="+str(y)
        if z is not None:
            title = "Plot at z="+str(z)
    plt.title(title)
    if x is None:
        plt.xlabel("x")
        if y is None:
            plt.ylabel("y")
        else:
            plt.ylabel("z")
    else:
        plt.xlabel("y")
        plt.ylabel("z")
    if filename is not None:
        plt.savefig(filename)
    if show:
        plt.show()

    plt.clf()
