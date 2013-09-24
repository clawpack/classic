def setplot(plotdata):
    plotdata.clearfigures()

    # Figures corresponding to Figure 9.5 of LeVeque, "Finite Volume
    # Methods for Hyperbolic Problems," 2002 (though more of them)

    # Tuples of (variable name, variable number)
    figdata = [('Pressure', 0),
               ('Velocity', 1)]

    # Afteraxes function: draw a vertical dashed line at the interface
    # between different media
    def draw_interface(current_data):
        import pylab
        pylab.plot([0., 0.], [-1000., 1000.], 'k--')

    for varname, varid in figdata:
        plotfigure = plotdata.new_plotfigure(name=varname, figno=varid)

        plotaxes = plotfigure.new_plotaxes()
        plotaxes.xlimits = [-5., 5.]
        plotaxes.ylimits = [-0.5, 1.5]    # Good for both vars because of near-unit impedance
        plotaxes.title = varname
        plotaxes.afteraxes = draw_interface

        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.plot_var = varid
        plotitem.color = 'b'

    plotdata.printfigs = True          # Whether to output figures
    plotdata.print_format = 'png'      # What type of output format
    plotdata.print_framenos = 'all'    # Which frames to output
    plotdata.print_fignos = 'all'      # Which figures to print
    plotdata.html = True               # Whether to create HTML files
    plotdata.latex = False             # Whether to make LaTeX output

    return plotdata

