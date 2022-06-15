#!/usr/bin/env python
''' Quick plot for lfric_atm global output '''

# Need to set a non-interactive backend for suites
from __future__ import absolute_import
from __future__ import print_function
import matplotlib
matplotlib.use('Agg')

# Note non-PEP8 collecting of imports as the backend needs to be
# set before we import iris.
import iris
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Index to the fields
varname=0
colbar_min=1
colbar_max=2

# Fields which are available to plot
theta           = ['theta',           265,  300]
grid_surface_temperature = ['grid_surface_temperature',265,300]
m_v             = ['m_v',             1e-3, 15e-3]
m_cl            = ['m_cl',            0,    5e-5]
m_ci            = ['m_ci',            0,    1e-4]
u_in_w3         = ['u_in_w3',        -15,  15]
sw_temperature_incr = ['sw_temperature_incr', 0, 0.1]
lw_temperature_incr = ['lw_temperature_incr', -2, 1]
sw_heating_rate = ['sw_heating_rate', 0,    7e-5]
trop_level      = ['trop_level', 20, 50]
grid_snow_mass  = ['grid_snow_mass',0,1000]
total_prec      = ['total_prec',0,1e-3]
ls_prec         = ['ls_prec',0,1e-3]
sw_direct_toa   = ['sw_direct_toa', 0, 1450]
sw_up_toa       = ['sw_up_toa',0,600]
sw_down_surf    = ['sw_down_surf', 0, 1400]
lw_down_surf    = ['lw_down_surf', 100, 600]
lw_up_surf      = ['lw_up_surf', 100, 600]
lw_up_toa       = ['lw_up_toa', 100, 600]
cloud_amount_maxrnd = ['cloud_amount_maxrnd',0,1]
cloud_cover_rts = ['cloud_cover_rts', 0, 1]
cloud_fraction_rts = ['cloud_fraction_rts', 0, 1]
cloud_droplet_re_rts = ['cloud_droplet_re_rts', 0, 20e-6]
warm_cloud_top_re_rts = ['warm_cloud_top_re_rts', 0, 20e-6]
warm_cloud_top_weight_rts = ['warm_cloud_top_weight_rts', 0, 1]
liq_cloud_frac_rts = ['liq_cloud_frac_rts', 0, 1]
ice_cloud_frac_rts = ['ice_cloud_frac_rts', 0, 1]
liq_cloud_mmr_rts  = ['liq_cloud_mmr_rts', 0, 1e-4]
ice_cloud_mmr_rts  = ['ice_cloud_mmr_rts', 0, 1e-4]
liq_cloud_path_rts  = ['liq_cloud_path_rts', 0, 5e-1]
ice_cloud_path_rts  = ['ice_cloud_path_rts', 0, 5e-1]
sw_aod_rts = ['sw_aer_optical_depth_rts', 0, 1e-3]
sw_net_surf_rts = ['sw_net_surf_rts', 0, 1200]
sw_direct_toa_rts = ['sw_direct_toa_rts', 0, 1450]
sw_up_toa_rts = ['sw_up_toa_rts', 0, 600]
sw_up_rts = ['sw_up_rts', 0, 600]
sw_up_clear_toa_rts = ['sw_up_clear_toa_rts', 0, 600]
sw_down_clear_surf_rts = ['sw_down_clear_surf_rts', 0, 1400]
sw_up_clear_surf_rts = ['sw_up_clear_surf_rts', 0, 600]
lw_net_surf_rts = ['lw_net_surf_rts', -300, 100]
lw_up_toa_rts = ['lw_up_toa_rts', 100, 400]
lw_up_clear_toa_rts = ['lw_up_clear_toa_rts', 100, 400]
lw_down_clear_surf_rts = ['lw_down_clear_surf_rts', 100, 600]
lw_up_clear_surf_rts = ['lw_up_clear_surf_rts', 100, 600]

def load_cube_by_varname(filename, var):
   variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == var))
   return iris.load_cube(filename, constraint=variable_constraint)

def gen_markersize(n_points):
    """
    The "points" based plot of the LFRic grid needs the size of the
    dots to be set sensibly or plotting artifacts are introduced.
    This method uses a rough exponential fit to pick a size.
    """
    c_value = int((n_points/6)**0.5)

    def exponential(x, a, b, c):
        return a*np.exp(-b*x)+c

    # C Value sample points
    x = [48, 108, 168, 192]
    # Markersize values
    y = [40, 12, 6, 1]

    # Create the fit
    pcoeffs, _ = curve_fit(exponential, x, y, p0=(1, 1e-6, 1))

    # Now derive the value (don't return a value less than 1)
    return max(round(exponential(c_value, *pcoeffs)), 1)

def do_plot(datapath, plotfield, plotpath='.', plotlevel=0):
    ''' Do the plotting using data from datapath. Send output to plotpath '''

    lfric = load_cube_by_varname(datapath, plotfield[varname])
    if lfric.ndim == 2:
        lfric = lfric[-1]
    else:
        lfric = lfric[-1, plotlevel]

    # Get the x and y co-ordinates
    x_coord = np.around(lfric.coord('longitude').points, decimals=5)
    y_coord = np.around(lfric.coord('latitude').points,  decimals=5)

    # Save the min and max of the data
    field_min = np.around(np.min(lfric.data), decimals=7)
    field_max = np.around(np.max(lfric.data), decimals=7)

    # Set up the colourbar
    plt.set_cmap(plt.cm.RdYlBu_r)

    plt.figure(figsize=(8, 5))
    markersize = gen_markersize(lfric.data.shape[0])
    plot = plt.scatter(x_coord, y_coord, c = lfric.data,
                edgecolor = "none", s = markersize,
                vmin = plotfield[colbar_min],
                vmax = plotfield[colbar_max])
    plt.colorbar(plot,orientation='vertical')

    plt.title(plotfield[varname]+', min = '+str(field_min)
                                +', max = '+str(field_max) )
    plt.xlim([np.min(x_coord), np.max(x_coord)])
    plt.ylim([-90, 90])

    plt.savefig(plotpath+'/'+plotfield[varname]+'.png', bbox_inches='tight')


if __name__ == "__main__":

    import sys
    try:
        opts = [opt for opt in sys.argv[1:] if opt.startswith('-')]
        args = [arg for arg in sys.argv[1:] if not arg.startswith('-')]
        datapath, plotpath = args[0:2]
        rts_plots = '-rts' in opts
        ral_plots = '-ral' in opts
    except ValueError:
        print("Usage: {0} <datapath> <plotpath>".format(sys.argv[0]))
        exit(1)
    do_plot(datapath, theta,           plotpath)
    do_plot(datapath, m_v,             plotpath)
    do_plot(datapath, m_cl,            plotpath)
    do_plot(datapath, m_ci,            plotpath)
    do_plot(datapath, u_in_w3,         plotpath)
    do_plot(datapath, cloud_amount_maxrnd, plotpath)
    do_plot(datapath, grid_surface_temperature, plotpath)
    do_plot(datapath, grid_snow_mass,  plotpath)
    do_plot(datapath, sw_up_toa,       plotpath)
    do_plot(datapath, sw_down_surf,    plotpath)
    do_plot(datapath, lw_down_surf,    plotpath)
    do_plot(datapath, lw_up_toa,       plotpath)
    do_plot(datapath, trop_level,      plotpath)
    if rts_plots:
        do_plot(datapath, cloud_cover_rts,        plotpath)
        do_plot(datapath, cloud_fraction_rts,     plotpath, plotlevel=17)
        do_plot(datapath, cloud_droplet_re_rts,   plotpath, plotlevel=17)
        do_plot(datapath, sw_aod_rts,             plotpath, plotlevel=38)
        do_plot(datapath, sw_net_surf_rts,        plotpath)
        do_plot(datapath, sw_direct_toa_rts,      plotpath)
        do_plot(datapath, sw_up_toa_rts,          plotpath)
        do_plot(datapath, sw_up_rts,              plotpath, plotlevel=70)
        do_plot(datapath, sw_up_clear_toa_rts,    plotpath)
        do_plot(datapath, sw_down_clear_surf_rts, plotpath)
        do_plot(datapath, sw_up_clear_surf_rts,   plotpath)
        do_plot(datapath, lw_net_surf_rts,        plotpath)
        do_plot(datapath, lw_up_toa_rts,          plotpath)
        do_plot(datapath, lw_up_clear_toa_rts,    plotpath)
        do_plot(datapath, lw_down_clear_surf_rts, plotpath)
        do_plot(datapath, lw_up_clear_surf_rts,   plotpath)
    if ral_plots:
        do_plot(datapath, ls_prec,      plotpath)
    else:
        do_plot(datapath, total_prec,   plotpath)
