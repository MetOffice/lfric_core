#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
'''
Basic python script to plot the x-z profile of variables for the orography test cases
on subplots.

This version takes nodal format output files and
interpolates onto a regular grid.

Levels are determined from the data

This version stitches together a directory of files
and extracts all levels so it can work in the serial
case where there is one file or the parallel case where
there is a file for each processor.

'''


import numpy as np

# Need to set a non-interactive backend for suites
import matplotlib
matplotlib.use('Agg')

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import sys

from read_data import read_nodal_data

levels = None
levels0 = None
data = None
data0 = None

def make_figure(plotpath, nx, ny, field, component, timestep, plotlevel_x, plotlevel_y, case, zoom, cntrs, sp_cols, iii, slice_fig):

  val_col = 'c' + str(component)

  # get min and max of x,y data for plot axes

  min_lev = min(levels)

  xmin = data.loc[data['level'] == min_lev]['x'].min()
  xmax = data.loc[data['level'] == min_lev]['x'].max()
  ymin = data.loc[data['level'] == min_lev]['y'].min()
  ymax = data.loc[data['level'] == min_lev]['y'].max()
  
  # zmin is min of bottom level

  zmin = data['z'].min()

  # zmax is a test-dependent value
  if zoom == 'zoom_1': # zooming as in paper
    if case != ('cosine' and 'hat'):
      zmax = 12000.0
    else:
      print 'Error: cosine and hat orography do not have the zoom option'
  elif zoom == 'zoom_2': # zooming up to bottom of damping layer
    if case == 'nhmw':
      zmax = 25000.0
    elif case == 'hmw':
      zmax = 30000.0
    elif case == 'schar':
      zmax = 20000.0
    else:
      print 'Error: cosine and hat orography do not have the zoom option'
  elif zoom == 'no_zoom': # plotting the entire domain depth
    if case == 'nhmw':
      zmax = 35000.0
    elif case == 'hmw':
      zmax = 50000.0
    elif case == 'schar':
      zmax = 30000.0
    elif case == 'cosine':
      zmax = 6000.0
    elif case == 'hat':
      zmax = 2000.0


  val_col = 'c' + str(component)

  nx = int(nx)
  ny = int(ny)
  nl = len(levels)

 
  #create 2D plot
  val_i = np.zeros([ny,nx,nl])
  val_i0 = np.zeros([ny,nx,nl])
  x_i  = np.zeros([ny,nx,nl])
  height_i  = np.zeros([ny,nx,nl])

  # interpolate field values and heights onto xy for each level 
  for p in xrange(len(levels)):
    p_data = data.loc[data['level'] == levels[p]]

    # Using reshape of numpy array
    val_i[:,:,p] = (p_data[val_col].values).reshape((ny, nx))
    height_i[:,:,p] = (p_data['z'].values).reshape((ny, nx))
    x_i[:,:,p]  = (p_data['x'].values).reshape((ny, nx))
    
  if field in ('theta', 'm_c', 'rho'):
    for p in xrange(len(levels0)):
      p_data = data0.loc[data0['level'] == levels0[p]] 
      # Using reshape of numpy array
      val_i0[:,:,p] = (p_data[val_col].values).reshape((ny, nx))
      

  if field in ('theta', 'm_c', 'rho'): # subtracting entire background profile (profile at initial time) if plotting theta
    val_i -= val_i0

  # Take actual orography heights as z
  zi_adj = height_i[0,:,:]


  # extract the slice for plotting
  dval = np.zeros([nx,len(levels)])
  dval = val_i[1,:,:] 

  x_plt = x_i[0,:,:] 

  ax1 = slice_fig.add_subplot(1,sp_cols,iii+1)
  
  # Setting contour limits and intervals for vertical velocity
  
  if case == 'nhmw':
    cmin_w = -0.0048
    cmax_w = 0.0048
    nc_w   = 17
    cmin_t = -0.002
    cmax_t = 0.002
    nc_t = 9
  elif case == 'hmw':
    cmin_w = -0.004
    cmax_w = 0.004
    nc_w   = 17
  elif case == 'schar':
    cmin_w = -0.4
    cmax_w = 0.4
    nc_w   = 17
  elif case == 'cosine':
    cmin_w = -3.0 
    cmax_w = 3.0
    nc_w   = 13
  elif case == 'hat':
    cmin_t = -0.35
    cmax_t = 0.35
    nc_t   = 15
      
  # Plots 
  if (field=='w3projection_u1' or (field=='u' and component=='1')):
    if case == 'nhmw':
      if cntrs == 'lines':
        cc = np.linspace(9.99, 10.01, 21)
        cf = ax1.contour(x_plt, zi_adj, np.round(dval, 10), cc, cmap=cm.spectral, linewidths=3)
      elif cntrs == 'colours': 
        cf = ax1.pcolormesh(x_plt, zi_adj, np.round(dval, 10), vmin=9.99, vmax=10.01, cmap=cm.coolwarm)
    elif case == 'cosine':
      if cntrs == 'lines':
        cc = np.linspace(7.0, 12.6, 15)
        cf = ax1.contour(x_plt, zi_adj, np.round(dval, 10), cc, colors='k', linewidths=3)
      elif cntrs == 'colours': 
        cf = ax1.pcolormesh(x_plt, zi_adj, np.round(dval, 10), vmin=7.0, vmax=12.6, cmap=cm.coolwarm)
    else:
      cf = ax1.pcolormesh(x_plt, zi_adj, np.round(dval, 10))
      plt.colorbar(cf,  cmap=cm.coolwarm)
  elif (field=='w3projection_u3' or (field=='u' and component=='3')):
    if cntrs == 'lines':
      cc = np.linspace(cmin_w, cmax_w, nc_w)
      if case == 'cosine':
        cf = ax1.contour(x_plt, zi_adj, np.round(dval, 10), cc, colors='k', linewidths=3)
      else:
        cf = ax1.contour(x_plt, zi_adj, np.round(dval, 10), cc, cmap=cm.spectral, linewidths=3)
    elif cntrs == 'colours': 
      cf = ax1.pcolormesh(x_plt, zi_adj, np.round(dval, 10), vmin=cmin_w, vmax=cmax_w, cmap=cm.coolwarm)
  elif field=='theta':
    if cntrs == 'lines':
      cc = np.linspace(cmin_t, cmax_t, nc_t)
      cf = ax1.contour(x_plt, zi_adj, np.round(dval, 10), cc, cmap=cm.spectral, linewidths=3)
      plt.colorbar(cf,  cmap=cm.spectral)
    elif cntrs == 'colours': 
      cf = ax1.pcolormesh(x_plt, zi_adj, np.round(dval, 10), vmin=cmin_t, vmax=cmax_t, cmap=cm.coolwarm)
      plt.colorbar(cf,  cmap=cm.coolwarm)
  elif field=='rho':
    if cntrs == 'lines':
      cc = np.linspace(-0.00001, 0.00001, 21)
      cf = ax1.contour(x_plt, zi_adj, np.round(dval, 10), cc, cmap=cm.spectral, linewidths=3)
      plt.colorbar(cf,  cmap=cm.spectral)
    elif cntrs == 'colours': 
      cf = ax1.pcolormesh(x_plt, zi_adj, np.round(dval, 10), vmin=-0.00001, vmax=-0.00001, cmap=cm.coolwarm)
      plt.colorbar(cf,  cmap=cm.coolwarm)
  else:
      if cntrs == 'lines':
        cc = np.linspace(-0.0005, 0.0005, 21)
        cf = ax1.contour(x_plt, zi_adj, np.round(dval, 10), cc, cmap=cm.spectral, linewidths=3)
        plt.colorbar(cf,  cmap=cm.spectral)
      elif cntrs == 'colours':
        cf = ax1.pcolormesh(x_plt, zi_adj, np.round(dval, 10), cmap=cm.coolwarm)
        plt.colorbar(cf,  cmap=cm.coolwarm)

  ax1.set_title('max: %2.4e, min: %2.4e'%(np.max(dval),np.min(dval)))
  ax1.set_xlabel("x(m)", fontsize=24)
  ax1.set_ylabel("z(m)", fontsize=24)

  if zoom == 'no_zoom': # plotting the entire domain horizontal length
    if case == 'nhmw':
      xmin_lim = -72000
      xmax_lim  = 72000
    elif case == 'hmw':
      xmin = -200000
      xmax_lim = 200000
    elif case == 'schar':
      xmin_lim = -50000
      xmax_lim = 50000
    elif case == 'hat':
      xmin_lim = -2000
      xmax_lim = 2000
      ax1.set_xticks(np.arange(-2000, 2100, 100))
    else: # cosine case
      xmin_lim = -2000
      xmax_lim = 2000
  else: # plotting only a portion of horizontal domain length
    if case == 'nhmw':
      xmin_lim = -12500
      xmax_lim = 32500
      ax1.set_xticks(np.arange(-10000, 35000, 5000))
    elif case == 'hmw':
      xmin_lim = -40000
      xmax_lim = 40000
      ax1.set_xticks(np.arange(-40000, 45000, 10000))
    elif case == 'schar':
      xmin_lim = -20000
      xmax_lim = 20000
      ax1.set_xticks(np.arange(-20000, 25000, 5000))
    else:
      print 'Error: cosine orography does not have the zoom option.'

  if zoom == 'zoom_1': # Ticks as in paper 
    ax1.set_yticks(np.arange(0, 14000, 2000))
  else:
    if case == 'cosine': 
      ax1.set_yticks(np.arange(0, 7000, 1000))
    elif case == 'hat': 
      ax1.set_yticks(np.arange(0, 2100, 100))

  ax1.set_xlim(xmin_lim, xmax_lim)
  if case == 'hat':
    ax1.tick_params(axis='x', labelsize=8)
    ax1.tick_params(axis='y', labelsize=14)
    ax1.grid(b=True, which='major', color='k', linestyle='--')
    ax1.set_aspect(1.0)
  else:
    ax1.tick_params(axis='both', labelsize=24)
  ax1.set_ylim(zmin, zmax)


# end for loop on fields

if __name__ == "__main__":

  try:
    datapath, nx, ny, fields, component, timesteps, plotlevel_x, plotlevel_y, case, zoom, cntrs, plotpath = sys.argv[1:13]
  except ValueError:
    print("Usage: {0} <datapath> <nx> <ny> <field_names> <component> <timestep_list> <plotlevel_x> <plotlevel_y> <case> <zoom> <cntrs> <plotpath>".format(sys.argv[0]))
    exit(1)


  # Split out the list of fields
  field_list = fields.split(':')

  # Split out the list of timesteps
  ts_list = timesteps.split(':')


  for field in field_list:
    if((field!='u') & (field!='xi')):
      if (component!='1'):
        print "Scalars have only one component!"
        exit(1)
    if field in ('theta', 'm_c', 'rho'):
      # Create initial data for theta
      filestem =  datapath + "/diagDynamo_nodal_" + field + "_T000000" + "*"
      data0 = read_nodal_data(filestem, 1, component)

      # Sort the data (needed to be able to reshape and not regrid)
      data0 = data0.sort(['y','x','z'])

      levels0 = np.sort(data0.level.unique())


  for ts in ts_list:

    # Making space for more fields in the plot. 

    slice_fig = plt.figure(figsize=(15*(len(field_list)),10))

    for i, field in enumerate(field_list):

      filestem =  datapath + "/diagDynamo_nodal_" + field + "_" + ts + "*"

      if field in ['u']:
        data = read_nodal_data(filestem, 3, component)
      else:
        data = read_nodal_data(filestem, 1, component)

      # Sort the data (needed to be able to reshape and not regrid)
      data = data.sort_values(['y','x','z'])

      levels = np.sort(data.level.unique())


      # Only try to plot if we found some files for this timestep
      if len(levels) > 0:
        make_figure(plotpath, nx, ny, field, component, ts, plotlevel_x, plotlevel_y, case, zoom, cntrs, len(field_list), i, slice_fig)

    out_file_name = plotpath + "/" "nodal_slices_" + field + '_' + ts +  ".png"
    plt.tight_layout()
    plt.savefig(out_file_name)

