#!/usr/bin/env python2.7
''' Quick plot of for lfric_atm scm output '''

# Need to set a non-interactive backend for suites
import matplotlib
matplotlib.use('Agg')

# Note non-PEP8 collecting of imports as the backend needs to be
# set before we import iris.
import iris
import matplotlib.pyplot as plt


def do_plot(datapath, plotpath='.'):
    ''' Do the plotting using data from datapath. Send output to plotpath '''

    lfric = iris.load_cube(datapath,
                           iris.Constraint(name='air_potential_temperature'))
    lfric = lfric[:, :, 0]

    plt.figure(figsize=(15, 10))
    for n, time in enumerate([0, 10, 20, 30, 40, 47]):
        plt.subplot(2, 3, n+1)
        plt.plot(lfric.data[time, 1:],
                 lfric.coord('full_levels_faces').points[1:],
                 label='LFRic',
                 linewidth=2)

        plt.legend(loc='best')
        plt.xlabel('theta')
        plt.ylabel('Model Level Number')
        plt.xlim([290, 315])

        plt.title('Timestep = '+str(time+1))

    plt.savefig(plotpath+'/lfric_scm_theta.png', bbox_inches='tight')


if __name__ == "__main__":

    import sys
    try:
        datapath, plotpath = sys.argv[1:3]
    except ValueError:
        print "Usage: {0} <datapath> <plotpath>".format(sys.argv[0])
        exit(1)
    do_plot(datapath, plotpath)
