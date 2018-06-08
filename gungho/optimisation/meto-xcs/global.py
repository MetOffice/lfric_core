##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################


''' PSyclone transformation script for the Dynamo0p3 API
to apply colouring and OpenMP generically. '''
from psyclone.transformations import Dynamo0p3ColourTrans, \
                                     Dynamo0p3OMPLoopTrans, \
                                     OMPParallelTrans
from psyclone.dynamo0p3 import DISCONTINUOUS_FUNCTION_SPACES


def trans(psy):
    ''' Applies PSyclone colouring and OpenMP transformations. '''
    ctrans = Dynamo0p3ColourTrans()
    otrans = Dynamo0p3OMPLoopTrans()
    oregtrans = OMPParallelTrans()

    # Loop over all of the Invokes in the PSy object
    for invoke in psy.invokes.invoke_list:

        print "Transforming invoke '{0}' ...".format(invoke.name)
        schedule = invoke.schedule

        # Colour loops over cells unless they are on discontinuous
        # spaces (W3, WTHETA and W2V) or over dofs
        for loop in schedule.loops():
            if loop.iteration_space == "cells" \
                and loop.field_space.orig_name \
                    not in DISCONTINUOUS_FUNCTION_SPACES:
                schedule, _ = ctrans.apply(loop)

        # Add OpenMP to loops unless they are over colours
        for loop in schedule.loops():
            if loop.loop_type != "colours":
                schedule, _ = oregtrans.apply(loop)
                schedule, _ = otrans.apply(loop, reprod=True)

        # Take a look at what we've done
        schedule.view()

    return psy
