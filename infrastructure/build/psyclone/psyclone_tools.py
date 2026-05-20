# -----------------------------------------------------------------------------
#  (C) Crown copyright 2023 Met Office. All rights reserved.
#  The file LICENCE, distributed with this code, contains details of the terms
#  under which the code may be used.
# -----------------------------------------------------------------------------


"""
This file contains frequently used transformations to simplify
their application in PSyclone optimisations scripts.

"""

from psyclone.domain.lfric import LFRicConstants
from psyclone.psyGen import InvokeSchedule
from psyclone.psyir.nodes import Loop, Routine, Directive
try:
    from psyclone.psyir.transformations import OMPParallelTrans
except ImportError:
    # Support for psyclone < 3.3
    from psyclone.transformations import OMPParallelTrans
try:
    from psyclone.domain.lfric.transformations import (
        LFRicRedundantComputationTrans)
except ImportError:
    # Support for psyclone < 3.3
    from psyclone.transformations import LFRicRedundantComputationTrans

from psyclone.transformations import (
    LFRicColourTrans,
    LFRicOMPLoopTrans,
)

# List of allowed 'setval_*' built-ins for redundant computation transformation
SETVAL_BUILTINS = ["setval_c"]


# -----------------------------------------------------------------------------
def redundant_computation_setval(psyir):
    """
    Applies the redundant computation transformation to loops over DoFs
    for the initialision built-ins, 'setval_*'.

    To reduce MPI communications, current PSyclone-LFRic strategy does not
    apply halo swaps on input arguments to kernels with increment
    operations on continuous fields such as 'GH_INC'. For such kernels,
    PSy-layer code needs to loop into the halo to correctly compute owned
    DoFs on the boundary between the halo and the domain. Therefore values
    of the remaining DoFs in the first halo cell need to be initialised to
    values that will not induce numerical errors.

    By default, the initialisation 'setval_*' built-ins do not initialise
    into the halos. This transform causes them to do so, and so permits
    developers to set safe values in halos.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    :raises Exception: if there is more than one built-in call per DoF loop.

    """
    # Import redundant computation transformation
    rtrans = LFRicRedundantComputationTrans()

    # Loop over all the InvokeSchedule in the PSyIR object
    for subroutine in psyir.walk(InvokeSchedule):
        # Make setval_* built-ins compute redundantly to the level-1 halo
        # if they are in their own loop
        for loop in subroutine.loops():
            if loop.iteration_space == "dof":
                if len(loop.kernels()) != 1:
                    raise Exception(
                        f"Expecting loop to contain 1 call but found "
                        f"'{len(loop.kernels())}'"
                    )
                if loop.kernels()[0].name in SETVAL_BUILTINS:
                    rtrans.apply(loop, options={"depth": 1})


# -----------------------------------------------------------------------------
def colour_loops(psyir, enable_tiling=False):
    """
    Applies the colouring transformation to all applicable loops and optionally
    enables tiling.
    It creates the instance of `LFRicColourTrans` only once.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    """
    const = LFRicConstants()
    ctrans = LFRicColourTrans()

    # Loop over all the subroutines in the PSyIR object
    for subroutine in psyir.walk(Routine):
        # Colour loops over cells unless they are on discontinuous
        # spaces or over DoFs
        for child in subroutine.children:
            if (
                isinstance(child, Loop)
                and child.iteration_space.endswith("cell_column")
                and child.field_space.orig_name
                not in const.VALID_DISCONTINUOUS_NAMES
            ):
                ctrans.apply(child, options={"tiling": enable_tiling})


# -----------------------------------------------------------------------------
def openmp_parallelise_loops(psyir):
    """
    Applies OpenMP Loop transformation to each applicable loop.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    """
    otrans = LFRicOMPLoopTrans()
    oregtrans = OMPParallelTrans()

    # Loop over all the InvokeSchedule in the PSyIR object
    for subroutine in psyir.walk(InvokeSchedule):
        # Add OpenMP to loops unless they are over colours, are null,
        # or if an outer loop is already parallelised (OpenMP is applied
        # to loop over tiles instead of cells if tiling is enabled)
        for loop in subroutine.loops():
            if loop.loop_type not in ["colours", "null"] and \
               not loop.ancestor(Directive):
                oregtrans.apply(loop)
                otrans.apply(loop, options={"reprod": True})


# -----------------------------------------------------------------------------
def view_transformed_schedule(psyir):
    """
    Provides view of transformed Invoke schedule in the PSy-layer.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    """
    setval_count = 0

    # Loop over all the Invokes in the PSyIR object
    for subroutine in psyir.walk(InvokeSchedule):
        print(f"Transformed invoke '{subroutine.name}' ...")

        # Count instances of setval_* built-ins
        for loop in subroutine.loops():
            if loop.iteration_space == "dof":
                if loop.kernels()[0].name in SETVAL_BUILTINS:
                    setval_count += 1

        # Take a look at what we have done
        print(f"Found {setval_count} {SETVAL_BUILTINS} calls")
        print(subroutine.view())
