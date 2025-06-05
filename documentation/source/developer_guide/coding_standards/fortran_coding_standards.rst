.. ------------------------------------------------------------------------------
     (c) Crown copyright 2023 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------

.. _fortran coding standards:

Fortran Coding Standards for LFRic
==================================

Rules for coding standards are driven by the need for readability and
consistency (which aids readability). While some people are happy to read
inconsistent code, other people find inconsistency to be distracting and their
needs should be respected.

Some of the rules are required to meet the technical needs of the LFRic
core and application code structure and organisation.

LFRic coding standards start from the `UM standards
<https://code.metoffice.gov.uk/doc/um/latest/umdp.html#003>`_ which cover
Fortran 95. These rules should be followed unless the LFRic rules on this page
override them.

When can I break the rules
--------------------------

In the following, if the word **must** is used, the rule must not be broken unless
it is completely unavoidable.

If the word must is not used, then the standard *should* be followed unless it
can be argued that breaking the standard is better in the particular context of
the code. Routinely breaking the standard because you prefer a different style
is not a sufficient argument. The code reviewer's judgement is final.

Copyright
---------

The copyright statement references the LFRic licence, and must be included in
all new LFRic code. A Fortran example for 2024 is:

.. code-block:: rst

 !-----------------------------------------------------------------------------
 ! (C) Crown copyright 2024 Met Office. All rights reserved.
 ! The file LICENCE, distributed with this code, contains details of the terms
 ! under which the code may be used.
 !-----------------------------------------------------------------------------

While the date should be correct for new code, updating the date for each change
to the code is not important. Note that some older code has a different form of
copyright referencing the Queen's Printer. These must be left as they are.

Quick List of most-commonly forgotten things
--------------------------------------------

* File names must match the name of the module they contain.
* ``implicit none`` must be included in every module and every
  subroutine/function it contains.
* ``use`` statements must have ``only`` statements and must only use things that
  are actually used.
* Procedures must have Doxygen comments with a short (typically one-line)
  ``@brief``; an optional ``@details`` and a description and intent of each
  input ``@param`` to subroutines and functions.
* For readability and for use of diff tools, aim to limit lines to 80
  characters.
* Every allocate must have a matching deallocate.
* A ``kind`` must be attached to all real variables, real literals and literal
  arguments in subroutine calls.
* When removing variables from code, you must clean-up related variable
  declarations and comments.
* You must not have any trailing white-space characters at the end of
  lines, including comment lines.

Calling hierararchy
-------------------

The diagram below gives a high level overview of how the various parts of a
model system should relate to each other. Enforcement of the rules best ensures
that parts of the code-base remain independent from other parts.

.. figure:: images/BroadCodeStructure.svg
   :width: 80%

   Calls can be made "down" but must never be made "up" the hierarchy. In some
   circumstances, calling horizontally across the diagram (for example, from one
   component to another) is acceptable, but caution should be used.

General syntax and style rules
------------------------------

Other than exceptions noted below:

* Lower-case must be used for all code. Comments and text output should follow
  normal grammar, so will allow upper-case.
* Meaningful names should be used for variables or program units. Where names
  have more than one word, use the under-score character to separate them. For
  example ``cell_number``.

Exceptions to the above two rules for naming variables and program units are as
follows:

* Parameters used to enumerate one option among several can use upper case. For
  example, the parameters used to define the level of a log message (see below)
  or the parameters used to define the function space or argument types in LFRic
  (such as ``W0`` and ``GH_READ``).
* Variables that would be widely recognised as representing a variable in a
  mathematical formula can use upper-case. For instance ``Cp`` for the heat
  capacity of dry air or ``Rd`` for the gas constant.
* Where LFRic calls routines in libraries that break the LFRic rules, follow the
  rule commonly used for the library (e.g. in the documentation for the
  library). For example, ESMF uses camel case such as
  ``ESMF_VMGetCurrent``. Follow ESMF's use of camel case when calling their
  functions. Do not write ``call esmf_vmgetcurrent``.
* The suffixes ``_mod``, ``_type`` must be used for modules and Fortran types. A
  module whose key role is to define a type and which has a type constructor
  must use type names with the same prefix. For example, the operator_mod module
  may include an operator_type and an operator_constructor.
* Constructors and destructors must be given a suffix that identifies them
  (either ``_constructor``, ``_destructor`` or ``_init``, ``_final`` can be
  appropriate).
* Where more than one constructor exists, the main constructor must have
  just the chosen constructor suffix. Additional constructors must have the same
  name but with a further descriptive suffix e.g. if a URL type has a
  ``url_constructor`` then a constructor that constructs a URL by copying
  features of another URL may be called ``url_constructor_copy``.

Use comments and Fortran labelling appropriately to clarify the structure of
heavily-nested loops, long ``if``- or loop- blocks, or if-blocks with many
``else if`` statements.

Preferred spelling and spacing
------------------------------

It makes sense to standardise spellings for words that have more than one
accepted spelling, as it makes searching the code-base for names easier.

* Use British English spellings - so use "-ise" rather than "-ize", for example
  in "initialise"
* The plural of "halo" has two accepted spellings in British English: "halos"
  and "haloes". We have made the, rather arbitrary, decision to standardise on
  "halos".

UM style-guide requests that spaces and blank lines are used "where appropriate"
to improve code readability.

Within LFRic, the following style has been adopted for spaces. For long and
complex lines, some spaces may be omitted if it is felt that they impact
readability. But omit spaces in a consistent way.

* Spaces after commas, except within declarations with ``:`` such as ``dofmap(:,:)``.
* Spaces either side of symbols for arithmetic expressions (``+``, ``-``, ``*``, ``/``),
  logical expressions (``==``, ``<=`` etc.), assignments (``=``, ``=>``) and other
  separators(``::``).
* Spaces outside parenthesis containing logical statements such as ``if (...)
  then`` statement expressions.
* No space between subroutine names, function names, array variable names and
  the open-parenthesis that may follow.
* There must be no "trailing whitespace": spaces at the end of lines. This rule
  applies to comment lines too.
* Spaces at the beginning and at the end of a comma-separated list within
  parentheses.

  * Where there is a single item within parenthesis, spaces are optional, but be
    consistent with nearby code
  * Spaces after the comment symbols (``!``, ``!>``)
  * Spaces before continuation-line markers (``&``)

A short illustration of acceptable spacing:

.. code-block:: fortran

      integer(i_def), allocatable :: dofmap(:,:)
      type(field_type), pointer :: exner_theta
      call invoke( my_kernel_type( field1, field2 ) )
      varbeta = 1.0 - varalpha
      if (use_wavedynamics) then
        call invoke( aX_plus_bY( u_adv, varbeta, state_after_slow(igh_u), &
                     varalpha, state(igh_u) ) )

Robust Coding
-------------

For Fortran 2003, modules, and types within them, must be private by
default. Where the syntax is supported, private should be specified as an
attribute.

Procedure variables and pointers must not be initialised when they are declared
as this gives them the save attribute, making the code unsafe for use within
shared-memory parallel regions. Note, it is safe to initialise pointers and
variables declared within a derived type.

.. code-block:: fortran

 subroutine my_sub()
   type(my_type), pointer :: my_pointer => null() ! Wrong!

Unassociated pointers can be unpredictable. Ensure that pointers are nullified
or initialised early in a procedure to reduce the risk of using unassociated
pointers later in the routine.

.. code-block:: fortran

 subroutine my_sub()
   type(my_type), pointer :: my_pointer
   nullify(my_pointer) ! Right!

C code must be called only using the ISO Fortran C interoperability features.

Rules about character variables:

* Character variables that are inputs (``intent(in)`` or
  ``intent(inout)``) to a procedure must be declared
  ``character(*)``. If the declaration specifies a length, problems
  can occur when passing in strings of a different length.
* The ``trim`` function should be used when passing a character
  variable to a procedure.


LFRic-specific standards - the basics
-------------------------------------

These rules apply to the LFRic infrastructure code and to science code that
runs within LFRic model configurations.

General rules
^^^^^^^^^^^^^

* LFRic comments use markup so that the interface documentation can be generated
  for rendering in a browser or a PDF document. Therefore, comment the LFRic
  interface and algorithm-layer code with appropriate markup.

  * Program units must all have, at the very least, a one line description that
    is prefixed with the Doxygen directive ``@brief``
  * If appropriate, more detailed descriptions can be added using the
    ``@details`` directive.
  * Each input argument must be described using
    the ``@param`` directive with the intent and name of variable.

    .. code-block:: fortran

       !> @brief A brief description of the program unit.
       !> @details A longer description of the program unit where a brief one is
       !!          insufficient. The longer description can go over several lines.
       !> @param[out]    output_arg   Description of an output argument
       !> @param[in,out] inoutput_arg Description of an updated argument
       !> @param[in]     input_arg    Description of an input argument

* A :ref:`field collection <field collection>` must be declared ``intent(in)`` as
  it is the set of fields in the collection, and not the field collection object
  itself, that is modified.
* Do not use ``write`` or ``print`` statements to write text to standard
  output. Use the LFRic :ref:`logger <logging>`.
* Do not hold any variables or objects in module scope as this can prevent the
  module being used in two different concurrent processes in which different
  processes want to configure the variable or object differently.

Calling kernels
^^^^^^^^^^^^^^^

* In algorithms, group kernel and built-in calls into a single ``call invoke``
  where possible. Review code, and if you can safely move things around to group
  more kernels, then do so. Grouping allows PSyclone to be used to generate
  optimised code that takes into account dependencies between operations. Note,
  PSyclone cannot see dependencies between separate invoke calls.

  For example, if:

  .. code-block:: fortran

      call invoke( setval_c(dtheta, 0.0_r_def) )
      call invoke( setval_c(du, 0.0_r_def) )
      call invoke( X_divideby_Y(u_physics, u, dA) )
      call invoke( extract_w_kernel_type(w_physics, u_physics) )

  becomes:

  .. code-block:: fortran

      call invoke( setval_c(dtheta, 0.0_r_def),                 &
                   setval_c(du, 0.0_r_def),                     &
                   X_divideby_Y(u_physics, u, dA)               &
                   extract_w_kernel_type(w_physics, u_physics) )

* ``invoke`` calls can include a name. The name must be a valid subroutine name. The
  subroutine name generated by PSyclone in the PSy-layer will be the name
  prefixed with ``invoke_`` making it easy to find the link between the invoke
  call and the PSy layer subroutine it relates to.

  .. code-block:: fortran

     call invoke( name = "map_to_physics",                     &
                  setval_c(dtheta, 0.0_r_def),                 &
                  setval_c(du, 0.0_r_def),                     &
                  X_divideby_Y(u_physics, u, dA)               &
                  extract_w_kernel_type(w_physics, u_physics) )

* By convention, kernels and built-ins list output or input output arguments
  first in their argument list followed by input arguments. Aim to follow this
  convention in other code too.

  * Kernels are designed to be called multiple times per invoke, each call
    modifying only a subset of the field. Therefore, kernel output arguments
    need always to have intent(inout) and not intent(out). To understand the
    true intent, examine the kernel metadata.
* PSyclone must be able to access field declaration for fields being passed into
  an invoke. Therefore do not pass fields obtained from modules or from function
  calls. Rather, declare a field pointer of the correct type, point it to your
  module field or function call, then pass the field pointer.

PSyclone and psykal-lite code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :ref:`PSyKAl <psykal and datamodel>` design implemented as part of LFRic has
several standards which define the directory and code structure of science code,
and rules about what can be done within kernels and algorithms. Documentation
should be referred to or advice sought from the LFRic team.

Where scientific requirements of a kernel are not met by the existing PSyclone
LFRic API, so-called "psykal-lite" code can be written to provide the link
between algorithm and kernel that would otherwise be generated by the PSyclone
code generator. Such code must be drawn to the attention of the LFRic team and
PSyclone developers. It must be plausible that PSyclone can be updated to
generate code meeting the new requirements. LFRic tickets or PSyclone issues
must be raised to extend support to the new requirement, and the tickets
must be accepted by an experienced PSyclone developer.

The psykal-lite code must include comments that reference the tickets and issues.

Even if a kernel cannot be called from the generated PSy layer and
requires psykal-lite code, it must still be a well-formed kernel,
complete with appropriate kernel metadata (including place-holder
metadata if there are inputs that are not currently supported). That
way, it will require minimal rewriting when PSyclone is updated.

Declaring precision
^^^^^^^^^^^^^^^^^^^

* LFRic supports mixed-precision codes where an application may combine
  substantial amounts of codes running at different real precisions.

  * All real and integer variables and parameters must be declared
    with a specified ``kind``.
  * All literal real variables must be given a kind using the following syntax,
    where ``r_mykind`` references a specified kind parameter:

   .. code-block:: fortran

      my_var = 1.23_r_mykind

   This rule exists because the following alternative can give a numerically-different result:

   .. code-block:: fortran

      my_var = real( 1.23, r_mykind )

  * Literal integers do not need a kind unless their value is too
    large to hold in a native integer. However, if a literal integer
    is an argument in a call to a procedure in which the dummy
    argument is declared with a specific kind, then it must be given
    a matching kind.
  * Unless there is good reason, avoid mixing different precisions within a
    program unit.

Unit Testing
------------

Some key standards for writing tests are as follows:

* Name the test after the unit being tested as follows: the test file for a unit
  mycalc_kernel_mod.f90 will be called mycalc_kernel_mod_test.pf.

  * Large test files can be split by adding a further suffix: for example
    ``mycalc_kernel_mod_test_cubedsphere.pf`` and
    ``mycalc_kernel_mod_test_planar.pf`` might logically split tests into
    testing data from a cubed-sphere mesh and data from planar mesh. However, do
    consider whether testing could be simplified by breaking the module down
    into smaller components.
* The unit_test directory tree should mirror the source directory tree to make
  finding test files easy.
* Use named arguments when initialising configuration items with the feign
  functions, as the order of variables in the configuration cannot be guaranteed.
* Do not create LFRic fields in tests. Where field data is required simplified
  canned fields have been created.
* Use a class extending TestCase only if you need to pass fixtures to the
  tests. Otherwise use the ``@before`` and ``@after`` directives to garnish your
  set-up and tear-down subroutines.
