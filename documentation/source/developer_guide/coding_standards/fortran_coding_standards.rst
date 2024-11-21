.. ------------------------------------------------------------------------------
     (c) Crown copyright 2023 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------

.. _fortran coding standards:

Fortran Coding Standards for LFRic
==================================

Rules for coding standards are driven by the need for readability and
consistency (which aids readability).

While some people are happy to read inconsistent code, other people find
inconsistency to be distracting and their needs should be respected.

Some of the rules are required to meet the technical needs of the LFRic
core and application code structure and organisation.

LFRic coding standards start from the UM standards which cover Fortran 95. These
rules should be followed unless the LFRic rules on this page override them.

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
all new LFRic code. A Fortran example for 2022 is:

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
model system should relate to each other.

Of particular note is the requirement that calls can be made "down" but must
never be made "up" the hierarchy. In some circumstances, calling horizontally
across the diagram (for example, from one component to another) is acceptable,
but caution should be used.

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

Preferred spellings
-------------------

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
* Spaces att he beginning and at the end of a comma-separated list within
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

Fortran 2003 related aspects:

* The suffixes ``_mod``, ``_type`` must be used for modules and Fortran types. A
  module whose key role is to define a type and which has a type constructor
  must use type names with the same prefix. For example, the operator_mod module
  may include an operator_type and an operator_constructor.
* Constructors and destructors must be given a suffix that identifies them
  (either ``_constructor``, ``_destructor`` or ``_init``, ``_final`` can be
  appropriate).
* Where more than one constructor exists, the main constructor must have the
  just the chosen constructor suffix. Additional constructors must have the same
  name but with a further descriptive suffix e.g. if a URL type has a
  ``url_constructor`` then a constructor that constructs a URL by copying
  features of another URL may be called ``url_constructor_copy``.

Use comments and Fortran labelling appropriately to clarify the structure of
heavily-nested loops, long ``if``- or loop- blocks, or if-blocks with many
``else if`` statements.

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

LFRic-specific standards - the basics
=====================================

These rules apply to the LFRic infrastructure code and to science code that
runs within LFRic model configurations.

* LFRic uses Doxygen to generate interface documentation for the LFRic interface
  and for the algorithm layer of science code. Therefore, comment the LFRic
  interface and algorithm-layer code with Doxygen markup.

  * Program units must all have, at the very least, a one line description that
    is prefixed with the Doxygen directive ``@brief``
  * If appropriate more detailed description uses the ``@details`` directive.
  * Each input argument must (with the following exception) be described using
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
  output. Use the LFRic :ref:`logger <logger>`.
*
