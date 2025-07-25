.. -----------------------------------------------------------------------------
    (c) Crown copyright 2024 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _extended rose metadata:

Configuration Metadata File Extensions
======================================

`Rose`_ uses a modified Windows ``.ini`` file syntax to store the
metadata used to configure its GUI to generate namelists. The
:ref:`Configurator<configurator>` is a tool that generates Fortran
source to read such namelists, but requires extensions to the standard
Rose metadata to do so.

This section describes these additional key-value pairs. To prevent
this additional information disrupting Rose, within the Rose metadata
files they are commented out using the normal exclamation mark: ``!``.

.. _rose: https://metomi.github.io/rose/doc/html/index.html

Extensions to Namelist Definitions
----------------------------------

Commonly, a namelist appears once in a namelist configuration
file. However, some applications need some namelists to appear
multiple times. To mark a namelist as suitable for multiple instances
set the ``duplicate`` value (using standard Rose metadata).

When there are multiple instances of a namelist a means is needed to
allow the Configurator to distinguish each one. The extended-metadata
``instance_key_member`` key identifies one of the members of the
namelist as an index key. The member variable must be one that is set
to a unique value in each instance of the namelist in the
configuration file.

.. admonition:: Example

 An example of metadata for a namelist that may have multiple
 configurations within the same application::

  [namelist:partitioning]
  duplicate=true
  !instance_key_member=mesh_type

  [namelist:partitioning=mesh_type]

 The namelist configuration file might look as follows::

  &partitioning
  mesh_type             = 'destination',
  partitioner           = 'cubedsphere',
  panel_decomposition   = 'auto',
  /

  &partitioning
  mesh_type             = 'source',
  partitioner           = 'cubedsphere',
  panel_decomposition   = 'auto',
  /

 The instance key ``mesh_type`` variables are set to different values:
 ``destination`` in one instance and ``source`` in the other. The
 different setting will be the key to distinguishing the two sets of
 values.

Extension to Namelist Member Definitions
----------------------------------------

Although Rose understands the data type a namelist member must have, it
does not understand the concept of Fortran "kinds". A ``kind``
property provides this information::

    [namelist:thing=value]
    type=real
    !kind=default

The kind value maps on to one of the LFRic kind values defined in the
``constants_mod`` module as shown in this table:

+---------+---------------+--------------+
| Type    | Namelist Kind | Fortran Kind |
+=========+===============+==============+
| logical | default       | l_def        |
|         +---------------+--------------+
|         | native        | l_native     |
+---------+---------------+--------------+
| integer | default       | i_def        |
|         +---------------+--------------+
|         | short         | i_short      |
|         +---------------+--------------+
|         | medium        | i_medium     |
|         +---------------+--------------+
|         | long          | i_long       |
+---------+---------------+--------------+
| real    | default       | r_def        |
|         +---------------+--------------+
|         | native        | r_native     |
|         +---------------+--------------+
|         | single        | r_single     |
|         +---------------+--------------+
|         | double        | r_double     |
|         +---------------+--------------+
|         | second        | r_second     |
+---------+---------------+--------------+

While a namelist string may have any length, a Fortran character array has a
fixed length. The ``string_length`` property allows the Configurator
to define an appropriate length::

    [namelist:thing=string]
    type=character
    !string_length=filename

The length is taken from a set of values defined in the LFRic
``constants_mod`` module:

+----------------+-------------------+
| Property Value | Fortran Parameter |
+================+===================+
| default        | str_def           |
+----------------+-------------------+
| filename       | str_max_filename  |
+----------------+-------------------+

A data field in Rose can have one of a fixed set of values. Marking
such a field with the ``enumeration`` property results in the
Configurator defining an enumerator for use in the code, which can be
safer and more readable than referring to particular values or string
settings.

.. admonition:: Example

 For the following metadata defining an enumerator::

    [namelist:thing=choice]
    value-titles=First, The Second One, Other
    values='first', 'second', 'third'
    !enumeration=true

 the Configurator will create an enumeration object or variable
 ``choice`` which will be set to one of the values ``choice_first``,
 ``choice_second`` or ``choice_third`` depending on whether the
 namelist variable is set to the string ``first``, ``second`` or
 ``third``. It will be set to none of these values if the variable
 is undefined by the namelist configuration.

Array type structures are supported by Rose but it only understands that values
may be scalar or vector, both fixed and variable length. It can be useful to
limit the size of a vector value by a user specified amount. This is the
purpose of the ``bounds`` property::

    [namelist:thing=vector]
    type=real
    !kind=default
    length=:
    !bounds=namelist:thing=vector_size

Sometimes it is useful for a configuration setting to be specified by
the user in one set of units, but used in the code in another. A
classic example of this allowing a user to specify coordinates in
degrees in the namelist (where, commonly, values are in round
numbers), but use radians in the code::

    [namelist:thing=longitude_radians]
    type=real
    !kind=default
    !expression=namelist:thing=longitude_degrees:constants_mod=PI / 180.0_r_def

The code using this namelist will be able to access the value in
radians or the value in degrees.
