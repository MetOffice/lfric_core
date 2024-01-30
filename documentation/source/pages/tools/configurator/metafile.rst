.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

Configuration Metadata File Extensions
======================================

`Rose`_ uses a modified Windows ``.ini`` file syntax to store the metadata used
to describe the GUI it presents and thence the namelists it creates. In order
to effectively load those namelists further information is required. Therefore
the Configurator tools look for various additional key/value pairs.

In order to prevent this additional information disrupting Rose they are
commented out using the normal exclamation mark: ``!``.

.. _rose: https://metomi.github.io/rose/doc/html/index.html

Extensions to Namelist Definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Normally a namelist may only appear once in a namelist file however it can be
appropriate for them to appear multiple times. In order to mark a namelist as
suitible for multiple instances set the ``duplicate`` value.

When there are multiple instances of a namelist a means is needed to
distinguish them. The ``instance_key_member`` identifies one of the members of
the namelist as that index key.

An example showing both of these::

    [namelist:instance]
    duplicate = .true.
    instance_key_member = index

    [namelist:instance=index]
    ...

Extension to Namelist Member Definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Although Rose understands the data type a namelist member must have it does
not understand the concept of Fortran "kinds." To add this knowledge a ``kind``
property is added::

    [namelist:thing=value]
    type=real
    !kind=default

The value maps on to Fortran kinds as shown in this table:

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

While a namelist string may have any length a Fortran character array has a
fixed length. In order to let the configurator know what length of string to
read a particular member into use the ``string_length`` property::

    [namelist:thing=string]
    type=character
    !string_length=filename

The length is taken from a set of predefined values:

+----------------+-------------------+
| Property Value | Fortran Parameter |
+================+===================+
| default        | str_def           |
+----------------+-------------------+
| filename       | str_max_filename  |
+----------------+-------------------+

A data field in Rose can have a fixed set of values. We may well want such
fields to be treated as enumerated types. Marking such a field with the
``enumeration`` property allows this to happen::

    [namelist:thing=choice]
    value-titles=First, The Second One, Other
    values='first', 'second', 'third'
    !enumeration=true

Array type structures are supported by Rose but it only understands that values
may be scalar or vector, both fixed and variable length. It can be useful to
limit the size of a vector value by a user specified amount. This is the
purpose of the ``bounds`` property::

    [namelist:thing=vector]
    type=real
    !kind=default
    length=:
    !bounds=namelist:thing=vector_size

Sometimes it is useful for a piece of configuration to not be specified by the
user but derived from others. A classic example of this is converting a user
specified angle in degrees into radians::

    [namelist:thing=radians]
    type=real
    !kind=default
    !expression=namelist:thing=degreessource:constants_mod=PI / 180.0_r_def
