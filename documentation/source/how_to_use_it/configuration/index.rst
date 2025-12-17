.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _core configuration:

Application Configuration
=========================

Application configuration options can be described by Rose metadata
with some :ref:`extensions <extended rose metadata>` to Rose
metadata. With the extensions, the metadata enables a tool called the
:ref:`Configurator <configurator>` to generate the code required to
read configuration namelists into the application.

Once a namelist is read into the application, the generated code
stores it and makes it available. This section summarises how the code
within the application can access and use these configuration options.

.. toctree::
    :maxdepth: 1

    using_configuration
