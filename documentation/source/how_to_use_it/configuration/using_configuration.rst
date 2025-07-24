.. -----------------------------------------------------------------------------
    (c) Crown copyright 2024 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _using configuration:

Introduction
============

Applications that use the :ref:`extended Rose metadata <extended rose
metadata>` can run the LFRic :ref:`Configurator <configurator>` as
part of the application build process. The Configurator generates all
the code required to read namelist configuration files. The code also
makes the configuration information available in a user-friendly format.

This section describes how to load an application configuration into
the application, and how code can use the various types of application
configuration.

Loading the configuration
-------------------------
