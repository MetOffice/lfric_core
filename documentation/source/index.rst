.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

:html_theme.sidebar_secondary.remove: true

==========
LFRic Core
==========

The `LFRic Core <https://github.com/MetOffice/lfric_core>`_ project
develops a software infrastructure whose prime requirement is to
support the development of the Momentum atmosphere model. The LFRic
core software also underpins a range of other earth system modelling
requirements and related support tools. The LFRic core development is
led by the `Core Capability Development Team
<CoreCapabilityDevelopmentTeam@metoffice.gov.uk>`_ within the Science
IT group at the Met Office.

.. grid:: 3

    .. grid-item-card::
        :text-align: center

        Information on getting going, from software stacks to testing

        +++
        .. button-ref:: getting_started_index
            :ref-type: ref
            :color: primary
            :outline:
            :expand:

                Getting Started

    .. grid-item-card::
        :text-align: center

        Guide on the code for users and developers

        +++
        .. button-ref:: user_guide_index
            :ref-type: ref
            :color: primary
            :outline:
            :expand:

                User Guide

    .. grid-item-card::
        :text-align: center

        Guide on details of the development process

        +++
        .. button-ref:: developer_guide_index
            :ref-type: ref
            :color: primary
            :outline:
            :expand:

                Developer Guide

Development of the LFRic core infrastructure and the new atmosphere
model are being done within the Momentum Partnership. The Momentum
atmosphere application is developed in a separate repository
accessible to Met Office partners. Key initial aims for the Momentum
atmosphere model are as follows:

- The model will be scientifically as good as the UM atmosphere.
- The model will scale better on future exascale platforms.
- The infrastructure will be flexible enough to support future
  evolutions of the science.

LFRic core has a role to deliver for all of these aims: it has been
written to support the GungHo mixed finite element scheme that is key
to delivering the scientific performance of the Momentum atmosphere
model when running on the cubed-sphere grid that will be used for
global simulations; it is written with scalability and performance in
mind, particularly by being developed alongside the PSyclone Domain
Specific Language (DSL) tool; it follows modern software engineering
practices that aims to separate concerns between scientific and
technical aspects of the code.

.. toctree::
    :maxdepth: 1
    :hidden:

    getting_started/index
    user_guide/index
    developer_guide/index
    API/index