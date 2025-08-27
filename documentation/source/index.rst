.. -----------------------------------------------------------------------------
    (c) Crown copyright 2024 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

:html_theme.sidebar_secondary.remove: true

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

.. grid:: 2 2 4 4
    :gutter: 2

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

        Guide on how to use LFRic Core in your applications

        +++
        .. button-ref:: how_to_use_it_index
            :ref-type: ref
            :color: primary
            :outline:
            :expand:

                How to use it

    .. grid-item-card::
        :text-align: center

        Guide on how LFRic Core works under the hood

        +++
        .. button-ref:: how_it_works_index
            :ref-type: ref
            :color: primary
            :outline:
            :expand:

                How it works

    .. grid-item-card::
        :text-align: center

        Guide on what to do when contributing to LFRic Core

        +++
        .. button-ref:: how_to_contribute_index
            :ref-type: ref
            :color: primary
            :outline:
            :expand:

                How to contribute

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

.. grid:: 2 2 4 4
    :gutter: 5
    :padding: 0 0 5 5

    .. grid-item-card::
        :text-align: center

        :material-round:`menu_book;2em`

        +++
        .. button-ref:: glossary_of_terms
            :ref-type: ref
            :color: primary
            :outline:
            :expand:

               Glossary

    .. grid-item-card::
        :text-align: center

        :material-round:`help_center;2em`

        +++
        .. button-ref:: faqs
            :ref-type: ref
            :color: primary
            :outline:
            :expand:

                FAQs

    .. grid-item-card::
        :text-align: center

        :fab:`github;fa-xl`

        +++
        .. button-link:: https://github.com/MetOffice/lfric_core
            :ref-type: ref
            :color: primary
            :outline:
            :expand:

                GitHub

    .. grid-item-card::
        :text-align: center

        :far:`comments;fa-xl`

        +++
        .. button-link:: https://github.com/MetOffice/simulation-systems/discussions
            :ref-type: ref
            :color: primary
            :outline:
            :expand:

                Discussions


.. toctree::
    :maxdepth: 1
    :hidden:

    getting_started/index
    how_to_use_it/index
    how_it_works/index
    how_to_contribute/index
