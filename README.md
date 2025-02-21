LFRic Core
==========
[![Docs](https://github.com/MetOffice/lfric_core/actions/workflows/deploy-docs.yml/badge.svg?branch=main)](https://github.com/MetOffice/lfric_core/actions/workflows/deploy-docs.yml)

Location for LFRic infrastructure documentation. The source code will be migrated here at a later date.

To build the sphinx docs you will need to use the `lfric-core-docs/1.0.0` environment made specifically for this (it contains Sphinx and the required packages), be aware that this module is incompatible with the main lfric environment modules due to different Python dependencies. To build use `make html` in the documentation directory. `make help` will give you the other options available.
