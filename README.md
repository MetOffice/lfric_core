LFRic Core
==========
[![Docs](https://github.com/MetOffice/lfric_core/actions/workflows/deploy-docs.yml/badge.svg?branch=main)](https://github.com/MetOffice/lfric_core/actions/workflows/deploy-docs.yml)

Location for LFRic infrastructure documentation. The LFRic infrastructure 
source code will be migrated from its current Subversion repository to the 
``main`` branch in November 2025.

On the Met Office Azure Spice machine the main LFRic module environment 
contains all the required packages to build the documentation. To build use 
`make html` in the documentation directory. `make help` will give you the other
options available. Additionally, `make deploy` will build a copy of the 
documentation and deploy it to a directory in `$(HOME)/public_html` named after
the git branch.

Any changes should be developed on a fork of the `main` branch. Do not 
target the `trunk` branch as this is currently used to synchronise subversion 
`trunk`, and we will remove this branch following git migration in November 2025.
