# LFRic Core

[![ci](https://github.com/MetOffice/lfric_core/actions/workflows/ci.yml/badge.svg)](https://github.com/MetOffice/lfric_core/actions/workflows/ci.yml)

Location for LFRic infrastructure source code and documentation

On the Met Office Azure Spice machine the main LFRic module environment
contains all the required packages to build the documentation. To build use
`make html` in the documentation directory. `make help` will give you the other
options available. Additionally, `make deploy` will build a copy of the
documentation and deploy it to a directory in `$(HOME)/public_html` named after
the git branch.
