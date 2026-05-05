# LFRic Core Fab Build Scripts

Make sure you have Fab version 2.0.1 or later installed (in addition to all
LFRic core requirements of course).

## Setting up Site- and Platform-specific Settings
Site- and platform-specific settings are contained in 
```$LFRIC_CORE/infrastructure/build/site_specific/${SITE}-${PLATFORM}```
The default settings are in ```site_specific/default``` (and at this
stage each other site-specific setup inherits the values set in the default,
and then adds or modifies settings). The Fab build system provides various
callbacks to the ```config.py``` file in the corresponding directory (details
are in the [Fab documentation](https://metoffice.github.io/fab/fab_base/config.html).

If there is no existing site-specific setup, it is recommended to copy an existing
configuration file (e.g. from ```nci_gadi/config.py```). This act as a template
to indicate where you can specify linker information, select a default compiler
suite etc.

The default setup contains compiler flags for Cray, GNU, Intel-classic (ifort),
Intel-LLVM (ifx), and NVIDIA. For modularity's sake (and to keep the file length
shorter), the default configuration will get the settings from the corresponding
```setup_...py``` script. There is no need for a site to replicate this structure,
existing ```config.py``` scripts show how this can be done.


## Building the Skeleton Apps

In order to build the skeleton apps, change into the directory
```$LFRIC_CORE/applications/skeleton```,
and use the following command:

```
./fab_skeleton.py --nprocs 4 --suite gnu
```
Select an appropriate number of processes to run in parallel, and a compiler
suite, e.g. one of ```gnu```, ```cray```, ```intel-classic```, ```intel-llvm```
or ```nvidia```. Once the process is finished, you should have a binary in the
directory ```./fab-workspace/skeleton-full-debug-COMPILER``` (where
```COMPILER``` is the compiler used, e.g. ```mpif90-gfortran```).

Likely, you will have to setup corresponding compiler and linker options, e.g.
include paths, library paths and libraries to link. The directory ```site_specific```
contains site-specific setup files, and one called ```default``` (which is
used in the above example if no site is selected). If your site does not exist,
create a corresponding directory for your site and platform in the format
```site_platform``` by using an existing site as template. See also the
[Fab documentation](https://metoffice.github.io/fab/fab_base/config.html)
for details about specifying compiler and linker options).
You can then use the command line options ```--site``` and ```--platform```
to pick your setup, e.g.:

```
./fab_skeleton.py --nprocs 4 --site nci --platform gadi --suite intel-llvm
```

This would use the file ```site_specific/nci_gadi/config.py```, and all additional
compiler and linker options defined there.

Using ```./fab_skeleton.py -h``` will show a help message with all supported command line
options (and their default value). If a default value is listed using an environment
variables (```(default: $SITE or 'default')```), the corresponding environment variable
is used if no command line option has been specified.

A different compilation profile can be specified using ```--profile``` option. Note
that the available compilation profiles can vary from site to site (see
[Fab documentation](https://metoffice.github.io/fab/fab_base/config.html) for details).

If Fab has issues finding a compiler, you can use the Fab debug command line option
```--available-compilers```, which will list all compilers and linkers Fab has
identified as being available.
