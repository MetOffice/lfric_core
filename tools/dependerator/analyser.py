#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
# Generate a make file snippet holding build information about a source file.
#
# This snippet consists of a dependency list built up from which modules a
# file "use"s.
#
# Additionally a separate snippet file is created for each "program" unit found
# listing all the modules which go to make that program.
#
# These snippets may then be "include"ed into other make files.

from __future__ import print_function;

from abc import ABCMeta, abstractmethod
import os
import os.path
import re
import shutil
import subprocess

'''
Examine Fortran source and build dependency information for use by "make".
'''

###############################################################################
# Interface for analysers
#
class Analyser(object):
    __metaclass__ = ABCMeta

    ###########################################################################
    # Examine a source file and store dependency information in the database.
    #
    # Arguments:
    #   sourceFilename - The name of the object to scan.
    #
    @abstractmethod
    def analyse( self, sourceFilename ):
        pass

###############################################################################
# Examine a namelist description file for dependencies.
#
class NamelistDescriptionAnalyser(Analyser):
    ###########################################################################
    # Constructor
    # Arguments:
    #   logger   - A logger object to write output to.
    #   database - FileDependencies object to hold details.
    #
    def __init__( self, logger, database ):
        self._logger   = logger
        self._database = database

    ###########################################################################
    # Scan a namelist description file and harvest dependency information.
    #
    # Arguments:
    #   sourceFilename - File object to scan.
    #
    def analyse( self, sourceFilename ):
        if not sourceFilename.endswith( '.nld' ):
            raise Exception( 'File doesn''t look like a namelist description' \
                             + ' file: ' + sourceFilename )

        self._logger.logEvent( '  Scanning ' + sourceFilename )
        self._database.removeFile( sourceFilename )

        with open( sourceFilename, 'r' ) as sourceFile:
            for line in sourceFile:
                match = re.match( r'^\s*namelist\s+(\S+)', line, \
                                  flags=re.IGNORECASE )
                if match is not None:
                    fortranFilename = '{}_configuration_mod.f90' \
                                      .format( match.group( 1 ) )
                    self._database.addFileDependency( fortranFilename, \
                                                      sourceFilename )

###############################################################################
# Examine Fortran source code for module dependencies.
#
class FortranAnalyser(Analyser):
    ###########################################################################
    # Constructor
    #
    # Arguments:
    #   logger        - A Logger object to write output to.
    #   ignoreModules - A list of module names to ignore.
    #   database      - FortranDatabase object to hold details.
    #
    def __init__( self, logger, ignoreModules, database ):
        self._logger          = logger
        self._ignoreModules   = ignoreModules
        self._database        = database

        self._fpp = os.getenv( 'FPP', None )
        if not self._fpp:
            raise Exception( 'No Fortran preprocessor provided in $FPP' )
        self._fpp = self._fpp.split()

    ###########################################################################
    # Scan a Fortran source file and harvest dependency information.
    #
    # Arguments:
    #   sourceFilename - The name of the object to scan.
    #
    def analyse( self, sourceFilename ):
        if sourceFilename.endswith( '.F90' ):
            self._logger.logEvent( '  Preprocessing ' + sourceFilename )
            preprocessCommand = self._fpp + [sourceFilename]
            preprocessor = subprocess.Popen( preprocessCommand, \
                                             stdout=subprocess.PIPE, \
                                             stderr=subprocess.PIPE )
            processedSource, errors = preprocessor.communicate()
            if preprocessor.returncode:
                raise subprocess.CalledProcessError( preprocessor.returncode, \
                                              " ".join( preprocessCommand ) )
        elif sourceFilename.endswith( '.f90' ):
            with open( sourceFilename, 'r' ) as sourceFile:
                processedSource = sourceFile.read()
        else:
            raise Exception( 'File doesn''t look like a Fortran file: ' \
                            + sourceFilename )

        self._logger.logEvent( '  Scanning ' + sourceFilename )
        self._database.removeFile( sourceFilename )

        def addDependency( programUnit, prerequisiteUnit, reverseLink=False ):
            if prerequisiteUnit in self._ignoreModules:
                self._logger.logEvent( '      - Ignored 3rd party prerequisite' )
            elif prerequisiteUnit in modules:
                self._logger.logEvent( '      - Ignored self' )
            elif prerequisiteUnit in dependencies:
                self._logger.logEvent( '      - Ignore duplicate prerequisite' )
            else:
                dependencies.append( prerequisiteUnit )
                self._database.addModuleCompileDependency( programUnit, \
                                                           prerequisiteUnit )
                if reverseLink:
                    self._database.addModuleLinkDependency( prerequisiteUnit, \
                                                            programUnit )
                else:
                    self._database.addModuleLinkDependency( programUnit, \
                                                            prerequisiteUnit )

        programUnit = None
        modules = []
        dependencies = []
        for line in processedSource.splitlines():
            match = re.match( r'^\s*PROGRAM\s+(\S+)', line, \
                                flags=re.IGNORECASE)
            if match is not None:
                programUnit = match.group( 1 )
                self._logger.logEvent( '    Contains program: ' + programUnit )
                self._database.addProgram( programUnit, sourceFilename )

            match = re.match( \
                     r'^\s*MODULE\s+(?!(?:PROCEDURE|SUBROUTINE|FUNCTION)\s+)(\S+)', \
                              line, flags=re.IGNORECASE)
            if match is not None:
                programUnit = match.group( 1 )
                self._logger.logEvent( '    Contains module ' + programUnit )
                modules.append( programUnit )
                self._database.addModule( programUnit, sourceFilename )

            match = re.match( r'^\s*SUBMODULE\s+\((?:([^:]+):)?([^)]+)\)\s+(\S+)', line, \
                              flags=re.IGNORECASE )
            if match is not None:
                ancestorUnit = match.group(1)
                parentUnit  = match.group(2)
                programUnit = match.group(3)

                message = '{}Contains submodule {} of {}'.format(            \
                                                                ' ' * 4,     \
                                                                programUnit, \
                                                                parentUnit )
                if ancestorUnit:
                    message = message + '({})'.format( ancestorUnit )
                self._logger.logEvent( message )

                self._database.addModule( programUnit, sourceFilename )
                # I don't think it's necessary to append this to "modules".
                # It's really just a dependency.
                addDependency( programUnit, parentUnit, True )

            match = re.match( r'^\s*USE\s+([^,\s]+)', line, \
                                flags=re.IGNORECASE)
            if match is not None:
                moduleName = match.group( 1 )
                self._logger.logEvent( '    Depends on module ' + moduleName )
                addDependency( programUnit, moduleName )
