#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
# Process previously analysed dependency database. For fun and profit!

from __future__ import print_function;

import os.path

from utilities.path import replaceExtension

###############################################################################
# Process dependency database.
#
class FortranProcessor():
    ###########################################################################
    # Constructor.
    #
    # Arguments:
    #   logger    - A Logger object to write output to.
    #   database  - FortranDatabase object holding details.
    #   objectdir - The directory which holds .o files.
    #   moduledir - The directory which holds .mod files.
    #
    def __init__( self, logger, database, objectDirectory, moduleDirectory):
        self._logger          = logger
        self._database        = database
        self._objectDirectory = objectDirectory
        self._moduleDirectory = moduleDirectory

    ###########################################################################
    # Examine the program unit dependecies and work out the file dependencies.
    #
    # Arguments:
    #   fileStore - FileDependencies object to accept computed dependencies.
    #
    def determineCompileDependencies( self, fileStore ):
        for unit, unitFilename, prerequisite, prerequisiteFilename \
                in self._database.getCompileDependencies():
            self._logger.logEvent( '{} depends on {}'.format( unit, \
                                                              prerequisite ) )

            objectFilename = replaceExtension( unitFilename, 'o' )
            objectPathname = os.path.join( self._objectDirectory, \
                                           objectFilename )

            moduleFilename = replaceExtension( prerequisiteFilename, 'mod' )
            modulePathname = os.path.join( self._moduleDirectory, \
                                           moduleFilename )

            fileStore.addFileDependency( objectPathname, modulePathname )

    ###########################################################################
    # Determine all program units needed to build each program.
    #
    # TODO: Once we have a more recent version of SQLite we could consider
    # doing this at the database level.
    #
    def determineLinkDependencies( self ):
        for program in self._database.getPrograms():
            self._logger.logEvent( 'Program {}'.format( program ) )

            prerequisites = set()
            self._descend( program, prerequisites )
            unit, unit_file, prereq, prereq_file = self._database.getLinkDependencies( program ).next()
            program_object_file = os.path.join( self._objectDirectory, \
                                          replaceExtension( unit_file, 'o' ) )
            prerequisites.add( program_object_file )
            yield os.path.join( self._objectDirectory, program ), \
                  list(prerequisites)

    ##########################################################################
    def _descend( self, programUnit, prerequisites ):
        for unit, unit_file, prereq, prereq_file \
                         in self._database.getLinkDependencies( programUnit ):
            self._logger.logEvent( '  Requires {}'.format( unit ) )

            prereq_object_file = os.path.join( self._objectDirectory, \
                                        replaceExtension( prereq_file, 'o' ) )

            if prereq_object_file in prerequisites:
                self._logger.logEvent( '    Seen already, stopping descent' )
            else:
                prerequisites.add( prereq_object_file )
                self._descend( prereq, prerequisites )
