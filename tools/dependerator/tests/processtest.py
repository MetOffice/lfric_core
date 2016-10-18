#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

from __future__ import print_function

import os.path
import shutil
import tempfile
import unittest

import dependerator.database
import dependerator.process
import utilities.logging

##############################################################################
class NamelistDescriptionAnalyserTest(unittest.TestCase):
    ##########################################################################
    def setUp( self ):
        self._logger   = utilities.logging.NoLogger()
        self._scratchDirectory = tempfile.mkdtemp()
        dbFilename = os.path.join( self._scratchDirectory, 'fortran.db' )
        database = dependerator.database.SQLiteDatabase( dbFilename )
        self._fortranDatabase = dependerator.database.FortranDependencies( database )
        self._fileDatabase = dependerator.database.FileDependencies( database )

    ##########################################################################
    def tearDown( self ):
        del self._fileDatabase
        del self._fortranDatabase
        del self._logger
        shutil.rmtree( self._scratchDirectory )

    ##########################################################################
    def testCompileDependencies( self ):
        self.populateDatabase()

        uut = dependerator.process.FortranProcessor( self._logger,          \
                                                     self._fortranDatabase, \
                                                     "objects", "modules" )
        uut.determineCompileDependencies( self._fileDatabase )

        self.assertEqual( [(u'objects/bits/bar.o', [u'modules/bits/baz.mod']), \
                           (u'objects/bobs/qux.o', [u'modules/bits/baz.mod']), \
                           (u'objects/foo.o', [u'modules/bits/bar.mod'])],     \
                          list(self._fileDatabase.getDependencies()))

    ##########################################################################
    def testLinkDependencies( self ):
        self.populateDatabase()

        import sys
        uut = dependerator.process.FortranProcessor( self._logger,          \
                                                     self._fortranDatabase, \
                                                     "objects", "modules" )
        result = list(uut.determineLinkDependencies())

        self.assertEqual( [(u'objects/foo', ['objects/bobs/qux.o', \
                                             'objects/bits/baz.o', \
                                             'objects/bits/bar.o', \
                                             'objects/foo.o'])],   \
                          result )

    ##########################################################################
    def populateDatabase( self ):
        self._fortranDatabase.addProgram( 'foo', 'foo.f90' )
        self._fortranDatabase.addModule( 'bar', 'bits/bar.f90' )
        self._fortranDatabase.addModule( 'baz', 'bits/baz.f90' )
        self._fortranDatabase.addModule( 'qux', 'bobs/qux.f90' )

        self._fortranDatabase.addModuleCompileDependency( 'foo', 'bar' )
        self._fortranDatabase.addModuleCompileDependency( 'bar', 'baz' )
        self._fortranDatabase.addModuleCompileDependency( 'qux', 'baz' )

        self._fortranDatabase.addModuleLinkDependency( 'foo', 'bar' )
        self._fortranDatabase.addModuleLinkDependency( 'bar', 'baz' )
        self._fortranDatabase.addModuleLinkDependency( 'baz', 'qux' )
