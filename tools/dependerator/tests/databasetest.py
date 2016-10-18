#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

import unittest
import os.path
import shutil
import sqlite3
import tempfile

import dependerator.database

##############################################################################
class DatabaseTest( unittest.TestCase):
    ##########################################################################
    def setUp( self ):
        self._scratchDirectory = tempfile.mkdtemp()
        self._dbFilename = os.path.join( self._scratchDirectory, 'test.db' )

    ##########################################################################
    def tearDown( self ):
        shutil.rmtree( self._scratchDirectory )

    ##########################################################################
    def testSQLiteDatabase( self ):
        uut = dependerator.database.SQLiteDatabase( self._dbFilename )

        uut.ensureTable( 'test', ["'alpha' integer primary key"] )

        result = uut.query( 'select * from test' )
        self.assertEqual( 0, len(result) )

        result = uut.query( "insert into 'test' values ( 1 )" )
        self.assertEqual( 0, len(result) )

        result = uut.query( 'select * from test' )
        self.assertEqual( 1, len(result) )
        self.assertEqual( 1, result[0][0] )

##############################################################################
class FileDependencyTest(unittest.TestCase):
    ##########################################################################
    def setUp( self ):
        self._scratchDirectory = tempfile.mkdtemp()
        self._dbFilename = os.path.join( self._scratchDirectory, 'file.db' )
        self._database = dependerator.database.SQLiteDatabase( self._dbFilename )

    ##########################################################################
    def tearDown( self ):
        del self._database
        shutil.rmtree( self._scratchDirectory )

    ##########################################################################
    def testAll( self ):
        uut = dependerator.database.FileDependencies( self._database )

        result = uut.getDependencies()
        self.assertEquals( [], list(result) )

        uut.addFileDependency( 'foo.f90', 'bar' )
        result = uut.getDependencies()
        self.assertEquals( [(u'foo.f90', [u'bar'])], list(result) )

        uut.addFileDependency( 'foo.f90', 'baz' )
        uut.addFileDependency( 'qux.f90', 'bar' )
        result = uut.getDependencies()
        self.assertEquals( [(u'foo.f90', [u'bar', u'baz']), \
                            (u'qux.f90', [u'bar'])], list(result) )

        uut.removeFile( 'foo.f90' )
        result = uut.getDependencies()
        self.assertEquals( [(u'qux.f90', [u'bar'])], list(result) )

##############################################################################
class FortranDependencyTest(unittest.TestCase):
    ##########################################################################
    def setUp( self ):
        self._scratchDirectory = tempfile.mkdtemp()
        self._dbFilename = os.path.join( self._scratchDirectory, 'fortran.db' )
        self._database = dependerator.database.SQLiteDatabase( self._dbFilename )

    ##########################################################################
    def tearDown( self ):
        del self._database
        shutil.rmtree( self._scratchDirectory )

    ###########################################################################
    def testAddProgram( self ):
        uut = dependerator.database.FortranDependencies( self._database )
        uut.addProgram( 'foo', 'bar.f90' )

        result = uut.getPrograms()
        self.assertEqual( [u'foo'], result )

        # Should test that filename is correctly stored.

    ###########################################################################
    def testRemoveSourceFile( self ):
        uut = dependerator.database.FortranDependencies( self._database )
        self._populateDB( uut )

        # Delete something from the middle of a dependency chain
        uut.removeFile( 'bar.f90' )

        result = uut.getPrograms()
        self.assertEqual( [u'foo', u'fred'], list(result) )

        result = uut.getCompileDependencies()
        self.assertEqual( [(u'foo', u'foo.f90', u'baz', u'baz.f90'), \
                           (u'fred', u'fred.f90', u'wilma', u'wilma.f90')], \
                          list(result) )

        self.assertEqual( [(u'baz', u'baz.f90', u'foo', u'foo.f90')], \
                          list(uut.getLinkDependencies( 'baz' )) )
        self.assertEqual( [(u'wilma', u'wilma.f90', u'fred', u'fred.f90')], \
                          list(uut.getLinkDependencies( 'wilma' )) )

    ###########################################################################
    def testPrograms( self ):
        uut = dependerator.database.FortranDependencies( self._database )
        self._populateDB( uut )

        programs = uut.getPrograms()
        self.assertEqual( [u'foo', u'fred'], list(programs) )

    ###########################################################################
    def testGetAllFileDependencies( self ):
        uut = dependerator.database.FortranDependencies( self._database )
        self._populateDB( uut )

        dependencies = list( uut.getCompileDependencies() )
        self.assertEqual( [(u'bar', u'bar.f90', u'qux', u'qux.f90'), \
                           (u'foo', u'foo.f90', u'bar', u'bar.f90'), \
                           (u'foo', u'foo.f90', u'baz', u'baz.f90'), \
                           (u'fred', u'fred.f90', u'wilma', u'wilma.f90')], \
                          dependencies )

        self.assertEqual( [(u'bar', u'bar.f90', u'foo', u'foo.f90')], \
                          list(uut.getLinkDependencies('bar')) )
        self.assertEqual( [(u'baz', u'baz.f90', u'foo', u'foo.f90')], \
                          list(uut.getLinkDependencies('baz')) )
        self.assertEqual( [(u'qux', u'qux.f90', u'bar', u'bar.f90')], \
                          list(uut.getLinkDependencies('qux')) )
        self.assertEqual( [(u'wilma', u'wilma.f90', u'fred', u'fred.f90')], \
                          list(uut.getLinkDependencies('wilma')) )

    ###########################################################################
    def _populateDB( self, database ):
        database.addProgram( 'foo', 'foo.f90' )
        database.addModule( 'bar', 'bar.f90' )
        database.addModule( 'baz', 'baz.f90' )
        database.addModule( 'qux', 'qux.f90' )
        database.addProgram( 'fred', 'fred.f90' )
        database.addModule( 'wilma', 'wilma.f90' )

        database.addModuleCompileDependency( 'foo', 'bar' )
        database.addModuleCompileDependency( 'foo', 'baz' )
        database.addModuleCompileDependency( 'bar', 'qux' )
        database.addModuleCompileDependency( 'fred', 'wilma' )

        database.addModuleLinkDependency( 'bar', 'foo' )
        database.addModuleLinkDependency( 'baz', 'foo' )
        database.addModuleLinkDependency( 'qux', 'bar' )
        database.addModuleLinkDependency( 'wilma', 'fred' )

if __name__ == '__main__':
    unittest.main()
