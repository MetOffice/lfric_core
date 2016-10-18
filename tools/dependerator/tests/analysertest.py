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

import dependerator.analyser
import dependerator.database
import utilities.logging

##############################################################################
class NamelistDescriptionAnalyserTest(unittest.TestCase):
    ##########################################################################
    def setUp( self ):
        self._scratchDirectory = tempfile.mkdtemp()
        self._logger       = utilities.logging.NoLogger()
        dbFilename = os.path.join( self._scratchDirectory, 'test.db' )
        self._database     = dependerator.database.SQLiteDatabase( dbFilename )
        self._dependencies = dependerator.database.FileDependencies( self._database )

    ##########################################################################
    def tearDown( self ):
        del self._dependencies
        del self._database
        del self._logger
        shutil.rmtree( self._scratchDirectory )

    ##########################################################################
    def testAnalysis( self ):
        testFilename = os.path.join( self._scratchDirectory, 'test.nld' )
        with open(testFilename, 'w') as nldFile:
            print( '''
namelist foo

  bar  : string( filename ) !gumph
  baz  : enumeration[ thing1, thing2 ]
  qux  : real
  fred : real[ 'qux * 2' ]

end namelist foo
                   '''.strip(), file=nldFile )

        uut = dependerator.analyser.NamelistDescriptionAnalyser( \
                                             self._logger, self._dependencies )
        uut.analyse( testFilename )
        dependencies = list(self._dependencies.getDependencies())
        self.assertEqual( [(u'foo_configuration_mod.f90', \
                           [testFilename])], dependencies )

##############################################################################
class FortranAnalyserTest(unittest.TestCase):
    ##########################################################################
    def setUp( self ):
        self._scratchDirectory = tempfile.mkdtemp()
        self._logger       = utilities.logging.NoLogger()
        dbFilename = os.path.join( self._scratchDirectory, 'test.db' )
        self._database     = dependerator.database.SQLiteDatabase( dbFilename )
        self._dependencies = dependerator.database.FortranDependencies( self._database )

    ##########################################################################
    def tearDown( self ):
        del self._dependencies
        del self._database
        del self._logger
        shutil.rmtree( self._scratchDirectory )

    ##########################################################################
    def testAnalyseProgram( self ):
        testFilename = os.path.join( self._scratchDirectory, 'test.f90' )
        with open(testFilename, 'w') as fortranFile:
            print( '''
program foo

  use constants_mod, only : i_def
  use trumpton_mod, only : hew, pew, barney, mcgrey, cuthbirt, dibble, grub

  implicit none

end program foo
                   '''.strip(), file=fortranFile )

        uut = dependerator.analyser.FortranAnalyser( self._logger, \
                                                     [],           \
                                                     self._dependencies )
        uut.analyse( testFilename )

        programs = list( self._dependencies.getPrograms() )
        self.assertEqual( [u'foo'], programs )

        dependencies = list(self._dependencies.getCompileDependencies())
        self.assertEqual( [], dependencies )

        dependencies = list(self._dependencies.getLinkDependencies( 'foo' ))
        self.assertEqual( [], dependencies )

    ##########################################################################
    def testAnalyseModule( self ):
        uut = dependerator.analyser.FortranAnalyser( self._logger, \
                                                     [],           \
                                                     self._dependencies )

        testFilename = os.path.join( self._scratchDirectory, 'test.f90' )
        with open(testFilename, 'w') as fortranFile:
            print( '''
module foo

  use constants_mod, only : i_def
  use trumpton_mod, only : hew, pew, barney, mcgrey, cuthbirt, dibble, grub

  implicit none

  private

contains

end module foo

module trumpton_mod

end module trumpton_mod
                   '''.strip(), file=fortranFile )
        uut.analyse( testFilename )

        otherFilename = os.path.join( self._scratchDirectory, 'other.f90' )
        with open(otherFilename, 'w') as otherFile:
            print( '''
module constants_mod

end module constants_mod
                   '''.strip(), file=otherFile )
        uut.analyse( otherFilename )


        programs = list( self._dependencies.getPrograms() )
        self.assertEqual( [], programs )

        dependencies = list(self._dependencies.getCompileDependencies())
        self.assertEqual( [(u'foo', testFilename,             \
                            u'constants_mod', otherFilename), \
                           (u'foo', testFilename,             \
                            u'trumpton_mod', testFilename)], dependencies )

        dependencies = list(self._dependencies.getLinkDependencies( 'foo' ))
        self.assertEqual( [(u'foo', testFilename,             \
                            u'constants_mod', otherFilename), \
                           (u'foo', testFilename,             \
                            u'trumpton_mod', testFilename)], dependencies )

    ##########################################################################
    def testAnalyseSubModule( self ):
        uut = dependerator.analyser.FortranAnalyser( self._logger, \
                                                     [],           \
                                                     self._dependencies )

        parentFilename = unicode(os.path.join( self._scratchDirectory, \
                                               'parent.f90' ), "utf-8")
        with open(parentFilename, 'w') as parentFile:
            print( '''
module parent

  implicit none

  private

  type, public :: test_type
  contains
    procedure foo
    procedure bar
    procedure baz
  end type test_type

  interface
    module subroutine foo( this, cheese )
      class(test_type), intent(inout) :: this
      real,             intent(in)    :: cheese
    end subroutine foo

    module subroutine bar( this, teapot )
      class(test_type), intent(inout) :: this
      character(*),     intent(in)    :: teapot
    end subroutine bar

    type(something) module function baz( this )
      class(test_type), intent(in) :: this
    end function baz
  end interface

end module parent

submodule (parent) child3

  implicit none

contains

  type(something) module function baz( this )
    class(test_type), intent(in) :: this
  end function baz

end submodule child3
                   '''.strip(), file=parentFile )
        uut.analyse( parentFilename )

        child1Filename = unicode(os.path.join( self._scratchDirectory, \
                                               'child1.f90' ), 'utf-8')
        with open(child1Filename, 'w') as child1File:
            print( '''
submodule (parent) child1

  implicit none

  type :: secondary_type
  contains
    procedure baz
  end type secondary_type

  interface
    module subroutine baz( this, wibble )
      class(secondary_type), intent(inout) :: this
      real,                  intent(in)    :: wibble
    end subroutine baz
  end interface

contains

  module subroutine foo( this, cheese )

    implicit none

    class(test_type), intent(inout) :: this
    real,             intent(in)    :: cheese

    type(secondary_type) :: thang

    thang = secondary_type()

    write(6, *) cheese
    call thang%baz( 12.7 )

  end subroutine foo

end submodule child1
                   '''.strip(), file=child1File )
        uut.analyse( child1Filename )

        child2Filename = unicode(os.path.join( self._scratchDirectory, \
                                 'child2.f90' ), 'utf-8')
        with open(child2Filename, 'w') as child2File:
            print( '''
submodule (parent) child2

  implicit none

contains

  module procedure bar

    implicit none

    write( 6, *) teapot

  end procedure bar

end submodule child2
                   '''.strip(), file=child2File )
        uut.analyse( child2Filename )

        child3Filename = unicode(os.path.join( self._scratchDirectory, \
                                 'child3.f90' ), 'utf-8')
        with open(child3Filename, 'w') as child3File:
            print( '''
submodule (parent:child1) grandchild

  implicit none

contains

  module subroutine baz( this, wibble )

    implicit none

    class(secondary_type), intent(inout) :: this
    real,                  intent(in)    :: wibble

    write(6, *) wibble

  end subroutine baz

end submodule grandchild
                   '''.strip(), file=child3File )
        uut.analyse( child3Filename )

        programs = list( self._dependencies.getPrograms() )
        self.assertEqual( [], programs )

        dependencies = list(self._dependencies.getCompileDependencies())
        self.assertEqual( [(u'child1', child1Filename,     \
                            u'parent', parentFilename),    \
                           (u'child2', child2Filename,     \
                            u'parent', parentFilename),    \
                           (u'grandchild', child3Filename, \
                            u'child1', child1Filename)], dependencies )

        dependencies = list(self._dependencies.getLinkDependencies( 'parent' ))
        self.assertEqual( [(u'parent', parentFilename,      \
                            u'child1', child1Filename),     \
                           (u'parent', parentFilename,      \
                            u'child2', child2Filename)], dependencies )
        dependencies = list(self._dependencies.getLinkDependencies( 'child1' ))
        self.assertEqual( [(u'child1', child1Filename,      \
                            u'grandchild', child3Filename)], dependencies)

    ##########################################################################
    def testFunctionInModuleName( self ):
        uut = dependerator.analyser.FortranAnalyser( self._logger, \
                                                     [],           \
                                                     self._dependencies )

        testFilename = os.path.join( self._scratchDirectory, 'test.f90' )
        with open(testFilename, 'w') as fortranFile:
            print( '''
module function_thing_mod

  use constants_mod, only : i_def

  implicit none

  private

contains

end module function_thing_mod
                   '''.strip(), file=fortranFile )
        uut.analyse( testFilename )

        otherFilename = os.path.join( self._scratchDirectory, 'other.f90' )
        with open(otherFilename, 'w') as otherFile:
            print( '''
module constants_mod

end module constants_mod
                   '''.strip(), file=otherFile )
        uut.analyse( otherFilename )

        programs = list( self._dependencies.getPrograms() )
        self.assertEqual( [], programs )

        dependencies = list(self._dependencies.getCompileDependencies())
        self.assertEqual( [(u'function_thing_mod', testFilename, \
                            u'constants_mod', otherFilename)], dependencies )

        dependencies = list(self._dependencies.getLinkDependencies( \
                                                       'function_thing_mod' ))
        self.assertEqual( [(u'function_thing_mod', testFilename, \
                            u'constants_mod', otherFilename)], dependencies )

if __name__ == '__main__':
    unittest.main()
