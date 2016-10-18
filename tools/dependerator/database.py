#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
# Manages a database of dependency information.
#

from __future__ import print_function

from abc import ABCMeta, abstractmethod
import sqlite3

##############################################################################
# Databases throw this exception.
#
class DatabaseException( Exception ):
    pass

##############################################################################
# Basic backend database functionality.
#
class _Database(object):
    __metaclass__ = ABCMeta

    ##########################################################################
    # Makes sure a described table exists in the database.
    #
    # Arguments:
    #   name    - String by which the table will be known.
    #   columns - List of name/type/modifiers tuples.
    #
    @abstractmethod
    def ensureTable( self, name, columns ):
        pass

    ##########################################################################
    # Passes an SQL query to the backend.
    #
    # \deprecated Clearly this breaks all encapsulation. It should be replaced
    #             with a series of objects from which to build an abstract
    #             query.
    #
    # Arguments:
    #   query - Potentially multi-line SQL query.
    #
    @abstractmethod
    def query( self, query ):
        pass

##############################################################################
# Database backend.
#
class SQLiteDatabase(_Database):
    ##########################################################################
    # Default constructor.
    #
    # Arguments:
    #   filename - The filename of the database.
    #
    def __init__( self, filename ):
        super(SQLiteDatabase, self).__init__()

        self._database = sqlite3.connect( filename, timeout=5.0 )
        self._database.row_factory     = sqlite3.Row

    ###########################################################################
    # Destructor.
    #
    def __del__( self ):
        self._database.commit()
        self._database.close()

    ###########################################################################
    def ensureTable( self, name, columns ):
        columnDefinitions = []
        for columnDetails in columns:
            columnDefinitions.append( ' '.join( columnDetails ) )
        query = 'CREATE TABLE IF NOT EXISTS {} ( {} )'
        with self._database:
            self._database.execute(query.format( name, \
                                              ', '.join( columnDefinitions ) ) )

    ###########################################################################
    def query( self, query ):
        cursor = self._database.cursor()
        if not isinstance( query, list):
            query = query.strip().replace( '\n', ' ' )
            cursor.execute( query )
        else:
            query = ';\n'.join( query ) + ';'
            cursor.executescript( query )
        return cursor.fetchall()

###############################################################################
# Basic file dependencies.
#
class FileDependencies(object):
    ###########################################################################
    # Default constructor.
    #
    # Arguments:
    #   database - Database object to use as backend.
    #
    def __init__( self, database ):
        self._database = database
        self._database.ensureTable( 'file_dependency',                \
                                    (('file', 'TEXT', 'NOT NULL'), \
                                     ('prerequisite', 'TEXT', 'NOT NULL')) )

    ###########################################################################
    def removeFile( self, filename ):
        query = 'DELETE FROM file_dependency WHERE file="{}"'.format( filename )
        self._database.query( query )

    ###########################################################################
    def addFileDependency( self, filename, prerequisite ):
        query = "INSERT INTO file_dependency VALUES ('{}', '{}')"
        self._database.query( query.format( filename, prerequisite ) )

    ###########################################################################
    # Get all the file dependency relationships.
    #
    # Arguments:
    # Return:
    #   A generator yielding (filename, [filename]) tuples.
    #
    def getDependencies( self ):
        query = 'SELECT * FROM file_dependency ORDER BY file'
        result = self._database.query( query )

        lastFile      = None
        prerequisites = []
        for row in result:
            if row['file'] != lastFile:
                if prerequisites:
                    yield lastFile, prerequisites
                lastFile = row['file']
                prerequisites = []

            prerequisites.append( row['prerequisite'] )

        if prerequisites:
            yield lastFile, prerequisites

###############################################################################
# Fortran dependencies.
#
class FortranDependencies(object):
    ###########################################################################
    # Default constructor.
    #
    # Arguments:
    #   database - Database object to use as backend.
    #
    def __init__( self, database ):
        self._database = database

        self._database.ensureTable( 'fortran_unit_type', \
                                    [('type', 'TEXT', 'PRIMARY KEY')] )
        self._database.query( 'INSERT OR IGNORE INTO fortran_unit_type VALUES("program")' )
        self._database.query( 'INSERT OR IGNORE INTO fortran_unit_type VALUES("module")' )

        self._database.ensureTable( 'fortran_dependency_type', \
                                    [('type', 'TEXT', 'PRIMARY KEY')] )
        self._database.query( 'INSERT OR IGNORE INTO fortran_dependency_type VALUES("compile")' )
        self._database.query( 'INSERT OR IGNORE INTO fortran_dependency_type VALUES("link")' )

        self._database.ensureTable( 'fortran_program_unit',           \
                                    (('unit', 'TEXT', 'PRIMARY KEY'), \
                                     ('file', 'TEXT', 'NOT NULL'),    \
                                     ('type', 'REFERENCES fortran_unit_type(type)')) )
        self._database.ensureTable( 'fortran_unit_dependency',             \
                                    (('unit', 'TEXT', 'NOT NULL'),         \
                                     ('prerequisite', 'TEXT', 'NOT NULL'), \
                                     ('type', 'REFERENCES fortran_dependency_type(type)')) )

    ###########################################################################
    def removeFile( self, filename ):
        query = ['CREATE TEMPORARY TABLE _units AS SELECT unit '              \
                 + 'FROM fortran_program_unit '                               \
                 + 'WHERE file="{}"'.format( filename ),                      \
                 'DELETE FROM fortran_unit_dependency WHERE unit IN _units', \
                 'DROP TABLE _units',                                         \
                 'DELETE FROM fortran_program_unit WHERE file="{}"'.format( filename )]
        self._database.query( query )

    ###########################################################################
    # Add a program to the database.
    #
    # Arguments:
    #   name     - Name of the program's program unit.
    #   filename - Source file in which program was found.
    #
    def addProgram( self, name, filename ):
        # Changes are transacted to ensure other processes can't find the
        # database with half a program.
        query = "INSERT INTO fortran_program_unit VALUES ( '{}', '{}', 'program' )"
        self._database.query( query.format( name, filename ) )

    ###########################################################################
    # Add a module to the database.
    #
    # Arguments:
    #   name     - The name of the module's program unit.
    #   filename - The source file in which the modules was found.
    #
    def addModule( self, name, filename ):
        query = "INSERT INTO fortran_program_unit VALUES ( '{}', '{}', 'module' )"
        self._database.query( query.format( name, filename ) )

    ###########################################################################
    # Add a compile dependency to the database.
    #
    # Arguments:
    #   module       - module name
    #   prerequisite - module name
    #
    def addModuleCompileDependency( self, module, prerequisite ):
        query = "INSERT INTO fortran_unit_dependency VALUES ( '{}', '{}', 'compile' )"
        self._database.query( query.format( module, prerequisite ) )

    ###########################################################################
    # Add a link dependency to the database.
    #
    # Arguments:
    #   module       - module name
    #   prerequisite - module name
    #
    def addModuleLinkDependency( self, module, prerequisite ):
        query = "INSERT INTO fortran_unit_dependency VALUES ( '{}', '{}', 'link' )"
        self._database.query( query.format( module, prerequisite ) )

    ###########################################################################
    # Gets all the programs from the database
    #
    # Arguments:
    # Return:
    #   List of program names.
    #
    def getPrograms( self ):
        query = "SELECT unit FROM fortran_program_unit " \
                + " WHERE type='program'"
        rows = self._database.query( query )
        return [row['unit']for row in rows]

    ###########################################################################
    # Gets program unit dependencies for a program unit.
    #
    # Arguments:
    #   programUnit - Discover prerequisites for this.
    #
    # Return:
    #   Yields (unit, unit filename, prerequisite, prerequisite filename) tuples
    #
    def getLinkDependencies( self, programUnit ):
        query = '''
SELECT d.unit, df.file AS unit_file, d.prerequisite, pf.file AS prerequisite_file
FROM fortran_unit_dependency AS d, fortran_program_unit AS df, fortran_program_unit as pf
WHERE df.unit=d.unit AND pf.unit=d.prerequisite AND d.unit="{}" AND d.type='link'
ORDER BY d.unit
                '''.strip()
        rows = self._database.query( query.format( programUnit ) )

        for unit, unit_file, prerequisite, prerequisite_file in rows:
            yield unit, unit_file, prerequisite, prerequisite_file

    ###########################################################################
    # Gets all program unit dependencies from database for compilation.
    #
    # Arguments:
    #
    # Return:
    #   Yields (unit, unit filename, prerequisite, prerequisite filename) tuples
    #
    def getCompileDependencies( self ):
        query = '''
SELECT d.unit, df.file AS unit_file, d.prerequisite, pf.file AS prerequisite_file
FROM fortran_unit_dependency AS d, fortran_program_unit AS df, fortran_program_unit as pf
WHERE df.unit=d.unit AND pf.unit=d.prerequisite AND d.type='compile'
ORDER BY d.unit
                '''.strip()
        rows = self._database.query( query )

        for unit, unit_file, prerequisite, prerequisite_file in rows:
            yield unit, unit_file, prerequisite, prerequisite_file
