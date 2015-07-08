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

import sqlite3

##############################################################################
# Databases throw this exception.
#
class DatabaseException( Exception ):
    pass

##############################################################################
# Store and retrieve dependency information.
#
class Dependencies():
    ##########################################################################
    # Constructor.
    #
    # Arguments:
    #   filename - The filename of the database.
    #
    def __init__( self, filename ):
        self._database = sqlite3.connect( filename, timeout=5.0 )

        self._database.isolation_level = 'EXCLUSIVE'
        self._database.row_factory     = sqlite3.Row

        with self._database:
            self._database.executescript( '''
CREATE TABLE IF NOT EXISTS provides (file TEXT NOT NULL, program_unit TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS dependencies (dependor TEXT NOT NULL, dependee TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS programs (name TEXT PRIMARY KEY);
                                          ''' )

    ###########################################################################
    # Destructor.
    #
    def __del__( self ):
        self._database.close()

    ###########################################################################
    # Remove a source file from the database.
    #
    # Arguments:
    #   filename: The filename of the source file to be removed.
    #
    def removeSourceFile( self, filename ):
        # These changes are made a transaction so that the file is removed, or
        # it isn't. No other query can find the database inconsistent.
        with self._database:
            self._database.executescript( '''
CREATE TEMPORARY TABLE _units AS SELECT  program_unit FROM provides WHERE file = "{filename}";
DELETE FROM dependencies WHERE dependor IN _units OR dependee IN _units;
DELETE FROM programs WHERE name IN _units;
DELETE FROM provides WHERE file = "{filename}";
DROP TABLE _units;
                                          '''.format( filename=filename ) )

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
        with self._database:
            self._database.executescript( '''
INSERT INTO provides VALUES ( "{filename}", "{name}" );
INSERT OR REPLACE INTO programs VALUES ( "{name}" );
                                          '''.format( filename=filename, \
                                                      name=name ) )

    ###########################################################################
    # Add a module to the database.
    #
    # Arguments:
    #   name     - The name of the module's program unit.
    #   filename - The source file in which the modules was found.
    #
    def addModule( self, name, filename ):
        with self._database:
            self._database.execute( 'INSERT INTO provides VALUES ( ?, ? )', \
                                    [filename, name] )

    ###########################################################################
    # Add a dependency to the database.
    #
    # Arguments:
    #   dependor - The module which depends on a module.
    #   dependee - The module depended on by a module.
    #
    def addDependency( self, dependor, dependee ):
        with self._database:
            self._database.execute( 'INSERT INTO dependencies VALUES ( ?, ? )', \
                                    [dependor, dependee] )

    ###########################################################################
    # Get a list of dependencies for a source file.
    #
    # Arguments:
    #   filename - The source file for which dependencies are sought.
    # Return:
    #   Iterator over result map.
    #
    def getDependencySources( self, filename ):
        cursor = self._database.cursor()
        cursor.execute( 'SELECT DISTINCT p.file AS dependent_file, dp.file AS dependee_file, d.dependee AS dependee FROM dependencies AS d, provides AS p, provides as dp' \
                        + ' WHERE p.file = ? AND d.dependor = p.program_unit AND dp.program_unit = d.dependee ORDER BY d.dependee', \
                        [filename] )

        for row in cursor.fetchall():
            yield row

    ###########################################################################
    # Get the source files each program requires.
    #
    # Return:
    #   Iterator over tuples: program name, list of source filenames
    #
    # Arguments:
    #
    # TODO:
    #   Currently the recursion is handled in Python. As of 3.8.3 sqlite
    #   supports recursive queries. Once performance becomes an issue and
    #   this version of the DBMS is more prevelent a move should be made.
    #
    def getProgramSourceList( self ):
        cursor = self._database.cursor()
        cursor.execute( 'SELECT file, program_unit FROM programs, provides WHERE program_unit = name' )

        for row in cursor.fetchall():
            sources = self._recurseDependencies( row['file'], row['program_unit'] )
            uniqueSources = set( sources )
            yield row['program_unit'], list(uniqueSources)

    ###########################################################################
    # Recursively descend the dependency tree.
    #
    def _recurseDependencies( self, filename, program_unit ):
        dependencies = [filename]

        cursor = self._database.cursor()
        cursor.execute( 'SELECT file, dependee FROM dependencies, provides WHERE program_unit = dependee AND dependor = ?', \
                        [program_unit] )
        for row in cursor.fetchall():
            dependencies.extend( self._recurseDependencies( row['file'], \
                                                            row['dependee'] ) )

        return dependencies

    ###########################################################################
    # Get the source filenames for all programs.
    #
    # Return:
    #   Map of program name to source filename.
    #
    def getProgramSources( self ):
        cursor = self._database.cursor()
        cursor.execute( 'SELECT name, file FROM programs, provides WHERE program_unit = name' )

        programSource = {}
        for row in cursor.fetchall():
            programSource[row['name']] = row['file']

        return programSource

    ###########################################################################
    # Get all the dependency relationships.
    #
    # Arguments:
    # Return:
    #   Yields dependor, dependees tuples
    #
    # Return:
    #   A generator yielding (dependor, [dependee]) tuples.
    #
    def getAllFileDependencies( self ):
        cursor = self._database.cursor()
        cursor.execute( 'SELECT d.dependor AS program_unit, do.file AS dependor, de.file AS dependee' \
                        + ' FROM dependencies AS d, provides AS do, provides AS de' \
                        + ' WHERE do.program_unit = d.dependor AND de.program_unit = d.dependee' \
                        + ' ORDER BY d.dependor' )

        dependor = None
        for row in cursor.fetchall():
            if row['dependor'] != dependor:
                if dependor:
                    yield dependor, dependees
                dependor = row['dependor']
                dependees = []
            dependees.append( row['dependee'] )
        if dependor:
            yield dependor, dependees
