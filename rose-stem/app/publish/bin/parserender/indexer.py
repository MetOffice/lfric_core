#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
'''
Examine the output directory from a nightly suite run and determine things
about its contents.
'''
from __future__ import print_function

from abc import ABCMeta
import datetime
import glob
import hashlib
import os
import os.path
import xml.etree.ElementTree as et

##############################################################################
class PageDetails(object):
    '''
    Holds details about a page file.
    '''
    def __init__( self, url, name, timestamp, digest ):
        self.url       = url
        self.name      = name
        self.timestamp = timestamp
        self.digest    = digest

##############################################################################
class SubsiteDetails(object):
    '''
    Holds details about a subsite directory.
    '''
    def __init__( self, url, name, timestamp ):
        self.url       = url
        self.name      = name
        self.timestamp = timestamp

##############################################################################
class Indexer(object):
    '''
    Indexers inherit from this class.
    '''
    __MetaClass__ = ABCMeta

##############################################################################
class DynamoIndexer(Indexer):
    '''
    Examine a Dynamo site.
    '''
    def __init__( self, pathname ):
        self.subsites = {}

        self._rootPath = pathname

    def examine( self ):
        '''
        Whizz through the directory finding out what's in it and commenting
        on the result.
        '''
        cronLogFilename = os.path.join( self._rootPath, 'cron.out' )
        if os.path.exists( cronLogFilename ):
            self.cronOut = 'cron.out'
            timestamp = os.path.getmtime( cronLogFilename )
            self.cronTimestamp = datetime.datetime.utcfromtimestamp( timestamp )
        else:
            self.cronOut = None
            self.cronTimestamp = None
        self.subsites = self.findDocumentation( self._rootPath )
        self.pages = self.findPages( self._rootPath )

    def findDocumentation( self, rootPath ):
        '''
        Walks through the directory tree from the bottom upwards, examining
        the contents of each directory. If an index.html file is found
        it generates a SubsiteDetails object linking to that. If no
        index.html is found then the directory is linked to instead
        '''
        subsites = {}
        for path, directorynames, filenames in os.walk( rootPath, topdown=False ):

            relativePath = os.path.relpath( path, rootPath )

            if relativePath != '.': relativePath = relativePath + os.sep


            if 'index.html' in filenames and relativePath != '.':
                # Found an index.html file 
                # Determine the timestamp of index.html
                timestamp = os.path.getmtime( os.path.join( path, \
                                                            'index.html' ) )
                timestamp = datetime.datetime.utcfromtimestamp( timestamp )
                # Add it to the list of subsites
                subsites[relativePath] = SubsiteDetails( relativePath, \
                                relativePath.replace( os.sep, ' ' ).title(), \
                                                         timestamp )

            else:
                if (len(directorynames) > 0):
                    if(relativePath != '.'):
                      for d in directorynames:
                        # Need to check that we are not duplicating top level of paths
                        # that have index.html and have already been added
                        if not(self.siteExists(relativePath, subsites)):
                          timestamp = os.path.getmtime( os.path.join( path, d ) )
                          timestamp = datetime.datetime.utcfromtimestamp( timestamp )
                          subsites[relativePath] = SubsiteDetails( relativePath, \
                                     relativePath.replace( os.sep, ' ' ).title(), \
                                                                      timestamp )

        return subsites

   
    def siteExists(self, relPath, subsites):
       # Checks if a path is a subpath of an entry in the current subsites list
       exists = False
       # Check if relPath is a substring of each entry in subsites
       for s in subsites:
         if relPath in s:
           exists = True
       return exists

    def findPages( self, rootPath ):
        '''
        Looks for HTML files which are not "index.html" and pulls details from
        them. Basically it is assumed that they are compile or run reports
        generated by other scripts in this package.
        '''
        pages = {}
        for filename in glob.iglob( '{}/*.html'.format( rootPath ) ):
            relativeFilename = os.path.relpath( filename, rootPath )
            if relativeFilename == 'index.html': continue
            leafname = os.path.basename( filename )
            document = et.parse( filename )
            compiler = document.getroot().find( './/span[@id=\'compiler\']' ).text
            contextNode  = document.getroot().find( './/span[@id=\'context\']' )
            context = contextNode.text if contextNode is not None else None
            timestampNode = document.getroot().find( './/span[@id=\'timestamp\']' )
            timestamp = datetime.datetime.strptime( timestampNode.text, \
                                                    '%Y-%m-%dT%H:%M:%SZ' )

            hasher = hashlib.sha1()
            hasher.update( compiler )
            if context: hasher.update( context )
            for eventsNode in document.getroot().find( './/div[@id=\'events\']' ).iter():
                for node in eventsNode:
                    if node.text and node.text.strip() != '':
                        hasher.update( node.text.strip() )

            name, extension = leafname.rsplit( '.', 1 )
            pages[relativeFilename] = PageDetails( relativeFilename,         \
                                           name.replace( '.', ' ' ).title(), \
                                                   timestamp,                \
                                                   hasher.hexdigest() )

        return pages
