#!/usr/bin/env python3
##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
XIOS Iodef.xml metadata consistency validation tests.
"""

import glob
import os
import re
import xml.etree.ElementTree as ET

import pytest


# security pattern to check whether `src` links are local and link to known
# controlled facets of the local source tree.  Do not load from unknown sources.
# allows strings that start with:
# 'metadata/' etc/' './' '../metadata' '$SOURCE_ROOT/'
# followed by alphabetic characters only
security_re = re.compile(r'^(metadata/|etc/|\./|../metadata/|\$SOURCE_ROOT/)'
                         r'[a-zA-Z]+\w')

def src_replace(parent, path):
    """
    Recursively update input parent element, resolVing `src` content and
    replacing this with parsed XML elements from referenced files.
    """
    i = 0
    for elem in parent:
        if elem.attrib.get('src'):
            inf = os.path.join(path, elem.attrib['src'])
            print(f'parsing: {inf}')
            newelem = ET.parse(inf).getroot()
            parent[i] = newelem
        else:
            src_replace(elem, path)
        i += 1


def load_source_xml(fname):
    """
    Load a source XML file and recursively populated `src` links from
    identified safe local paths.
    """
    path = os.path.dirname(fname)
    load_tree = ET.parse(fname)
    load_root = load_tree.getroot()

    # Validate that external elements, defined using `src=` are only accessing
    # known internal paths to minimise rish from "XML external entity attack"
    for elem in load_root.findall('.//*[@src]'):
        if security_re.match(elem.attrib['src']) is None:
            raise ValueError('only `src` attributes from local, `metadata/`, '
                             '`etc/` or `$SOURCE_ROOT/`paths are supported, not '
                             f'{elem.attrib["src"]}')
    src_replace(load_root, path)
    return load_tree

# Generator for the pytest parametrize fixture.
root_dir = os.environ.get('SOURCE_PATH', '')
print(f'root_dir = {root_dir}')
iodef_likes = []
for iodef_like in glob.glob('**/iodef*.xml',
                            root_dir=root_dir,
                            recursive=True):
    infile = os.path.join(root_dir, iodef_like)
    tree = load_source_xml(infile)
    root = tree.getroot()
    iodef_likes.append((root, infile))

if len(iodef_likes) == 0:
    raise ValueError('No `iodef*.xml` files found.')

@pytest.mark.parametrize("aroot, test_file", iodef_likes)
def test_unique_field_ids_within_context(aroot, test_file):
    """
    Pytest test to verify that `field` element `id`s are unique within
    a given XIOS context.
    """
    print(f'validating XIOS XML for {aroot} ...')

    for context in aroot.findall('.//context'):
        field_ids = set()
        for elem in context.findall('.//field[@id]'):
            eid = elem.attrib.get('id')
            if eid is not None:
                err_str = (f'Within {test_file}\n:'
                           '\tThe `field` element:\n'
                           f'{elem.attrib}\n'
                           f'\thas an `id`: "{eid}" '
                           'which is already defined for an existing field.\n'
                           '\tXIOS will not distinguish betwen these fields, '
                           'operations using this id will affect all fields '
                           'with this `id` leading to indeterminate behaviour')
                assert eid not in field_ids, err_str
                field_ids.add(eid)
