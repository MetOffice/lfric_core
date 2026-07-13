import glob
import os
import xml.etree.ElementTree as ET

import pytest


def srcReplace(parent, path):
    """
    Recursively update input parent element, resolVing `src` content and replacing this
    with parsed XML elements from referenced files.
    """
    i = 0
    for elem in parent:
        if elem.attrib.get('src'):
            inf = os.path.join(path, elem.attrib['src'])
            print(f'parsing: {inf}')
            newelem = ET.parse(inf).getroot()
            parent[i] = newelem
        else:
            srcReplace(elem, path)
        i += 1


def load_source_xml(fname):
    path = os.path.dirname(fname)
    tree = ET.parse(fname)
    root = tree.getroot()

    # Validate that external elements, defined using `src=` are only accessing known internal paths
    # to minimise rish from "XML external entity attack"
    for elem in root.findall('.//*[@src]'):
        if not (elem.attrib['src'].startswith('metadata/') or
                elem.attrib['src'].startswith('etc/') or
                elem.attrib['src'].startswith('./') or
                elem.attrib['src'].startswith('../') or
                elem.attrib['src'].startswith('$SOURCE_ROOT/')):
            raise ValueError('only `src` attributes from local `metadata/` or '
                             '`etc/` or `$SOURCE_ROOT/`paths are supported, not '
                             f'{elem.attrib["src"]}')
    srcReplace(root, path)
    return tree

root_dir = os.environ.get('$SOURCE_ROOT', '/home/users/mark.hedley/metofficegit/') + 'lfric_core'

iodef_likes = []
for iodef_like in glob.glob('**/iodef*.xml',
                            root_dir=root_dir,
                            recursive=True):
    infile = os.path.join(root_dir, iodef_like)
    tree = load_source_xml(infile)
    root = tree.getroot()
    iodef_likes.append((root, infile))

# tree = load_source_xml(os.path.join(hard_path, iodef))
# root = tree.getroot()

@pytest.mark.parametrize("aroot, infile", iodef_likes)
def test_unique_field_ids_within_context(aroot, infile):
    print(f'validating XIOS XML for {aroot} ...')

    for context in aroot.findall('.//context'):
        field_ids = set()    
        for elem in context.findall('.//field[@id]'):
    # for elem in aroot.findall('.//field[@id]'):
            eid = elem.attrib.get('id')
            if eid is not None:
                err_str = (f'Within {infile}\n:'
                           '\tThe `field` element:\n'
                           f'{elem.attrib}\n'
                           f'\thas an `id`: "{eid}" '
                           'which is already defined for an existing field.\n'
                           '\tXIOS will not distinguish betwen these fields, operations '
                           'using this id will affect all fields with this `id`')
                assert eid not in field_ids, err_str
                field_ids.add(eid)
