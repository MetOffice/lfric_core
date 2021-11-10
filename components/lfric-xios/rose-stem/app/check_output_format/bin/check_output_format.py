#!/usr/bin/env python
##############################################################################
# (c) Crown copyright 2021 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Check that diagnostic output from apps matches what is expected for UGRID
output with XIOS.
'''
import sys
import glob
import numpy as np
from netCDF4 import Dataset as ncload


def identify_mesh(filename):
    '''
    Identifies the mesh type and size based on the header of the diagnostic
    file
        filename:   a string holding the path to the file to be tested
    '''
    mesh_dict = {}
    links_dict = {}

    data = ncload(filename)

    n_mesh_faces = data.dimensions['nMesh2d_face_face'].size
    n_mesh_edges = data.dimensions['nMesh2d_edge_edge'].size
    n_mesh_nodes = data.dimensions['nMesh2d_node_node'].size

    if np.sqrt((n_mesh_faces/6)).is_integer():
        mesh_dict['type'] = "cubedsphere"
        mesh_dict['panel_length'] = int(np.sqrt((n_mesh_faces/6)))
        mesh_dict['name'] = "C{0}".format(mesh_dict['panel_length'])

        links_dict['face_face'] = n_mesh_faces
        links_dict['face_edge'] = n_mesh_faces * 2
        links_dict['face_node'] = n_mesh_faces
        links_dict['edge_edge'] = n_mesh_faces * 2
        links_dict['edge_node'] = n_mesh_faces
        links_dict['node_node'] = n_mesh_faces + 2

        mesh_dict['links'] = links_dict

    elif n_mesh_faces == n_mesh_nodes:
        mesh_dict['type'] = "planar-biperiodic"
        mesh_dict['name'] = "Biperiodic mesh of size {0}".format(n_mesh_faces)

        links_dict['face_face'] = n_mesh_faces
        links_dict['face_edge'] = n_mesh_faces * 2
        links_dict['face_node'] = n_mesh_nodes
        links_dict['edge_edge'] = n_mesh_faces * 2
        links_dict['edge_node'] = n_mesh_nodes
        links_dict['node_node'] = n_mesh_nodes

        mesh_dict['links'] = links_dict

    else:
        mesh_dict['type'] = "planar"
        mesh_dict['name'] = "Planar mesh of size {0}".format(n_mesh_faces)

        links_dict['face_edge'] = n_mesh_edges
        links_dict['face_node'] = n_mesh_nodes
        links_dict['edge_node'] = n_mesh_nodes

        mesh_dict['links'] = links_dict

    return mesh_dict


def test_header(filename, mesh_dict):
    '''
    Tests the mesh connectivity information from the NetCDF header against the
    expected values for the identified mesh
        filename:   a string holding the path to the file to be tested
        mesh_dict:  a dictionary holding information about the mesh
    '''
    header_errors = []

    mesh_links = mesh_dict['links']

    data = ncload(filename)

    for link in mesh_links:
        for dim in data.dimensions.values():
            if "Mesh2d" in dim.name:
                if "{0}".format(link) in dim.name:
                    if dim.size != mesh_links[link]:
                        header_errors.append([dim.name, dim.size,
                                              mesh_links[link]])

    return header_errors


if __name__ == "__main__":

    try:
        OUTPUT_PATH = sys.argv[1]
    except ValueError:
        sys.exit("Usage: {0} <output_path>".format(sys.argv[0]))

    TOTAL_ERRORS = 0

    DIAGNOSTIC_PATHS = glob.glob(OUTPUT_PATH + "*/*/*/*/*/*diag.nc")
    for data_path in DIAGNOSTIC_PATHS:
        print("Checking file: {0}".format(data_path))

        mesh = identify_mesh(data_path)
        errors = test_header(data_path, mesh)

        TOTAL_ERRORS = TOTAL_ERRORS + len(errors)

        if errors:
            print("\n{0} errors found in {1}:".format(len(errors), data_path))
            print("{0} mesh identified...".format(mesh['name']))

            for err in errors:
                print("{0}: expected {1}, got {2}".format(err[0], err[2],
                                                          err[1]))

    if TOTAL_ERRORS > 0:
        sys.exit("Errors found in output file mesh links - see std.out")
