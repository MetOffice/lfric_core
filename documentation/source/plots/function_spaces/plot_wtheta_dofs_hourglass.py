# (C) British Crown Copyright 2024, Met Office.
# Please see LICENSE for license details.
"""
Wtheta Function space
=====================

Plot showing the dof locations for the lowest order Wtheta Function space.
With cones to an hourglass?
"""

import plotly.graph_objects as go
import numpy as np

# sphinx_gallery_thumbnail_path = '_static/fs_thumbnails/wtheta_k0_hourglass.svg'


CUBE_EDGE_COLOUR = "darkgrey"


def add_cube(fig):
    zeros = np.zeros(2)
    ones = np.ones(2)
    zero_to_one = np.linspace(0, 1, 2, endpoint=True)

    cube_edges = []

    # bottom face edges
    cube_edges.append(["edge_1", [zero_to_one, zeros, zeros]])  # 0,0,0 to 1,0,0
    cube_edges.append(["edge_2", [ones, zero_to_one, zeros]])  # 1,0,0 to 1,1,0
    cube_edges.append(["edge_3", [zero_to_one, ones, zeros]])  # 0,1,0 to 1,1,0
    cube_edges.append(["edge_4", [zeros, zero_to_one, zeros]])  # 0,0,0 to 0,1,0

    # sides
    cube_edges.append(["edge_5", [zeros, zeros, zero_to_one]])  # 0,0,0 to 0,0,1
    cube_edges.append(["edge_6", [ones, zeros, zero_to_one]])  # 1,0,0 to 1,0,1
    cube_edges.append(["edge_7", [ones, ones, zero_to_one]])  # 1,1,0 to 1,1,1
    cube_edges.append(["edge_8", [zeros, ones, zero_to_one]])  # 0,1,0 to 0,1,1

    # top face edges
    cube_edges.append(["edge_9", [zero_to_one, zeros, ones]])  # 0,0,1 to 1,0,1
    cube_edges.append(["edge_10", [ones, zero_to_one, ones]])  # 1,0,1 to 1,1,1
    cube_edges.append(["edge_11", [zero_to_one, ones, ones]])  # 0,1,1 to 1,1,1
    cube_edges.append(["edge_12", [zeros, zero_to_one, ones]])  # 0,0,1 to 0,1,1

    for edge in cube_edges:
        fig.add_trace(
            go.Scatter3d(
                x=edge[1][0],
                y=edge[1][1],
                z=edge[1][2],
                mode="lines",
                hovertext=edge[0],
                name=edge[0],
                line=dict(
                    color=CUBE_EDGE_COLOUR,
                ),
            )
        )
        

# create figure with the cube
fig = go.Figure()
add_cube(fig)
fig.update_traces(line={"width": 10})

# add the dofs
dof_x = [0.5, 0.5]
dof_y = [0.5, 0.5]
dof_z = [0.0, 1.0]
dof_names = ["dof 1", "dof 2"]

BLUE = f"rgba(0, 0, 255, 1)"
LIGHTBLUE = f"rgba(100, 100, 255, 1)"
CONE_LENGTH = 0.3
CONE_SIZEREF = 1
CONE_OVERLAP = 0.5 # more than .5 and tips come out top
for i in range(len(dof_x)):
    
    u, v, w = (0, 0, 0)
    
    bottom_cone_tip_x = dof_x[i]
    top_cone_tip_x = dof_x[i]
    bottom_cone_tip_y = dof_y[i]
    top_cone_tip_y = dof_y[i]
    bottom_cone_tip_z = dof_z[i]
    top_cone_tip_z = dof_z[i]
    
    if dof_x[i] == dof_y[i]:
        # z
        w = 1
        bottom_cone_tip_z = dof_z[i] + CONE_LENGTH*CONE_OVERLAP
        top_cone_tip_z = dof_z[i] - CONE_LENGTH*CONE_OVERLAP
    if dof_x[i] == dof_z[i]:
        # y
        v = 1
        bottom_cone_tip_y = dof_y[i] + CONE_LENGTH*CONE_OVERLAP
        top_cone_tip_y = dof_y[i] - CONE_LENGTH*CONE_OVERLAP
    if dof_z[i] == dof_y[i]:
        # x
        u = 1
        bottom_cone_tip_x = dof_x[i] + CONE_LENGTH*CONE_OVERLAP
        top_cone_tip_x = dof_x[i] - CONE_LENGTH*CONE_OVERLAP
    
    fig.add_trace(
        go.Cone(
            x=[bottom_cone_tip_x],
            y=[bottom_cone_tip_y],
            z=[bottom_cone_tip_z],
            u=[u * CONE_LENGTH],
            v=[v * CONE_LENGTH],
            w=[w * CONE_LENGTH],
            anchor="tip",
            colorscale=[[0, BLUE], [1, BLUE]],
            showscale=False,
            sizemode="scaled",
            sizeref=CONE_SIZEREF,
            cmax=0,
            cmin=-1*CONE_LENGTH,
            lighting_diffuse=0.3,
        ))
    fig.add_trace(
        go.Cone(
            x=[top_cone_tip_x],
            y=[top_cone_tip_y],
            z=[top_cone_tip_z],
            u=[-1 * u * CONE_LENGTH],
            v=[-1 * v * CONE_LENGTH],
            w=[-1 * w * CONE_LENGTH],
            anchor="tip",
            colorscale=[[0, LIGHTBLUE], [1, LIGHTBLUE]],
            showscale=False,
            sizemode="scaled",
            sizeref=CONE_SIZEREF,
            cmax=1*CONE_LENGTH,
            cmin=0,
            lighting_diffuse=0.3,
        )
    )


# add any shaded planes
FACE_OPACITY = 0.1
TOP_FACE = [[1, 1], [1, 1]]
fig.add_trace(
    go.Surface(
        z=TOP_FACE,
        showscale=False,
        opacity=FACE_OPACITY,
        colorscale=[[0, "grey"], [1, "grey"]],
        name="Top Face",
    )
)
BOTTOM_FACE = [[0, 0], [0, 0]]
fig.add_trace(
    go.Surface(
        z=BOTTOM_FACE,
        showscale=False,
        opacity=FACE_OPACITY,
        colorscale=[[0, "grey"], [1, "grey"]],
        name="Bottom Face",
    )
)

# turn off legend, axis and show
fig.update_layout(
    showlegend=False,
    scene_camera_eye=dict(x=-1.2, y=1.8, z=0.8),
    scene=dict(
        xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False)
    ),
    autosize=False,
    width=700,
    height=500,
    margin=dict(l=50, r=50, b=50, t=50, pad=4),
    paper_bgcolor="rgba(250, 249, 246, 0.9)",
)

# Local only
# fig.show() # to test changes locally
# fig.write_image("source/_static/fs_thumbnails/wtheta_k0_hourglass.svg")  # for the thumbnail

# Correct output for sphinx-gallery
fig
