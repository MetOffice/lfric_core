# (C) British Crown Copyright 2024, Met Office.
# Please see LICENSE for license details.
"""
Wtheta Function space
=====================

Plot showing the dof locations for the lowest order Wtheta Function space.
"""

import plotly.graph_objects as go
import numpy as np

# sphinx_gallery_thumbnail_path = '_static/fs_thumbnails/wtheta_k0.png'


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
DOF_OPACITY = 0.75

dof_x = [0.5, 0.5]
dof_y = [0.5, 0.5]
dof_z = [0.0, 1.0]
dof_names = ["dof 1", "dof 2"]

for i in range(len(dof_x)):
    fig.add_trace(
        go.Scatter3d(
            x=[dof_x[i]],
            y=[dof_y[i]],
            z=[dof_z[i]],
            mode="markers",
            marker=dict(
                symbol="x",
                color=np.zeros(1),
                colorscale="picnic",
                cmin=0,  # the cmin and cmax make sure only the start colour of the map is used
                cmax=10000,
                opacity=DOF_OPACITY,
            ),
            name=dof_names[i],
            surfaceaxis=1,
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
# fig.show() # to test changes locally
fig  # correct output for sphinx-gallery
