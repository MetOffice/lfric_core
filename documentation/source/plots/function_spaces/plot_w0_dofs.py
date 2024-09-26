# (C) British Crown Copyright 2024, Met Office.
# Please see LICENSE for license details.
"""
W0 Function space
=====================

Plot showing the dof locations for the lowest order W0 Function space.
"""

import plotly.graph_objects as go
import numpy as np

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
dof_x = [0, 1, 1, 0, 0, 1, 1, 0]
dof_y = [0, 0, 1, 1, 0, 0, 1, 1]
dof_z = [0, 0, 0, 0, 1, 1, 1, 1]
dof_names = ["dof 1", "dof 2"]

LINE_LENGTH = 0.1
CROSS_OPACITY = 0.5
BLUE = f"rgba(0, 0, 255, {CROSS_OPACITY})"

for i in range(len(dof_x)):
    l1_start_x = dof_x[i] + LINE_LENGTH / 2
    l1_start_y = dof_y[i] + LINE_LENGTH / 2
    l1_start_z = dof_z[i] + LINE_LENGTH / 2
    l1_end_x = dof_x[i] - LINE_LENGTH / 2
    l1_end_y = dof_y[i] - LINE_LENGTH / 2
    l1_end_z = dof_z[i] - LINE_LENGTH / 2

    l2_start_x = dof_x[i] - LINE_LENGTH / 2
    l2_start_y = dof_y[i] + LINE_LENGTH / 2
    l2_start_z = dof_z[i] + LINE_LENGTH / 2
    l2_end_x = dof_x[i] + LINE_LENGTH / 2
    l2_end_y = dof_y[i] - LINE_LENGTH / 2
    l2_end_z = dof_z[i] - LINE_LENGTH / 2

    l3_start_x = dof_x[i] + LINE_LENGTH / 2
    l3_start_y = dof_y[i] - LINE_LENGTH / 2
    l3_start_z = dof_z[i] + LINE_LENGTH / 2
    l3_end_x = dof_x[i] - LINE_LENGTH / 2
    l3_end_y = dof_y[i] + LINE_LENGTH / 2
    l3_end_z = dof_z[i] - LINE_LENGTH / 2

    l4_start_x = dof_x[i] - LINE_LENGTH / 2
    l4_start_y = dof_y[i] - LINE_LENGTH / 2
    l4_start_z = dof_z[i] + LINE_LENGTH / 2
    l4_end_x = dof_x[i] + LINE_LENGTH / 2
    l4_end_y = dof_y[i] + LINE_LENGTH / 2
    l4_end_z = dof_z[i] - LINE_LENGTH / 2

    fig.add_trace(
        go.Scatter3d(
            x=[l1_start_x, l1_end_x],
            y=[l1_start_y, l1_end_y],
            z=[l1_start_z, l1_end_z],
            mode="lines",
            line=dict(
                color=BLUE,
                width=10,
            ),
        )
    )
    fig.add_trace(
        go.Scatter3d(
            x=[l2_start_x, l2_end_x],
            y=[l2_start_y, l2_end_y],
            z=[l2_start_z, l2_end_z],
            mode="lines",
            line=dict(
                color=BLUE,
                width=10,
            ),
        )
    )
    fig.add_trace(
        go.Scatter3d(
            x=[l3_start_x, l3_end_x],
            y=[l3_start_y, l3_end_y],
            z=[l3_start_z, l3_end_z],
            mode="lines",
            line=dict(
                color=BLUE,
                width=10,
            ),
        )
    )
    fig.add_trace(
        go.Scatter3d(
            x=[l4_start_x, l4_end_x],
            y=[l4_start_y, l4_end_y],
            z=[l4_start_z, l4_end_z],
            mode="lines",
            line=dict(
                color=BLUE,
                width=10,
            ),
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

# turn off hover tooltips
fig.update_traces(hoverinfo="skip", hovertemplate=None)

fig.update_scenes(
    xaxis_showspikes=False, yaxis_showspikes=False, zaxis_showspikes=False
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

fig.show()
output_html_path = r"./html/plot_w0_dofs.html"
fig.write_html(output_html_path, include_plotlyjs=False, full_html=False)
