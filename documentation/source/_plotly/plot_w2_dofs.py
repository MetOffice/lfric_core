# (C) British Crown Copyright 2024, Met Office.
# Please see LICENSE for license details.
"""
W2 Function space
=====================

Plot showing the dof locations for the lowest order W2 Function space.
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
                hoverinfo="none",
                hovertext=edge[0],
                name=edge[0],
                line=dict(
                    color=CUBE_EDGE_COLOUR,
                ),
            )
        )


fig = go.Figure()
add_cube(fig)
fig.update_traces(line={"width": 10})

# add the vectors - made up of a line and a cone
LINE_LENGTH = 0.3
CONE_LENGTH = 0.05
CONE_SIZEREF = 3
VECTOR_OPACITY = 0.8

# (x, y, z), (u, v, w)
# (x, y, z) start coordinates for end of arrow
# (u, v, w) vector for the arrow
vectors = [
    [(0.5, 0, 0.5), (0, -1, 0)],  # blue left
    [(0.5, 1, 0.5), (0, 1, 0)],  # blue right
    [(0, 0.5, 0.5), (-1, 0, 0)],  # green in/north
    [(1, 0.5, 0.5), (1, 0, 0)],  # green out/south
    [(0.5, 0.5, 1), (0, 0, 1)],  # red up
    [(0.5, 0.5, 0), (0, 0, -1)],  # red down
]
vector_names = ["West", "East", "North", "South", "Up", "Down"]

BLUE = f"rgba(0, 0, 255, {VECTOR_OPACITY})"
GREEN = f"rgba(0, 128, 0, {VECTOR_OPACITY})"
RED = f"rgba(255, 0, 0, {VECTOR_OPACITY})"
vector_colours = [BLUE, BLUE, GREEN, GREEN, RED, RED]

for i in range(len(vectors)):
    vector = vectors[i]
    x, y, z = vector[0]
    u, v, w = vector[1]

    line_end_x = x + u * LINE_LENGTH
    line_end_y = y + v * LINE_LENGTH
    line_end_z = z + w * LINE_LENGTH

    cone_tip_x = line_end_x + u * CONE_LENGTH
    cone_tip_y = line_end_y + v * CONE_LENGTH
    cone_tip_z = line_end_z + w * CONE_LENGTH

    fig.add_trace(
        go.Scatter3d(
            x=[x, line_end_x],
            y=[y, line_end_y],
            z=[z, line_end_z],
            mode="lines",
            line=dict(
                color=vector_colours[i],
                width=10,
            ),
            name=vector_names[i],
        )
    )
    fig.add_trace(
        go.Cone(
            x=[cone_tip_x],
            y=[cone_tip_y],
            z=[cone_tip_z],
            u=[u * CONE_LENGTH],
            v=[v * CONE_LENGTH],
            w=[w * CONE_LENGTH],
            anchor="tip",
            colorscale=[[0, vector_colours[i]], [1, vector_colours[i]]],
            name=vector_names[i],
            showscale=False,
            sizeref=CONE_SIZEREF,
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
    hovermode=False,
    scene_camera_eye=dict(x=-1.2, y=1.8, z=0.8),
    scene=dict(
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        zaxis=dict(visible=False),
        aspectmode="manual",
        aspectratio=go.layout.scene.Aspectratio(x=2, y=2, z=2),
    ),  # aspect set to 2 times zoom
    autosize=False,
    width=700,
    height=500,
    margin=dict(l=50, r=50, b=50, t=50, pad=4),
    paper_bgcolor="rgba(250, 249, 246, 0.9)",  # off white
)

output_html_path = r"generated/html/plots/plot_w2_dofs.html"
fig.write_html(output_html_path, include_plotlyjs=False, full_html=False)
