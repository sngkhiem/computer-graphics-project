import numpy as np
import matplotlib.pyplot as plt

# File path for input and output
file_path = "D:/CG/log.txt"
output_path = "D:/CG/after.png"

# Read the input file
with open(file_path, "r") as file:
    lines = file.readlines()

# Parse the input data
numbFaces = int(lines[0].strip())
faces = np.array([list(map(int, lines[i].strip().split())) for i in range(1, 1 + numbFaces)])
num_vertices = int(lines[1 + numbFaces].strip())
source_vertex = int(lines[2 + numbFaces].strip())
distances = np.array([float(lines[i].strip()) for i in range(3 + numbFaces, 3 + numbFaces + num_vertices)])
vertices = np.array([
    list(map(float, lines[i].strip().split())) 
    for i in range(3 + numbFaces + num_vertices, 3 + numbFaces + 2 * num_vertices)
])

# Assertions for consistency
assert faces.shape[1] == 3, "Each face must have exactly 3 vertices"
assert len(distances) == num_vertices, "Distances array length must match the number of vertices"
assert vertices.shape[1] == 3, "Each vertex must have exactly 3 coordinates"

# Prepare data for plotting
x = vertices[:, 0]
y = distances
z = vertices[:, 2]

# Create the plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 12))

# First plot: Geodesic distance
contour = ax1.tricontourf(x, z, faces, y, levels=20, cmap="viridis")
cbar = fig.colorbar(contour, ax=ax1, orientation="horizontal", label="Geodesic Distance")
ax1.set_title("Geodesic Distance")
ax1.set_xlabel("X")
ax1.set_ylabel("Z")
ax1.set_aspect("equal")

# Calculate angles of triangles
v0 = vertices[faces[:, 0]]
v1 = vertices[faces[:, 1]]
v2 = vertices[faces[:, 2]]
a = np.linalg.norm(v1 - v2, axis=1)
b = np.linalg.norm(v2 - v0, axis=1)
c = np.linalg.norm(v0 - v1, axis=1)

angle_0 = np.arccos((b**2 + c**2 - a**2) / (2 * b * c))
angle_1 = np.arccos((a**2 + c**2 - b**2) / (2 * a * c))
angle_2 = np.arccos((a**2 + b**2 - c**2) / (2 * a * b))
angles = np.degrees(np.stack((angle_0, angle_1, angle_2), axis=1))

# Second plot: Triangle angles
for i in range(3):  
    ax2.bar(
        np.arange(len(angles)) + i * 0.25,  
        angles[:, i],
        width=0.25,
        label=f"Angle {i + 1}"
    )
ax2.set_title("Triangle Angles for Each Face")
ax2.set_xlabel("Face Index")
ax2.set_ylabel("Angle (Degrees)")
ax2.legend()

# Save the plot as an image
plt.savefig(output_path, dpi=300)
plt.close()

# Calculate angle statistics
max_angle = np.max(angles)
min_angle = np.min(angles)
mean_angle = np.mean(angles)

