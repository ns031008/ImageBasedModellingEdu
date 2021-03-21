# workon dl4cv
# pip install open3d==0.9
import open3d as o3d
pcd = o3d.io.read_point_cloud("surface_mesh_clean.ply")
o3d.visualization.draw_geometries([pcd])

