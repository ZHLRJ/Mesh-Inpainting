//
// Created by haoliang on 10/28/25.
//
//
// Created by haoliang on 10/16/25.
//

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#include "geometrycentral/surface/signed_heat_method.h"
#include "geometrycentral/pointcloud/point_cloud_heat_solver.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/utilities/vector3.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "igl/readOBJ.h"
#include "igl/point_mesh_squared_distance.h"
#include "igl/AABB.h"
#include "igl/barycentric_coordinates.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/exact_geodesics.h"
#include "polyscope/curve_network.h"
#include <Vector>
#include <unordered_set>
#include <unordered_map>
#include "utils.h"
#include <fstream>
#include <iostream>
#include <format>
#include <unistd.h>
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> tmpmesh_uptr;
std::unique_ptr<VertexPositionGeometry> tmpgeometry_uptr;
// so we can more easily pass these to different classes while preserving syntax

ManifoldSurfaceMesh* tmpmesh;
VertexPositionGeometry* tmpgeometry;
// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;


int main() {
    // Initialize Polyscope
    polyscope::init();

    // Read the mesh
    Eigen::MatrixXd meshV;
    Eigen::MatrixXi meshF;
    std::string fileName = "/Users/haoliang/Downloads/ongoing/compareModel/result/4-1-ISO.obj";

    igl::readOBJ(fileName, meshV, meshF);

    std::unique_ptr<SurfaceMesh> mesh_uptr;
    std::unique_ptr<VertexPositionGeometry> geometry_uptr;
    SurfaceMesh* mesh;
    VertexPositionGeometry* geometry;


    std::tie(mesh_uptr, geometry_uptr) = makeSurfaceMeshAndGeometry(meshV, meshF);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();
    // std::tie(mesh_uptr, geometry_uptr) = makeManifoldSurfaceMeshAndGeometry(meshV, meshF);
    // mesh = mesh_uptr.release();
    // geometry = geometry_uptr.release();
    // Register the mesh with Polyscope
    polyscope::registerSurfaceMesh("input_mesh", meshV, meshF);
    // writeSurfaceMesh(*mesh, *geometry, "thisistest.obj");

    // Show the GUI
    polyscope::show();
}
