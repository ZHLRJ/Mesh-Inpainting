//
// Created by haoliang on 10/21/25.
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

#include <unistd.h>
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr,projmesh_uptr;
std::unique_ptr<VertexPositionGeometry> geometry_uptr,projgeometry_uptr;
// so we can more easily pass these to different classes while preserving syntax
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;
ManifoldSurfaceMesh* projmesh;
VertexPositionGeometry* projgeometry;
// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;
// Custom hash struct for std::pair

int main() {
    // Initialize Polyscope
    polyscope::init();
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    for (int fileIdx = 0; fileIdx < 4; ++fileIdx) {
        std::string dragonFile = "/Users/haoliang/Downloads/ongoing/HGSD/data/GOYLE_LOW_d" + std::to_string(fileIdx) + ".obj";
        Eigen::MatrixXd dragonV; Eigen::MatrixXi dragonF;
        igl::readOBJ(dragonFile, dragonV, dragonF);
        drawMesh("Dragon" + std::to_string(fileIdx), dragonV, dragonF,231,		205,		115);

    }


    polyscope::show();
}



