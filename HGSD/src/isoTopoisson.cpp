//
// Created by haoliang on 8/27/25.
//
//
// Created by haoliang on 8/1/25.
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
#include "geometrycentral/surface/flip_geodesics.h"
#include "geometrycentral/surface/mesh_graph_algorithms.h"
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
#include "cstdio"
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr,isomesh_uptr;
std::unique_ptr<VertexPositionGeometry> geometry_uptr,isogeometry_uptr;
// so we can more easily pass these to different classes while preserving syntax
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;
ManifoldSurfaceMesh* poissonmesh;
VertexPositionGeometry* poissongeometry;
// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;

polyscope::SurfaceMesh * loadmesh(std::string& filename,std::string meshName) {
    Eigen::MatrixXd isomeshV;
    Eigen::MatrixXi isomeshF;
    igl::readOBJ(filename, isomeshV, isomeshF);
    polyscope::SurfaceMesh * tmpmesh = polyscope::registerSurfaceMesh(meshName, isomeshV, isomeshF);
    return tmpmesh;
}
int edgewidth = 0.45;
int main(int argc, char **argv) {

    // Initialize Polyscope
    polyscope::init();
    std::string gtMeshfile = "/Users/haoliang/Downloads/ongoing/signed-heat-demo-3d/data/zhldataset/bunny_hole.obj";

    std::string isoMeshfile = "/Users/haoliang/Downloads/ongoing/result/isosurface.obj";
    std::string poissonDirectlyFile = "/Users/haoliang/Downloads/ongoing/result/poissonDirectly.obj";
    std::string poissonFinalFile = "/Users/haoliang/Downloads/ongoing/result/poissonMesh.obj";

    // Eigen::RowVector3d shift(0,0,0);
    polyscope::SurfaceMesh * gtMesh = loadmesh(gtMeshfile,"gtMesh");
    gtMesh->setSurfaceColor(glm::vec3(241.0/255,128.0/255,223.0/255));
    gtMesh->setEdgeWidth(0.45);
    polyscope::SurfaceMesh * isoMesh = loadmesh(isoMeshfile,"isoMesh");
    isoMesh->setSurfaceColor(glm::vec3(81.0/255,212.0/255,131.0/255));
    isoMesh->setEdgeWidth(0.45);
    polyscope::SurfaceMesh * poissonDense = loadmesh(poissonDirectlyFile,"poissonDense");
    poissonDense->setSurfaceColor(glm::vec3(231.0/255,205.0/255,115.0/255));
    poissonDense->setEdgeWidth(0.45);
    polyscope::SurfaceMesh * poissonFinal = loadmesh(poissonFinalFile,"poissonFinal");
    poissonFinal->setSurfaceColor(glm::vec3(231.0/255,205.0/255,115.0/255));
    poissonFinal->setEdgeWidth(0.45);

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::readOBJ(poissonDirectlyFile, V, F);


    polyscope::PointCloud* poissonpointCloud = polyscope::registerPointCloud("gtVertexPoints", V);
    // holeVerticespoints->setPointColor(glm::vec3(249./255,0.0/255,0.0/255));
    poissonpointCloud->setPointColor(glm::vec3(81.0/255,212.0/255,131.0/255));
    // poissonpointCloud->setPointRadius(0.3);
    poissonpointCloud->setEnabled(false);
    // set soft shadows on the ground
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    polyscope::options::ssaaFactor = 2;
    polyscope::show();
    }