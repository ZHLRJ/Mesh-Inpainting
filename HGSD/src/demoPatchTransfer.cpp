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
    Eigen::MatrixXd patchV,originalMeshV; Eigen::MatrixXi patchF,originalMeshF;
    std::string originalMeshFile = "/Users/haoliang/Downloads/ongoing/signed-heat-demo-3d/data/zhldataset/bunny_hole.obj";
    igl::readOBJ(originalMeshFile, originalMeshV, originalMeshF);
    // cut patch from given information

    for (int fileIdx = 0; fileIdx < 5; ++fileIdx) {
        std::string patchMeshFile = "/Users/haoliang/Downloads/ongoing/HGSD/bunnyCut/patch" + std::to_string(fileIdx) + ".obj";
        std::string inputMeshCutInfoFile = "/Users/haoliang/Downloads/ongoing/HGSD/bunnyCut/cutInfo_" + std::to_string(fileIdx) + ".txt";
        igl::readOBJ(patchMeshFile, patchV, patchF);
        drawMesh("Patch" + std::to_string(fileIdx), patchV, patchF,156,245,190);

        std::ifstream inputFile(inputMeshCutInfoFile);
        std::vector<std::array<size_t,3>> inputInfo;
        size_t a,b,c;
        while (inputFile >> a >> b>> c) {
            inputInfo.push_back({a,b,c});
        }
        inputFile.close(); // Close the file stream


        // while (inputFile >> val >> idx>> holeidx) {
        //     isProjPoints.push_back(idx); inputInfo[0]
        //     projPointsIdx.push_back(val); inputInfo[1]
        //     HoleboundaryPointsIdx.push_back(holeidx); inputInfo[2]



    }

    // std::string refMeshFile = "/Users/haoliang/Downloads/ongoing/HGSD/bunnyCut/targetMesh_" + std::to_string(0) + ".obj";
    // Eigen::MatrixXd refMeshV; Eigen::MatrixXi refMeshF;
    // igl::readOBJ(refMeshFile, refMeshV, refMeshF);
    // drawMesh("refMeshV", refMeshV, refMeshF,156,245,190	);

    polyscope::show();
}



