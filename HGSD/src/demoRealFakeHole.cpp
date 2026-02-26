//
// Created by haoliang on 10/22/25.

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
std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr;
std::unique_ptr<VertexPositionGeometry> geometry_uptr;
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

    std::string meshFile = "/Users/haoliang/Downloads/ongoing/HGSD/data/GOYLE_LOW_d0.obj";
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(meshFile);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();
    Eigen::MatrixXd meshV; Eigen::MatrixXi meshF;
    igl::readOBJ(meshFile, meshV, meshF);
    std::vector<std::array<double, 3>> fColor(meshF.rows());
    std::array<double, 3> baseColor = {241./255,245./255,255./255};
    std::array<double, 3> realHolecolor = {231./255,147./255, 251./255};
    std::array<double, 3> fakeHolecolor = {156./255,245./255, 190./255};
    std::fill(fColor.begin(), fColor.end(), baseColor);

    // 面扩展 ==========================================================================================
    // std::vector<size_t> queue{0},nextqueue;
    // std::unordered_set<size_t> visited;
    // int nlayers = 6;
    // for (int i=0;i<nlayers;i++) {
    //     while (queue.size()) {
    //         size_t Fidx = queue.back();
    //         queue.pop_back();
    //         visited.insert(Fidx);
    //         fColor[Fidx] = std::array<double, 3>{231./255,147./255, 251./255};
    //         for (auto f:mesh->face(Fidx).adjacentFaces()) {
    //             if (visited.count(f.getIndex())==0) {
    //                 nextqueue.push_back(f.getIndex());
    //                 visited.insert(f.getIndex());
    //             }
    //         }
    //     }
    //     queue = nextqueue;
    //     nextqueue.clear();
    // }
    //  ==========================================================================================
    // 点扩展

    // Define range
    int min = 0;
    int max = mesh->nVertices();
    unsigned int rseed = 11;
    // Initialize a random number generator
    std::random_device rd;
    std::mt19937 gen(rseed);
    std::uniform_int_distribution<> distrib(min, max);





    int cnthole = 70;
    for (int nh=0;nh<cnthole;nh++){
        // Generate random number in the range [min, max]

        size_t vseed = distrib(gen);
        // std::cout << "vseed = " << vseed << std::endl;
        std::vector<size_t> queue{vseed},nextqueue;
        std::unordered_set<size_t> visited;
        int nlayers = 6;
        for (int i=0;i<nlayers;i++) {
            while (queue.size()) {
                size_t vidx = queue.back();
                queue.pop_back();
                visited.insert(vidx);

                for (auto f:mesh->vertex(vidx).adjacentFaces()) {
                    size_t fIdx = f.getIndex();
                    if (nh<14){fColor[fIdx] = realHolecolor;}
                    else{fColor[fIdx] = baseColor;}

                }
                for (auto vnbr: mesh->vertex(vidx).adjacentVertices())
                    if (visited.count(vnbr.getIndex())==0) {
                        nextqueue.push_back(vnbr.getIndex());
                        visited.insert(vnbr.getIndex());
                    }
            }
            queue = nextqueue;
            nextqueue.clear();
        }

    }

    //
    //
    // for (int fileIdx = 0; fileIdx < 4; ++fileIdx) {
    //     std::string dragonFile = "/Users/haoliang/Downloads/ongoing/HGSD/data/GOYLE_LOW_d" + std::to_string(fileIdx) + ".obj";
    //     Eigen::MatrixXd dragonV; Eigen::MatrixXi dragonF;
    //     igl::readOBJ(dragonFile, dragonV, dragonF);
    //     drawMesh("Dragon" + std::to_string(fileIdx), dragonV, dragonF,231,		205,		115);
    //
    // }

    drawMesh("inputMesh", meshV, meshF,231,		205,		115);
    polyscope::getSurfaceMesh("inputMesh")->addFaceColorQuantity("fColor", fColor);
    polyscope::show();
}



