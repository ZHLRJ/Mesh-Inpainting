//
// Created by haoliang on 9/11/25.
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
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> holemesh_uptr,projmesh_uptr;
std::unique_ptr<VertexPositionGeometry> holegeometry_uptr,projgeometry_uptr;
// so we can more easily pass these to different classes while preserving syntax
ManifoldSurfaceMesh* holemesh;
VertexPositionGeometry* holegeometry;
ManifoldSurfaceMesh* projmesh;
VertexPositionGeometry* projgeometry;
// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;
int main() {
    // Initialize Polyscope

    std::string filepath = "/Users/haoliang/Downloads/ongoing/compareModel/2-1-Hole.obj";
    std::string projMeshfile = "/Users/haoliang/Downloads/ongoing/compareModel/2-1-ISO.obj";
    std::tie(holemesh_uptr, holegeometry_uptr) = readManifoldSurfaceMesh(filepath);
    holemesh = holemesh_uptr.release();
    holegeometry = holegeometry_uptr.release();
    std::vector<std::vector<int>> holeVindex;
    std::vector<Eigen::MatrixXd> holeVertices;
    std::tie(holeVindex,holeVertices) = MeshHoleFind(holemesh,holegeometry);

    // Read the projection mesh
    Eigen::MatrixXd projmeshV;
    Eigen::MatrixXi projmeshF;
    std::vector<geometrycentral::Vector3> cutTrace;
    igl::readOBJ(projMeshfile, projmeshV, projmeshF);
    std::tie(projmesh_uptr, projgeometry_uptr) = makeManifoldSurfaceMeshAndGeometry(projmeshV, projmeshF);
    projmesh = projmesh_uptr.release();
    projgeometry = projgeometry_uptr.release();

    // // point projection
    // std::vector<Eigen::MatrixXd> projectedPoints;std::vector<Eigen::VectorXi> projectedFaces;
    // std::tie(projectedPoints,projectedFaces) = ProjToMesh(holeVertices,projmeshV,projmeshF);

    printf("+=============");
    polyscope::init();

    PathInfo ans,input;
    input.V = projmeshV;
    input.F = projmeshF;
    input.mesh = projmesh;
    input.geometry = projgeometry;
    drawMesh("oripoissonMesh", input.V, input.F);
    // SurfacePoint px = SurfacePoint(mesh->face(1000),geometrycentral::Vector3{0.350274, 0.277123, 0.372602});
    // SurfacePoint py = SurfacePoint(mesh->face(832),geometrycentral::Vector3{0.6, 0.3, 0.1});
    std::vector<bool> isProjPoints;
    std::vector<size_t> projPointsIdx;
    std::vector<int> OneholeVindex;
    for (auto i=0; i < 1; ++i) {
        // Process Hole i: projectedPoints[i]
        i = 0;
        Eigen::MatrixXd selectedHole = holeVertices[i];
        OneholeVindex = holeVindex[i];
        printf(" hole : %d, # of hole vertices: %d -- %d \n",i,selectedHole.rows(),OneholeVindex.size());
        // for (auto j=0; j < selectedHole.rows(); ++j)
        bool begin = true;
        size_t theFirst;
        isProjPoints.clear();
        projPointsIdx.clear();

        for (auto j=0; j < selectedHole.rows(); ++j) {
            // point projection

            Eigen::MatrixXd twoPoints(2,3);
            SurfacePoint px,py;
            twoPoints.row(0) = selectedHole.row(j);
            twoPoints.row(1) = selectedHole.row((j+1)%selectedHole.rows());
            // std::cout << "twoPoints: " << twoPoints.row(1) << std::endl;

            Eigen::MatrixXd C;Eigen::VectorXi I;
            std::tie(C,I) = PointOnMesh(twoPoints,input.V,input.F);
            // get bc of twoPoints
            Eigen::MatrixXd bc = ToBCCoordinates(input.V,input.F,C,I);
            bcProcess(bc);
            if (begin){px = SurfacePoint(projmesh->face(I(0,0)),geometrycentral::Vector3{bc(0,0), bc(0,1), bc(0,2)});}
            else{px = SurfacePoint(projmesh->vertex(input.endPoint));}
            if (j==selectedHole.rows()-1){py = projmesh->vertex(theFirst);}
            else{py = SurfacePoint(projmesh->face(I(1,0)),geometrycentral::Vector3{bc(1,0), bc(1,1), bc(1,2)});}

            // std::cout << px << std::endl;
            // std::cout << py << std::endl;
            //
            // for (auto v:py.face.adjacentVertices()){std::cout<<v.getIndex()<<std::endl;}

            printf(" Arrive segement %d; \n", j);
            GeodesicAlgorithmExact mmp(*projmesh, *projgeometry);
            mmp.propagate(py);
            std::vector<SurfacePoint> path = mmp.traceBack(px);
            // for (auto p : path) {std::cout << p << std::endl;}
            input.inputPath = path;
            ans =  updateConnection(input);
            for (auto p:ans.inputPath) {
                // std::cout << p<<std::endl;
                cutTrace.push_back(toXYZ(p,projmesh,projgeometry));
            }
            for (auto k=0; k<ans.vLink.size()-1; ++k) {
                if (k==0){isProjPoints.push_back(true);}
                else{isProjPoints.push_back(false);}
                projPointsIdx.push_back(ans.vLink[k]);
                // std::cout << isProjPoints.back() << std::endl;
                // printf("vlink # %lu, status: %d \n",projPointsIdx.back(),isProjPoints.back());
            }
            input.V = ans.V;
            input.F = ans.F;
            std::tie(projmesh_uptr, projgeometry_uptr) = makeManifoldSurfaceMeshAndGeometry(ans.V, ans.F);
            projmesh = projmesh_uptr.release();
            projgeometry = projgeometry_uptr.release();
            input.mesh = projmesh;
            input.geometry = projgeometry;
            input.endPoint = ans.endPoint;
            begin = false;
            if (j==0){theFirst = ans.start;}
            // break;
        }
    }




    // store the intermediate result
    std::string loopinfofilename = "loopinfo-2-1.txt";

    std::ofstream outputFile(loopinfofilename);
    int holeVcnt = 0;
    for (int i = 0; i < projPointsIdx.size(); i++) {
        if (isProjPoints[i]) {
            outputFile << projPointsIdx[i] << "\t" << isProjPoints[i] << "\t" << OneholeVindex[holeVcnt] << std::endl;
            holeVcnt ++;
        }
        else{   outputFile << projPointsIdx[i] << "\t" << isProjPoints[i] << "\t" << 0 << std::endl;}

    }
    outputFile.close();
    // ==================
    writeSurfaceMesh(*projmesh, *projgeometry, "holeCutting-2-1.obj");



    glm::vec3 insidecolor = glm::vec3(1.,0.,0.);
    // // Load mesh
    polyscope::CurveNetwork* trace = polyscope::registerCurveNetworkLine("cutTrace ", cutTrace);
    trace->setColor(glm::vec3(247./255,5./255,5./255));
    trace->setRadius(0.0003);

    drawMesh("ISOMesh", ans.V, ans.F);
    // polyscope::getSurfaceMesh("poissonMesh")->addFaceColorQuantity("fColor", fColor);
    polyscope::show();
}
