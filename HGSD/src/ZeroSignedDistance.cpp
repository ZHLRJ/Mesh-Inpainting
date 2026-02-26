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
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr,isomesh_uptr;
std::unique_ptr<VertexPositionGeometry> geometry_uptr,isogeometry_uptr;
// so we can more easily pass these to different classes while preserving syntax
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;
ManifoldSurfaceMesh* isomesh;
VertexPositionGeometry* isogeometry;
// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;
std::vector<std::vector<int>> holeVerticesFind(ManifoldSurfaceMesh* mesh) {
    std::unordered_map<int, int> halfEdgeMap;
    for (Halfedge e :mesh->exteriorHalfedges()) {
        halfEdgeMap[e.sibling().vertex().getIndex()] = e.vertex().getIndex();
    }
    std::unordered_set<int> isVisited;
    std::vector<int> curLoop;
    std::vector<std::vector<int>> holeVertices;
    // std::vector<std::vector<Eigen::Vector3d>> holeVertices;
    for (std::pair<int,int> v: halfEdgeMap) {
        if (isVisited.find(v.first) == isVisited.end()) {
            curLoop.push_back(v.second);
            while (isVisited.find(halfEdgeMap[curLoop.back()]) == isVisited.end()){
                isVisited.insert(halfEdgeMap[curLoop.back()]);
                curLoop.push_back(halfEdgeMap[curLoop.back()]);
            }
            holeVertices.push_back(curLoop);
            curLoop.clear();
        }
    }
    return holeVertices;
}
std::vector<Eigen::MatrixXd> holeVerticesPosFind(ManifoldSurfaceMesh* mesh,VertexPositionGeometry* geometry) {
    std::unordered_map<int, int> halfEdgeMap;
    for (Halfedge e :mesh->exteriorHalfedges()) {
        halfEdgeMap[e.sibling().vertex().getIndex()] = e.vertex().getIndex();
    }
    std::unordered_set<int> isVisited;
    std::vector<int> curLoop;
    std::vector<Eigen::MatrixXd> holeVertices;
    for (std::pair<int,int> v: halfEdgeMap) {
        if (isVisited.find(v.first) == isVisited.end()) {
            curLoop.push_back(v.second);
            while (isVisited.find(halfEdgeMap[curLoop.back()]) == isVisited.end()){
                isVisited.insert(halfEdgeMap[curLoop.back()]);
                curLoop.push_back(halfEdgeMap[curLoop.back()]);
            }
            Eigen::MatrixXd curLoopPosMatrix(curLoop.size(),3);
            for (int i = 0; i < curLoop.size(); i++) {
                curLoopPosMatrix.row(i) <<geometry->inputVertexPositions[curLoop[i]].x, geometry->inputVertexPositions[curLoop[i]].y,geometry->inputVertexPositions[curLoop[i]].z;
            }
            holeVertices.push_back(curLoopPosMatrix);
            curLoop.clear();
        }
    }
    return holeVertices;
}
std::pair<std::vector<Eigen::MatrixXd>,std::vector<Eigen::VectorXi>> closestPointsOnMesh(std::vector<Eigen::MatrixXd> holeVerticesPos,Eigen::MatrixXd isomeshV,Eigen::MatrixXi isomeshF) {
    std::vector<Eigen::MatrixXd> closestPoints;
    std::vector<Eigen::VectorXi> closestFaces;
    igl::AABB<Eigen::MatrixXd,3> tree;
    tree.init(isomeshV,isomeshF);
    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;
    for (auto P: holeVerticesPos) {

        // igl::point_mesh_squared_distance(P,isomeshV,isomeshF,sqrD,I,C);
        tree.squared_distance(isomeshV,isomeshF,P,sqrD,I,C);
        closestPoints.push_back(C);
        closestFaces.push_back(I);
    }
    return std::make_pair(closestPoints,closestFaces);

}
Eigen::MatrixXd barycentricCoordinates(Eigen::MatrixXd isomeshV,Eigen::MatrixXi isomeshF,Eigen::MatrixXd points,Eigen::VectorXi closestFaces) {
    Eigen::MatrixXd L,Va,Vb,Vc; // barycentric coordinates
    Va = isomeshV(isomeshF(closestFaces,0),Eigen::all);
    Vb = isomeshV(isomeshF(closestFaces,1),Eigen::all);
    Vc = isomeshV(isomeshF(closestFaces,2),Eigen::all);
    // igl::barycentric_coordinates(closestPoints[0],faceVertices.col(0),faceVertices.col(1),faceVertices.col(2),L);
    // std::cout << "Barycentric coordinates : " << Va.row(0) <<Vb.row(0) <<Vc.row(0) << std::endl;
    igl::barycentric_coordinates(points,Va,Vb,Vc,L);
    return L;
}
// a list of your data for each frame
std::vector<std::vector<glm::vec3>> perFramePoints;

// UI state
Eigen::RowVector3d shift(0.6,0.03,0.00);
float x,y,z;
std::vector<Eigen::MatrixXd> closestPoints;
std::vector<Eigen::VectorXi> closestFaces;

// shift the to

int main(int argc, char **argv) {


    std::string filepath = "/Users/haoliang/Downloads/ongoing/signed-heat-demo-3d/data/zhldataset/bunny_hole.obj";
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();
    std::cout << "Vertices #: " << mesh->nVertices() << " Faces #:  "   << mesh->nFaces()<<std::endl;
    // some mesh info
    std::cout << "# of boundary vertices : "<<mesh->nExteriorHalfedges()<<std::endl;
    std::cout << "# of boundaries : "<<mesh->nBoundaryLoops()<<std::endl;
    // std::vector<std::vector<int>> holeVertices = holeVerticesFind(mesh);
    std::vector<Eigen::MatrixXd> holeVerticesPos = holeVerticesPosFind(mesh,geometry);

    // Read the complete mesh
    std::string isoMeshfile = "/Users/haoliang/Downloads/ongoing/result/isosurface.obj";
    Eigen::MatrixXd isomeshV;
    Eigen::MatrixXi isomeshF;
    igl::readOBJ(isoMeshfile, isomeshV, isomeshF);

    // // shift the to
    // Eigen::RowVector3d shift(0.2,0,5);
    // Initialize Polyscope
    polyscope::init();
    // // Load mesh
    polyscope::SurfaceMesh * holeMeshColor =  polyscope::registerSurfaceMesh("holeMesh", geometry->inputVertexPositions,mesh->getFaceVertexList());
    holeMeshColor->setSurfaceColor(glm::vec3(241.0/255,128.0/255,223.0/255));
    polyscope::SurfaceMesh * isoMeshColor = polyscope::registerSurfaceMesh("isoMesh", isomeshV.rowwise()+ shift, isomeshF);
    isoMeshColor->setSurfaceColor(glm::vec3(81.0/255,212.0/255,131.0/255));
    // closest point projection
    std::vector<Eigen::MatrixXd> closestPoints;
    std::vector<Eigen::VectorXi> closestFaces;
    std::tie(closestPoints,closestFaces) = closestPointsOnMesh(holeVerticesPos,isomeshV,isomeshF);
    //
    for (auto tmp : closestFaces) std::cout << tmp.size()<<" -- ";
    int select_idx = 4;

    /*
    for (int select_idx = 0; select_idx < closestFaces.size(); select_idx++) {
        // Add the curve network
        polyscope::CurveNetwork* curveHoleHandle = polyscope::registerCurveNetworkLine("holeLoop" + std::to_string(select_idx), holeVerticesPos[select_idx]);
        curveHoleHandle->setColor(glm::vec3(0,0,0));
        curveHoleHandle->setRadius(0.001);
        // polyscope::CurveNetwork* isoHoleHandle = polyscope::registerCurveNetworkLine("isoholeLoop"+ std::to_string(select_idx), closestPoints[select_idx].rowwise()+ shift);
        // isoHoleHandle->setRadius(0.001);
        // isoHoleHandle->setColor(glm::vec3(0,0,0));

        // add point cloud
        polyscope::PointCloud* holeVerticespoints = polyscope::registerPointCloud("holeVertices points"+ std::to_string(select_idx), holeVerticesPos[select_idx]);
        // holeVerticespoints->setPointColor(glm::vec3(249./255,0.0/255,0.0/255));
        holeVerticespoints->setPointColor(glm::vec3(0,0,0));
        polyscope::PointCloud* projVerticespoints = polyscope::registerPointCloud("projVertices points"+ std::to_string(select_idx), closestPoints[select_idx].rowwise()+ shift);
        // projVerticespoints->setPointColor(glm::vec3(0/255,249./255,0.0/255));
        projVerticespoints->setPointColor(glm::vec3(0,0,0));
        // set some options
        holeVerticespoints->setPointRadius(0.003);
        projVerticespoints->setPointRadius(0.003);
    }
    */
    // Add the curve network
    polyscope::CurveNetwork* curveHoleHandle = polyscope::registerCurveNetworkLine("holeLoop", holeVerticesPos[select_idx]);
    curveHoleHandle->setColor(glm::vec3(0,0,0));
    curveHoleHandle->setRadius(0.001);
    polyscope::CurveNetwork* isoHoleHandle = polyscope::registerCurveNetworkLine("isoholeLoop", closestPoints[select_idx].rowwise()+ shift);
    isoHoleHandle->setRadius(0.001);
    isoHoleHandle->setColor(glm::vec3(0,0,0));

    // add point cloud
    polyscope::PointCloud* holeVerticespoints = polyscope::registerPointCloud("holeVertices points", holeVerticesPos[select_idx]);
    // holeVerticespoints->setPointColor(glm::vec3(249./255,0.0/255,0.0/255));
    holeVerticespoints->setPointColor(glm::vec3(0,0,0));
    polyscope::PointCloud* projVerticespoints = polyscope::registerPointCloud("projVertices points", closestPoints[select_idx].rowwise()+ shift);
    // projVerticespoints->setPointColor(glm::vec3(0/255,249./255,0.0/255));
    projVerticespoints->setPointColor(glm::vec3(0,0,0));
    // set some options
    holeVerticespoints->setPointRadius(0.002);
    projVerticespoints->setPointRadius(0.002);

    // Connect map  poings
    std::vector<geometrycentral::Vector3> SegPoints;
    Eigen::MatrixXd closestPointsTmp = closestPoints[select_idx].rowwise()+ shift;
    for (int i = 0; i<closestFaces[select_idx].size(); i++) {
        SegPoints.push_back(geometrycentral::Vector3{holeVerticesPos[select_idx](i,0),holeVerticesPos[select_idx](i,1),holeVerticesPos[select_idx](i,2)});
        SegPoints.push_back(geometrycentral::Vector3{closestPointsTmp(i,0),closestPointsTmp(i,1),closestPointsTmp(i,2)});
    }
    polyscope::CurveNetwork* SegmentsLink = polyscope::registerCurveNetworkSegments("Segments", SegPoints);
    SegmentsLink->setRadius(0.001);

    Eigen::MatrixXd bc = barycentricCoordinates(isomeshV,isomeshF,closestPoints[select_idx],closestFaces[select_idx]);
    // std::cout << "Barycentric coordinates : " << Va.row(0) <<Vb.row(0) <<Vc.row(0) << std::endl;
    std::cout << "Barycentric coordinates example : " << bc.row(0) << std::endl;
    // std::cout << "Barycentric coordinates example : " << bc.row(0)[0]<<bc.row(0).y<<bc(0,2) << std::endl;
    //
    //
    // // search geodesic, Create the GeodesicAlgorithmExact
    // std::tie(isomesh_uptr, isogeometry_uptr) = makeManifoldSurfaceMeshAndGeometry(isomeshV, isomeshF);
    // isomesh = isomesh_uptr.release();
    // isogeometry = isogeometry_uptr.release();
    // GeodesicAlgorithmExact mmp(*isomesh, *isogeometry);
    //
    // // std::vector<geometrycentral::Vector3> geoPathPoints;
    // std::vector<geometrycentral::Vector3> geoPathPoints;
    // int interEdgeCount = 0;
    // int nPoj = closestFaces[select_idx].size();

    // for (int i = 0; i < nPoj; ++i) {
    //     Face sf = isomesh->face(closestFaces[select_idx][i]),tf = isomesh->face(closestFaces[select_idx][(i+1)%nPoj]);
    //     geometrycentral::Vector3 sfbc {bc(i%nPoj,0),bc(i%nPoj,1),bc(i%nPoj,2)},
    //     tfbc{bc((i+1)%nPoj,0),bc((i+1)%nPoj,1),bc((i+1)%nPoj,2)};
    //     mmp.propagate(SurfacePoint(sf, sfbc));
    //     std::vector<SurfacePoint> path = mmp.traceBack(SurfacePoint(tf, tfbc));
    //     geoPathPoints.push_back(geometrycentral::Vector3{closestPoints[select_idx](i,0),closestPoints[select_idx](i,1),closestPoints[select_idx](i,2)});
    //     std::vector<geometrycentral::Vector3> oneTrace;
    //     for (auto p : path) {
    //         // std::cout << p << std::endl;
    //         if (p.type == SurfacePointType::Edge) {
    //             interEdgeCount++;
    //             geometrycentral::Vector3 halfedgestartPos = isogeometry->inputVertexPositions[p.edge.firstVertex()];
    //             geometrycentral::Vector3 halfedgeendPos = isogeometry->inputVertexPositions[p.edge.secondVertex()];
    //             geometrycentral::Vector3 halfedgedir = halfedgeendPos - halfedgestartPos;
    //             geometrycentral::Vector3 interPos= halfedgestartPos + halfedgedir * p.tEdge ;
    //             oneTrace.push_back(interPos);
    //         }
    //         std::reverse(oneTrace.begin(), oneTrace.end());
    //         for (auto tmp : oneTrace) {
    //             geoPathPoints.push_back(tmp);
    //         }
    //     }
    //     std::cout << "search point at: " << i <<"  : total interEdgeCount # " <<interEdgeCount << std::endl;
    // }
    // std::cout << "The total # of intersection edge: "<< interEdgeCount << std::endl;
    // polyscope::CurveNetwork* isoGeoPath = polyscope::registerCurveNetworkLine("isoGeoPath", geoPathPoints);
    // isoGeoPath->setRadius(0.001);

    // polyscope::state::userCallback = myCallback; // specify the callback
    polyscope::show();
    }