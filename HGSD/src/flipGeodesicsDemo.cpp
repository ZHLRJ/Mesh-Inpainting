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
// Eigen::RowVector3d shift(0.6,0.03,0.00);
Eigen::RowVector3d shift(0,0,0);
float x,y,z;
std::vector<Eigen::MatrixXd> closestPoints;
std::vector<Eigen::VectorXi> closestFaces;

// shift the to

int main(int argc, char **argv) {


    std::string filepath = "/Users/haoliang/Downloads/ongoing/signed-heat-demo-3d/data/zhldataset/bunny_hole.obj";
    std::string isoMeshfile = "/Users/haoliang/Downloads/ongoing/result/isosurface.obj";
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();

    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(isoMeshfile);
    isomesh = mesh_uptr.release();
    isogeometry = geometry_uptr.release();

    std::cout << "Vertices #: " << mesh->nVertices() << " Faces #:  "   << mesh->nFaces()<<std::endl;
    // some mesh info
    std::cout << "# of boundary vertices : "<<mesh->nExteriorHalfedges()<<std::endl;
    std::cout << "# of boundaries : "<<mesh->nBoundaryLoops()<<std::endl;
    std::vector<Eigen::MatrixXd> holeVerticesPos = holeVerticesPosFind(mesh,geometry);

    // Read the complete mesh

    Eigen::MatrixXd isomeshV;
    Eigen::MatrixXi isomeshF;
    igl::readOBJ(isoMeshfile, isomeshV, isomeshF);

    // // shift the to
    // Eigen::RowVector3d shift(0.2,0,5);
    // Initialize Polyscope
    polyscope::init();
    // // Load mesh
    // polyscope::SurfaceMesh * holeMeshColor =  polyscope::registerSurfaceMesh("holeMesh", geometry->inputVertexPositions,mesh->getFaceVertexList());
    // holeMeshColor->setSurfaceColor(glm::vec3(241.0/255,128.0/255,223.0/255));
    polyscope::SurfaceMesh * isoMeshColor = polyscope::registerSurfaceMesh("isoMesh", isomeshV.rowwise()+ shift, isomeshF);
    isoMeshColor->setSurfaceColor(glm::vec3(81.0/255,212.0/255,131.0/255));
    // closest point projection
    std::vector<Eigen::MatrixXd> closestPoints;
    std::vector<Eigen::VectorXi> closestFaces;
    std::tie(closestPoints,closestFaces) = closestPointsOnMesh(holeVerticesPos,isomeshV,isomeshF);
    //
    for (auto tmp : closestFaces) std::cout << tmp.size()<<" -- ";
    int select_idx = 4;
    // approx point to closest vertex
    Eigen::MatrixXd closestPointsOnMeshbc = barycentricCoordinates(isomeshV,isomeshF,closestPoints[select_idx],closestFaces[select_idx]);
    std::vector<geometrycentral::surface::Vertex> closestPointsOnMeshVertex;
    for (int i = 0; i < closestPointsOnMeshbc.rows(); ++i) {

        Eigen::MatrixXd::Index minColIndex;
        closestPointsOnMeshbc.row(i).minCoeff(&minColIndex);
        int vp = isomeshF(closestFaces[select_idx][i],minColIndex);
        printf("The approx vertex %d . \n",vp);
        // printf("The current temperature in %s is %.1f degrees Celsius.\n", city.c_str(), temperature);
        closestPointsOnMeshVertex.push_back(isomesh->vertex(vp));
    }

    // traverse all projected vertices
    std::vector<geometrycentral::Vector3> EdgePathPoints;
    std::vector<geometrycentral::Vector3> onesidePoints;

    Vertex vStart,vEnd;

    int n = (int)closestPointsOnMeshVertex.size();
    for (int i = 0; i < n; ++i) {
        vStart = closestPointsOnMeshVertex[i%n];
        vEnd = closestPointsOnMeshVertex[(i+1)%n];
        std::vector<Halfedge> dijkstraPath = shortestEdgePath(*isogeometry, vStart, vEnd);
        // std::cout <<" i "<< dijkstraPath.size() << ", "<< std::endl;
        // if (dijkstraPath.size())
        printf("The path %d size :  %.lu.\n", i, dijkstraPath.size());
        for (int j = 0; j < dijkstraPath.size(); ++j) {
            Halfedge h = dijkstraPath[j];
            EdgePathPoints.push_back(isogeometry->inputVertexPositions[h.tailVertex()]);
            onesidePoints.push_back(isogeometry->inputVertexPositions[h.prevOrbitVertex().tailVertex()]);
        }



    }

    polyscope::CurveNetwork* isoEdgePath = polyscope::registerCurveNetworkLine("isoGeoPath", EdgePathPoints);
    isoEdgePath->setRadius(0.001);

    polyscope::PointCloud* oneside = polyscope::registerPointCloud("onesidePoints", onesidePoints);
    // holeVerticespoints->setPointColor(glm::vec3(249./255,0.0/255,0.0/255));
    oneside->setPointColor(glm::vec3(0,0,0));

    polyscope::PointCloud* gtVertexPoints = polyscope::registerPointCloud("gtVertexPoints", closestPoints[select_idx]);
    // holeVerticespoints->setPointColor(glm::vec3(249./255,0.0/255,0.0/255));
    gtVertexPoints->setPointColor(glm::vec3(1,0,0));

    // polyscope::show();
    }