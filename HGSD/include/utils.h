//
// Created by haoliang on 9/11/25.
//

#ifndef UTILS_H
#define UTILS_H

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
#include <unordered_set>
#include <unordered_map>
// Custom hash struct for std::pair
struct PairHash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        // A common way to combine hashes
        return h1 ^ (h2 << 1);
    }
};
struct remeshInfo {
    std::vector<std::vector<size_t>> addF;
    std::unordered_set<size_t> removedFaces;
    std::vector<geometrycentral::Vector3> addV;
};
struct  ConnectupdateInfo {
    geometrycentral::Vector3 s;
    geometrycentral::Vector3 t;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    std::vector<geometrycentral::Vector3> trace;
};

using namespace geometrycentral::surface;
struct  PathInfo {
    std::vector<SurfacePoint> inputPath;
    std::vector<size_t> vLink;
    size_t endPoint;
    size_t start;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;
};
std::vector<geometrycentral::Vector3> findGeopath(SurfacePoint source,SurfacePoint target,ManifoldSurfaceMesh* mesh,VertexPositionGeometry* geometry);
geometrycentral::Vector3 toXYZ(SurfacePoint p,ManifoldSurfaceMesh* mesh,VertexPositionGeometry* geometry);
// ConnectupdateInfo updateConnection(ConnectupdateInfo info);
PathInfo updateConnection(PathInfo info);
polyscope::SurfaceMesh* drawMesh(std::string name, Eigen::MatrixXd &V,Eigen::MatrixXi &F,
    double R=241,double G=128,double B=223) ;
std::pair<std::vector<std::vector<int>>,std::vector<Eigen::MatrixXd>> MeshHoleFind(ManifoldSurfaceMesh* mesh,VertexPositionGeometry* geometry);
std::pair<std::vector<Eigen::MatrixXd>,std::vector<Eigen::VectorXi>> ProjToMesh(std::vector<Eigen::MatrixXd> holeVerticesPos,
    Eigen::MatrixXd isomeshV,Eigen::MatrixXi isomeshF);
Eigen::MatrixXd ToBCCoordinates(Eigen::MatrixXd V,Eigen::MatrixXi F,Eigen::MatrixXd points,Eigen::VectorXi onFaces);
void bcProcess(Eigen::MatrixXd& bc, float alpha = 0.05,float beta = 0.95);
std::vector<SurfacePoint> robustPath(std::vector<SurfacePoint> path);
std::unordered_set<int> adjacentFaces(SurfacePoint p);
bool isVertex(SurfacePoint p);
int findCommonFace(std::unordered_set<int> a,std::unordered_set<int> b);
std::pair<Eigen::MatrixXd,Eigen::VectorXi> PointOnMesh(Eigen::MatrixXd Pts,Eigen::MatrixXd V,Eigen::MatrixXi F);

// remeshInfo reTriangulation(std::vector<SurfacePoint> tmppath);
std::pair<Eigen::MatrixXd,Eigen::MatrixXi> updateMesh(Eigen::MatrixXd V,Eigen::MatrixXi F,
    std::vector<std::vector<size_t>> addF,std::unordered_set<size_t> removedFaces,std::vector<geometrycentral::Vector3> addV);

class updateMesh {
public:
    Eigen::MatrixXd V,newV;
    Eigen::MatrixXi F,newF;
    // geometrycentral::Vector3 s{};
    // geometrycentral::Vector3 t{};
    updateMesh(Eigen::MatrixXd V,Eigen::MatrixXi F);
    bool updateConnection(geometrycentral::Vector3 s,geometrycentral::Vector3 t);


};
inline updateMesh::updateMesh(Eigen::MatrixXd V,Eigen::MatrixXi F) {
    this->V = V;
    this->F = F;
    this->newV = V;
    this->newF = F;
    // this->s = s;
    // this->t = t;
};

// bool updateMesh::updateConnection(geometrycentral::Vector3 s,geometrycentral::Vector3 t) {
//     if
//
//
// }






class MeshModify {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
public:
    // ManifoldSurfaceMesh* mesh;
    // VertexPositionGeometry* geometry;
    SurfaceMesh* mesh;
    VertexPositionGeometry* geometry;
    MeshModify(Eigen::MatrixXd V,Eigen::MatrixXi F);
    remeshInfo reTriangulation(std::vector<SurfacePoint> tmppath);
    std::vector<SurfacePoint> SeachPath2(SurfacePoint s,SurfacePoint t);
    std::vector<SurfacePoint> SeachPath(Face s,Eigen::MatrixXd sbc, Face t,Eigen::MatrixXd tbc);
    geometrycentral::Vector3 toXYZ(SurfacePoint p);

    // ~MeshModify();
};

inline MeshModify::MeshModify(Eigen::MatrixXd V,Eigen::MatrixXi F) {
    this->V = V;
    this->F = F;
    // std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr;
    // std::unique_ptr<VertexPositionGeometry> geometry_uptr;
    // std::tie(mesh_uptr, geometry_uptr) = makeManifoldSurfaceMeshAndGeometry(V, F);

    std::unique_ptr<SurfaceMesh> mesh_uptr;
    std::unique_ptr<VertexPositionGeometry> geometry_uptr;
    std::tie(mesh_uptr, geometry_uptr) = makeSurfaceMeshAndGeometry(V, F);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();
};
inline std::vector<SurfacePoint> MeshModify::SeachPath2(SurfacePoint s,SurfacePoint t) {
    GeodesicAlgorithmExact mmp(*mesh, *geometry);
    mmp.propagate(t);
    std::vector<SurfacePoint> path = mmp.traceBack(s);
    // the path is s -> t, if s is vertex delete duplicate
    if (s.faceCoords[0]==1 or s.faceCoords[1]==1 or s.faceCoords[2]==1) {path.erase(path.begin()+1);}

    return path;
};
inline std::vector<SurfacePoint> MeshModify::SeachPath(Face sf,Eigen::MatrixXd sbc, Face tf,Eigen::MatrixXd tbc) {
    GeodesicAlgorithmExact mmp(*mesh, *geometry);
    SurfacePoint s = SurfacePoint(sf,geometrycentral::Vector3{sbc(0),sbc(1),sbc(2)});
    SurfacePoint t = SurfacePoint(tf,geometrycentral::Vector3{tbc(0),tbc(1),tbc(2)});
    mmp.propagate(t);
    std::vector<SurfacePoint> path = mmp.traceBack(s);
    // the path is s -> t, if s is vertex delete duplicate
    if (s.faceCoords[0]==1 or s.faceCoords[1]==1 or s.faceCoords[2]==1) {path.erase(path.begin()+1);}
    // Eigen::VectorXi vec(4);
    // vec << 10, 0, 20, 0;
    // long nonZerosVec = (vec.array() != 0).count();
    return path;
};
inline geometrycentral::Vector3 MeshModify::toXYZ(SurfacePoint p) {
    geometrycentral::Vector3 xyz = geometrycentral::Vector3::zero();
    if (p.type == SurfacePointType::Edge) {
        geometrycentral::Vector3 halfedgestartPos = geometry->inputVertexPositions[p.edge.firstVertex()];
        geometrycentral::Vector3 halfedgeendPos = geometry->inputVertexPositions[p.edge.secondVertex()];
        geometrycentral::Vector3 halfedgedir = halfedgeendPos - halfedgestartPos;
        xyz = halfedgestartPos + halfedgedir * p.tEdge ;
    }
    else if (p.type == SurfacePointType::Face) {
        auto faces = mesh->getFaceVertexMatrix<double>();
        for (int i = 0; i < 3; i++) {
            xyz += geometry->inputVertexPositions[faces.row(p.face.getIndex())[i]] * p.faceCoords[i];
        }
    }
    else {
        xyz = geometry->inputVertexPositions[p.vertex.getIndex()];
    }
    return xyz;
}


#endif //UTILS_H
