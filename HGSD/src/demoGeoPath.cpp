//
// Created by haoliang on 8/11/25.
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

int main(int argc, char **argv) {
    std::string filepath = "/Users/haoliang/Downloads/ongoing/crumpleddevelopable/CrumpledDevelopablepart.obj";
    Eigen::MatrixXd projmeshV;
    Eigen::MatrixXi projmeshF;
    igl::readOBJ(filepath, projmeshV, projmeshF);


    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();
    polyscope::init();
     // Load mesh
    GeodesicAlgorithmExact mmp(*mesh, *geometry);
    polyscope::SurfaceMesh* meshdemo = polyscope::registerSurfaceMesh("meshdemo", geometry->inputVertexPositions,mesh->getFaceVertexList());
    meshdemo->setSurfaceColor(glm::vec3(255.0/255,233.0/255,184.0/255));
    Face sf = mesh->face(1),tf = mesh->face(1);
    for (auto v:tf.adjacentVertices()) std::cout << v << std::endl;
    geometrycentral::Vector3 sfbc {1.,0.,0}, tfbc{0.0,1,0.};
    // mmp.propagate(SurfacePoint(sf, sfbc));
    mmp.propagate(SurfacePoint(mesh->vertex(1)));

    std::vector<SurfacePoint> path = mmp.traceBack(SurfacePoint(tf, tfbc));
    std::vector<geometrycentral::Vector3> onTrace;
    std::vector<geometrycentral::Vector3> strightLink;
    std::vector<geometrycentral::Vector3> addV;
    std::unordered_set<size_t> removed;
    // path.erase(path.begin()+1);

    for (auto p : path) {
        std::cout <<"path p: "<< p <<  std::endl;
        if (p.type == SurfacePointType::Vertex) {
            geometrycentral::Vector3 a = geometry->inputVertexPositions[p.vertex.getIndex()];
            std::cout << a << std::endl;
        }
        // if (p.type == SurfacePointType::Face) std::cout << p.faceCoords << std::endl;
        // if  (p.type == SurfacePointType::Vertex)
        //
        // if (p.type == SurfacePointType::Edge) {
        //
        //     // std::cout <<"edge cross: "<< p.edge.firstVertex() << "--> "<<p.edge.secondVertex() << std::endl;
        //     // std::cout << "Pass face: "<< p.edge.halfedge().face() <<" " << p.edge.halfedge().twin().face() << std::endl;
        //
        //     geometrycentral::Vector3 halfedgestartPos = geometry->inputVertexPositions[p.edge.firstVertex()];
        //     geometrycentral::Vector3 halfedgeendPos = geometry->inputVertexPositions[p.edge.secondVertex()];
        //     geometrycentral::Vector3 halfedgedir = halfedgeendPos - halfedgestartPos;
        //     geometrycentral::Vector3 interPos= halfedgestartPos + halfedgedir * p.tEdge ;
        //     onTrace.push_back(interPos);
        //     addV.push_back(interPos);
        //
        //
        //     // add faces
        //     removed.insert(p.edge.halfedge().face().getIndex());
        //     removed.insert(p.edge.halfedge().twin().face().getIndex());
        //     // for (auto f:p.edge.adjacentFaces()) removed.insert(f.getIndex());
        // }
        // else {
        //     geometrycentral::Vector3 point = geometrycentral::Vector3::zero();
        //     auto faces = mesh->getFaceVertexMatrix<double>();
        //     for (int i = 0; i < 3; i++) {
        //         point += geometry->inputVertexPositions[faces.row(p.face.getIndex())[i]] * p.faceCoords[i];
        //     }
        //     // addV.push_back(point);
        //     // std::cout << p.face.adjacentHalfedges() << std::endl;
        //     removed.insert(p.face.getIndex());
        //     onTrace.push_back(point);
        //     strightLink.push_back(point);
        //     addV.push_back(point);
        //     // std::cout << mesh->getFaceVertexMatrix<double>().row(1) << std::endl;
        //     // std::cout << p.face.getIndex() << std::endl;
        // }
        }

    // add projection point

    //
    // std::vector<std::vector<size_t>> addF;
    // auto faces = mesh->getFaceVertexMatrix<double>();
    // size_t cntAddVertex = 0;
    // printf("total # of vertices %lu , path length : %lu  \n",mesh->nVertices(),path.size());
    // for (size_t i = 0; i < path.size() -1 ; i++) {
    //     cntAddVertex +=1;
    //     if (i==0) {
    //         std::cout << "first surface:  ";
    //         SurfacePoint p1 = path[i];
    //
    //         std::unordered_set<size_t> pointOnedge{path[i+1].edge.firstVertex().getIndex(),path[i+1].edge.secondVertex().getIndex()};
    //         // for (auto e:p1.face.getIndex())
    //         for (auto p : p1.face.adjacentVertices()) {std::cout << p << " - ";}
    //         for (int k = 0; k < 3; k++) {
    //             size_t v1 = faces.row(p1.face.getIndex())[k%3];
    //             size_t v2 = faces.row(p1.face.getIndex())[(k+1)%3];
    //             if (pointOnedge.find(v1)!=pointOnedge.end() && pointOnedge.find(v2)!=pointOnedge.end()) {
    //                 addF.push_back(std::vector<size_t>{v1,cntAddVertex+mesh->nVertices(),mesh->nVertices()});
    //                 addF.push_back(std::vector<size_t>{mesh->nVertices(),cntAddVertex+mesh->nVertices(),v2});
    //                 printf("add triangle type 1-0 - %lu, %lu,%lu; ",v1,cntAddVertex+mesh->nVertices(),mesh->nVertices());
    //                 printf("add triangle type 1-0 - %lu, %lu,%lu; ",mesh->nVertices(),cntAddVertex+mesh->nVertices(),v2);
    //             }
    //             else {
    //                 addF.push_back(std::vector<size_t>{v1,v2,mesh->nVertices()});
    //                 printf("add triangle type 1-1 - %lu, %lu,%lu; ",v1,v2,mesh->nVertices());
    //             }
    //         }
    //         std::cout << std::endl;
    //     }
    //     else if (i+1 == path.size() -1) {
    //         std::unordered_set<size_t> pointOnedge{path[i].edge.firstVertex().getIndex(),path[i].edge.secondVertex().getIndex()};
    //         for (int k = 0; k < 3; k++) {
    //             size_t v1 = faces.row(path[i+1].face.getIndex())[k%3];
    //             size_t v2 = faces.row(path[i+1].face.getIndex())[(k+1)%3];
    //             if (pointOnedge.find(v1)!=pointOnedge.end() && pointOnedge.find(v2)!=pointOnedge.end()) {
    //                 addF.push_back(std::vector<size_t>{v1,cntAddVertex+mesh->nVertices()-1,cntAddVertex+mesh->nVertices()});
    //                 addF.push_back(std::vector<size_t>{v2,cntAddVertex+mesh->nVertices(),cntAddVertex+mesh->nVertices()-1});
    //                 printf("add triangle type 3-0 - %lu, %lu,%lu; ",v1,cntAddVertex+mesh->nVertices()-1,cntAddVertex+mesh->nVertices());
    //                 printf("add triangle type 3-0 - %lu, %lu,%lu; ",v2,cntAddVertex+mesh->nVertices(),cntAddVertex+mesh->nVertices()-1);
    //                 break;
    //             }
    //         }
    //     }
    //     else {
    //         // path[i].edge.adjacentFaces()
    //         std::unordered_set<size_t> faceset;
    //         std::unordered_set<size_t> curEdgeV{path[i].edge.firstVertex().getIndex(),path[i].edge.secondVertex().getIndex()};
    //         std::unordered_set<size_t> nextEdgeV{path[i+1].edge.firstVertex().getIndex(),path[i+1].edge.secondVertex().getIndex()};
    //
    //         size_t target_face;
    //         for (auto f:path[i].edge.adjacentFaces()) {faceset.insert(f.getIndex());}
    //         std::cout << "surface : ";
    //         for (auto f:path[i+1].edge.adjacentFaces()) {
    //             if (faceset.find(f.getIndex()) != faceset.end()) {
    //                 target_face = f.getIndex();
    //                 for (auto p : f.adjacentVertices()) {
    //                     std::cout << p << " - ";
    //                 }
    //                 break;
    //             }
    //         }
    //
    //
    //         std::vector<size_t> halfedge1,halfedge2;
    //         for (int k = 0; k < 3; k++) {
    //             size_t v1 = faces.row(target_face)[k%3];
    //             size_t v2 = faces.row(target_face)[(k+1)%3];
    //             if (nextEdgeV.find(v1)!=nextEdgeV.end() && nextEdgeV.find(v2)!=nextEdgeV.end()) {
    //                 halfedge2 = std::vector<size_t>{v1,v2};
    //             }
    //             else if (curEdgeV.find(v1)!=curEdgeV.end() && curEdgeV.find(v2)!=curEdgeV.end()) {
    //                 halfedge1 = std::vector<size_t>{v1,v2};
    //             }
    //         }
    //         printf("halfedge1: %lu -> ,%lu,halfedge2: %lu -> ,%lu",halfedge1[0],halfedge1[1],halfedge2[0],halfedge2[1]);
    //         if (halfedge1[1] == halfedge2[0]) {
    //             addF.push_back(std::vector<size_t>{halfedge1[0],cntAddVertex+mesh->nVertices()-1,halfedge2[1]});
    //             addF.push_back(std::vector<size_t>{cntAddVertex+mesh->nVertices()-1,halfedge1[1],cntAddVertex+mesh->nVertices()});
    //             addF.push_back(std::vector<size_t>{cntAddVertex+mesh->nVertices()-1,cntAddVertex+mesh->nVertices(),halfedge2[1]});
    //
    //             printf("add triangle type 2-0 - %lu, %lu,%lu; ",halfedge1[0],cntAddVertex+mesh->nVertices()-1,halfedge2[1]);
    //             printf("add triangle type 2-0 - %lu, %lu,%lu; ",cntAddVertex+mesh->nVertices()-1,halfedge1[1],cntAddVertex+mesh->nVertices());
    //             printf("add triangle type 2-0 - %lu, %lu,%lu; ",cntAddVertex+mesh->nVertices()-1,cntAddVertex+mesh->nVertices(),halfedge2[1]);
    //
    //
    //         }
    //         else {
    //             addF.push_back(std::vector<size_t>{cntAddVertex+mesh->nVertices()-1,halfedge1[1],halfedge2[0]});
    //             addF.push_back(std::vector<size_t>{cntAddVertex+mesh->nVertices()-1,cntAddVertex+mesh->nVertices(),halfedge2[1]});
    //             addF.push_back(std::vector<size_t>{cntAddVertex+mesh->nVertices()-1,halfedge2[0],cntAddVertex+mesh->nVertices()});
    //             printf(" add triangle type 2-1 - %lu, %lu,%lu; ",cntAddVertex+mesh->nVertices()-1,halfedge1[1],halfedge2[0]);
    //             printf(" add triangle type 2-1 - %lu, %lu,%lu; ",cntAddVertex+mesh->nVertices()-1,cntAddVertex+mesh->nVertices(),halfedge2[1]);
    //             printf(" add triangle type 2-1 - %lu, %lu,%lu; ",cntAddVertex+mesh->nVertices()-1,halfedge2[0],cntAddVertex+mesh->nVertices());
    //
    //         }
    //         std::cout << std::endl;
    //     }
    //
    //
    // }
    //
    //
    // for (auto f:removed) std::cout << f << " - ";
    // std::reverse(onTrace.begin(), onTrace.end());
    // polyscope::CurveNetwork* trace = polyscope::registerCurveNetworkLine("onTrace", onTrace);
    // trace->setColor(glm::vec3(247./255,5./255,5./255));
    // // polyscope::CurveNetwork* strighttrace = polyscope::registerCurveNetworkLine("strightLink", strightLink);
    // // trace->addNodeScalarQuantity("sample value", xC);
    // trace->setRadius(0.0015);
    // // strighttrace->setRadius(0.0001);
    // // strighttrace->setColor(glm::vec3(5./255,5./255,244./255));
    // polyscope::PointCloud* psCloud = polyscope::registerPointCloud("really great points", strightLink);
    //
    // // set some options
    // psCloud->setPointRadius(0.003);
    //
    // // add mew vertices and faces
    // Eigen::MatrixXi newprojmeshF(0,3)  ;
    // Eigen::MatrixXd newprojmeshV = projmeshV;
    // printf("add # of vertices %lu \n",addV.size());
    // for (auto v :addV) {
    //     // std::cout << v << std::endl;
    //     newprojmeshV.conservativeResize(newprojmeshV.rows()+1, Eigen::NoChange);
    //     newprojmeshV.row(newprojmeshV.rows()-1)[0] = v.x;
    //     newprojmeshV.row(newprojmeshV.rows()-1)[1] = v.y;
    //     newprojmeshV.row(newprojmeshV.rows()-1)[2] = v.z;
    //
    // }
    // printf("total # of vertices in newmesh %lu \n",newprojmeshV.rows());
    // for (int i = 0;i<projmeshF.rows();i++) {
    //     if (removed.find(i) == removed.end()) {
    //         // std::cout << i << std::endl;
    //         newprojmeshF.conservativeResize(newprojmeshF.rows()+1, Eigen::NoChange);
    //         newprojmeshF.row(newprojmeshF.rows() - 1) = projmeshF.row(i);
    //     }
    // }
    // for (auto f :addF) {
    //     // std::cout<<f[0]<< f[1]<<f[2]<<std::endl;
    //     Eigen::Vector3i new_eigen_vector;
    //     new_eigen_vector << f[0], f[1], f[2];
    //     newprojmeshF.conservativeResize(newprojmeshF.rows()+1, Eigen::NoChange);
    //     newprojmeshF.row(newprojmeshF.rows()-1) = new_eigen_vector;
    // }
    //
    // polyscope::SurfaceMesh * projMeshColor = polyscope::registerSurfaceMesh("projMesh", newprojmeshV, newprojmeshF);
    // projMeshColor->setSurfaceColor(glm::vec3(231.0/255,205.0/255,115.0/255));

    // polyscope::show();


}
