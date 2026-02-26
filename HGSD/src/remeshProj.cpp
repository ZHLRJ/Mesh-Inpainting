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
ManifoldSurfaceMesh* projmesh;
VertexPositionGeometry* projgeometry;
// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;
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
Eigen::RowVector3d shift(0.,0.,0.00);
float x,y,z;
std::vector<Eigen::MatrixXd> closestPoints;
std::vector<Eigen::VectorXi> closestFaces;

// shift the to
bool isValidProj(std::vector<Eigen::MatrixXd> HP,Eigen::MatrixXd projmeshV,
    Eigen::MatrixXi projmeshF) {

    // closest point projection
    std::vector<Eigen::MatrixXd> closestPoints;
    std::vector<Eigen::VectorXi> closestFaces;
    std::tie(closestPoints,closestFaces) = closestPointsOnMesh(HP,projmeshV,projmeshF);

    for (int holeidx = 0; holeidx < closestFaces.size(); ++holeidx) {
        Eigen::VectorXi projFaces = closestFaces[holeidx];
        std::unordered_map<size_t,std::vector<size_t>> infoDic;
        for (int projFacesidx = 0; projFacesidx < projFaces.size(); ++projFacesidx) {
            // infoDic[projFaces[projFacesidx]].emplace_back(projV.row(projFacesidx).data(),projV.row(projFacesidx).data() + projV.cols());
            infoDic[projFaces[projFacesidx]].emplace_back(projFacesidx);
        }
        for (auto item:infoDic) {
            // std::cout << item.first << " " << item.second.size() << std::endl;
            if (item.second.size() > 1) {
                std::cout << item.first << " " << item.second.size() << std::endl;
                // std::cout << "Not prefect poisson mesh" ;return false;
            }
        }
    }
    return true;
}
int main(int argc, char **argv) {


    std::string filepath = "/Users/haoliang/Downloads/ongoing/signed-heat-demo-3d/data/zhldataset/bunny_hole.obj";
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();

    std::vector<Eigen::MatrixXd> holeVerticesPos = holeVerticesPosFind(mesh,geometry);

    // Read the complete mesh
    std::string projMeshfile = "/Users/haoliang/Downloads/ongoing/result/poissonMesh.obj";
    Eigen::MatrixXd projmeshV;
    Eigen::MatrixXi projmeshF;
    igl::readOBJ(projMeshfile, projmeshV, projmeshF);
    std::cout << "Total # of vertices: "<<projmeshV.rows()<<std::endl;

    // Initialize Polyscope
    polyscope::init();
    // // Load mesh
    polyscope::SurfaceMesh * holeMeshColor =  polyscope::registerSurfaceMesh("holeMesh", geometry->inputVertexPositions,mesh->getFaceVertexList());
    holeMeshColor->setSurfaceColor(glm::vec3(241.0/255,128.0/255,223.0/255));
    polyscope::SurfaceMesh * projMeshColor = polyscope::registerSurfaceMesh("projMesh", projmeshV.rowwise()+ shift, projmeshF);
    projMeshColor->setSurfaceColor(glm::vec3(231.0/255,205.0/255,115.0/255));
    // closest point projection
    std::vector<Eigen::MatrixXd> closestPoints;
    std::vector<Eigen::VectorXi> closestFaces;
    std::tie(closestPoints,closestFaces) = closestPointsOnMesh(holeVerticesPos,projmeshV,projmeshF);
    //
    // for (auto projFaces:closestFaces) {
    std::vector<std::vector<size_t>> addF;
    std::vector<geometrycentral::Vector3> addV;
    std::unordered_set<size_t> removed;
    size_t cntAddVertex = 0;

    for (int holeidx = 0; holeidx < closestFaces.size(); ++holeidx) {
        Eigen::VectorXi projFaces = closestFaces[holeidx];
        Eigen::MatrixXd projV = closestPoints[holeidx];
        // std::unordered_map<size_t,std::vector<std::vector<double>>> infoDic;
        std::unordered_map<size_t,std::vector<size_t>> infoDic;
        for (int projFacesidx = 0; projFacesidx < projFaces.size(); ++projFacesidx) {
            // infoDic[projFaces[projFacesidx]].emplace_back(projV.row(projFacesidx).data(),projV.row(projFacesidx).data() + projV.cols());
            infoDic[projFaces[projFacesidx]].emplace_back(projFacesidx);
        }
        for (auto item:infoDic) {
            if (item.second.size() > 2){ std::cout << "Error: Three points in a triangle !" ;break; }
            if (item.second.size()>1) {
                // std::cout << item.first << " " << item.second.size() << std::endl;
                std::cout << item.first << " " << item.second[0] << " - " << item.second[1] << std::endl;
                // std::cout << projV.row(item.second[0])+ projV.row(item.second[1]) << std::endl;

                // add new point
                Eigen::Vector3d midpoint = (projV.row(item.second[0])+ projV.row(item.second[1]))/2 ;
                addV.push_back(geometrycentral::Vector3{midpoint.x(),midpoint.y(),midpoint.z()});
                removed.insert(item.first);
                for (int i = 0; i < 3; i++) {
                        size_t v1 = projmeshF(item.first,i%3);
                        size_t v2 = projmeshF(item.first,(i+1)%3);
                        addF.push_back(std::vector<size_t>{v1,v2,projmeshV.rows()+cntAddVertex});
                        // printf("add triangle: - %lu, %lu,%lu; ",v1,v2,projmeshV.rows()+cntAddVertex);
                      }
                cntAddVertex ++;
                    }

            }
        }
    // add mew vertices and faces
    Eigen::MatrixXi newprojmeshF(0,3)  ;
    Eigen::MatrixXd newprojmeshV = projmeshV;
    printf("add # of vertices %lu \n",addV.size());
    for (auto v :addV) {
        // std::cout << v << std::endl;
        newprojmeshV.conservativeResize(newprojmeshV.rows()+1, Eigen::NoChange);
        newprojmeshV.row(newprojmeshV.rows()-1)[0] = v.x;
        newprojmeshV.row(newprojmeshV.rows()-1)[1] = v.y;
        newprojmeshV.row(newprojmeshV.rows()-1)[2] = v.z;

    }
    printf("total # of vertices in newmesh %lu \n",newprojmeshV.rows());
    for (int i = 0;i<projmeshF.rows();i++) {
        if (removed.find(i) == removed.end()) {
            // std::cout << i << std::endl;
            newprojmeshF.conservativeResize(newprojmeshF.rows()+1, Eigen::NoChange);
            newprojmeshF.row(newprojmeshF.rows() - 1) = projmeshF.row(i);
        }
    }
    for (auto f :addF) {
        // std::cout<<f[0]<< f[1]<<f[2]<<std::endl;
        Eigen::Vector3i new_eigen_vector;
        new_eigen_vector << f[0], f[1], f[2];
        newprojmeshF.conservativeResize(newprojmeshF.rows()+1, Eigen::NoChange);
        newprojmeshF.row(newprojmeshF.rows()-1) = new_eigen_vector;
    }

    // std::cout << "is Valid poisson mesh? " << isValidProj(holeVerticesPos,projmeshV,projmeshF) << std::endl;
    // std::cout << "is Valid poisson mesh? " << isValidProj(holeVerticesPos,newprojmeshV,newprojmeshF) << std::endl;
    polyscope::SurfaceMesh * NewprojMesh = polyscope::registerSurfaceMesh("NewprojMesh", newprojmeshV, newprojmeshF);
    projMeshColor->setSurfaceColor(glm::vec3(231.0/255,205.0/255,115.0/255));
    polyscope::show();
    }