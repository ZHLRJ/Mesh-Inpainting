//
// Created by haoliang on 10/8/25.
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
std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr,projmesh_uptr;
std::unique_ptr<VertexPositionGeometry> geometry_uptr,projgeometry_uptr;
// so we can more easily pass these to different classes while preserving syntax
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;
ManifoldSurfaceMesh* projmesh;
VertexPositionGeometry* projgeometry;
// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;
int main() {
    // Initialize Polyscope
    polyscope::init();
    std::string filepath = "/Users/haoliang/Downloads/ongoing/holeCutting.obj";
    std::string cutinfofile = "/Users/haoliang/Downloads/ongoing/loopinfo.txt";
    std::string filewithholepath = "/Users/haoliang/Downloads/ongoing/signed-heat-demo-3d/data/zhldataset/bunny_hole.obj";

    // Read the projection mesh
    Eigen::MatrixXd meshV,holemeshV;
    Eigen::MatrixXi meshF,holemeshF;
    std::vector<geometrycentral::Vector3> cutTrace;
    igl::readOBJ(filepath, meshV, meshF);
    igl::readOBJ(filewithholepath, holemeshV, holemeshF);
    std::tie(mesh_uptr, geometry_uptr) = makeManifoldSurfaceMeshAndGeometry(meshV, meshF);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();

    std::vector<bool> isProjPoints;
    std::vector<size_t> projPointsIdx;
    std::vector<size_t> HoleboundaryPointsIdx;
    std::unordered_map<size_t, size_t> patchVToholeVidx;
    std::ifstream inputFile(cutinfofile);

    int idx, val,holeidx;
    while (inputFile >> val >> idx>> holeidx) {
        isProjPoints.push_back(idx);
        projPointsIdx.push_back(val);
        HoleboundaryPointsIdx.push_back(holeidx);
        if (idx){patchVToholeVidx[val] = holeidx;}
        cutTrace.push_back(geometry->inputVertexPositions[mesh->vertex(val)]);
    }
    // add first for loop
    cutTrace.push_back(geometry->inputVertexPositions[mesh->vertex(projPointsIdx.front())]);
    inputFile.close(); // Close the file stream


    std::unordered_map<size_t, size_t> loopmap;
    std::unordered_set<size_t> loopset;
    std::unordered_set<size_t> insideFaces;
    std::unordered_set<size_t> insidevertices;
    std::deque<size_t> queue;
    std::unordered_set<size_t> visitedV;
    std::unordered_set<size_t> visitedF;

    size_t a,b,c;
    for (int i=0; i<projPointsIdx.size();i++) {

        loopmap[projPointsIdx[(i+1)%projPointsIdx.size()]] = projPointsIdx[i] ;
        // printf("DG : %lu -> %lu \n",projPointsIdx[i],projPointsIdx[(i+1)%projPointsIdx.size()]);
        insidevertices.insert(projPointsIdx[i]);
        loopset.insert(projPointsIdx[i]);
    }
    // for (auto p:loopmap){
    //     printf("parir (%lu --> %lu) \n",p.first,p.second);
    // }

    std::vector<std::array<double, 3>> fColor(meshF.rows());
    std::array<double, 3> common_row = {241./255,128./255,223./255};

    std::fill(fColor.begin(), fColor.end(), common_row);

    for (int iF=0; iF<meshF.rows();iF++) {
        //a -> b -> c
        a = meshF(iF,0);
        b = meshF(iF,1);
        c = meshF(iF,2);
        // fColor[iF] = std::array<double, 3>{241./255,128./255,223./255};
        // three points on hole
        if (loopset.find(a) != loopset.end() && loopset.find(b) != loopset.end() && loopset.find(c) != loopset.end()) {
            // keep orinitation
            if (loopmap[a]==b || loopmap[b]==c || loopmap[c]==a) {
                insideFaces.insert(iF);
                fColor[iF] = std::array<double, 3>{0.,1.,0.};
                // printf("face iF: %lu, (%lu,%lu,%lu) \n",iF,a,b,c);
            }
            // printf("Find boundary face: %lu (%lu,%lu,%lu) \n",iF,a,b,c);
        }
        else if (loopset.find(a) != loopset.end() && loopmap[a]==b) {
            insideFaces.insert(iF);
            fColor[iF] = std::array<double, 3>{0.,1.,0.};
            queue.push_back(c);
            // printf("v-a iF: %lu (%lu,%lu,%lu), inside v: %lu \n",iF,a,b,c,c);
            // printf("Find boundary face: %lu (%lu,%lu,%lu), edge: %lu -> %lu \n",iF,a,b,c,a,b);

        }
        else if (loopset.find(b) != loopset.end() && loopmap[b]==c) {
            insideFaces.insert(iF);
            fColor[iF] = std::array<double, 3>{0.,1.,0.};
            queue.push_back(a);
            // printf("v-b iF: %lu (%lu,%lu,%lu), inside v: %lu \n",iF,a,b,c,a);
            // printf("Find boundary face: %lu (%lu,%lu,%lu), edge: %lu -> %lu \n",iF,a,b,c,b,c);
        }
        else if (loopset.find(c) != loopset.end() && loopmap[c]==a) {
            insideFaces.insert(iF);
            fColor[iF] = std::array<double, 3>{0.,1.,0.};
            queue.push_back(b);
            // printf("v-c iF: %lu (%lu,%lu,%lu), inside v: %lu \n",iF,a,b,c,b);
            // printf("  Find boundary face: %lu (%lu,%lu,%lu), edge: %lu -> %lu  \n",iF,a,b,c,c,a);
        }
    }
    // for (auto v: queue){std::cout<<v<<std::endl;}
    while (!queue.empty()) {
        size_t vcur =queue.front();
        queue.pop_front();
        if (visitedV.find(vcur) == visitedV.end() && loopset.find(vcur)==loopset.end()) {
            visitedV.insert(vcur);
            for (auto f:mesh->vertex(vcur).adjacentFaces()) {
                if (visitedF.find(f.getIndex()) == visitedF.end()) {
                    visitedF.insert(f.getIndex());
                    size_t iF = f.getIndex();
                    insideFaces.insert(iF);
                    //a -> b -> c
                    a = meshF(iF,0);
                    b = meshF(iF,1);
                    c = meshF(iF,2);
                    fColor[iF] = std::array<double, 3>{0.,1.,0.};
                    if (visitedV.find(a)==visitedV.end() && loopset.find(a)==loopset.end()) {queue.push_back(a);}
                    if (visitedV.find(b)==visitedV.end() && loopset.find(b)==loopset.end()) {queue.push_back(b);}
                    if (visitedV.find(c)==visitedV.end() && loopset.find(c)==loopset.end()) {queue.push_back(c);}

                }
            }
        }
    }

    // cut off patch
    std::unordered_map<size_t, size_t> reorderMap;
    Eigen::MatrixXd patchV(loopset.size()+visitedV.size(),3);
    // printf("Total # of V: %lu  \n",loopset.size()+visitedV.size());
    Eigen::MatrixXi patchF(insideFaces.size(),3);
    size_t vCnt = 0;
    size_t FCnt = 0;
    for (auto iF:insideFaces) {
        //a -> b -> c
        for (int i=0;i<3;++i) {
            size_t idx = meshF(iF,i);
            if (reorderMap.find(idx) == reorderMap.end()) {
                reorderMap[idx] = vCnt++;
                printf("reordered map %lu --> %lu \n",idx,reorderMap[idx]);
                patchV.row(reorderMap[idx]) = meshV.row(idx);
            }
        }
        a = meshF(iF,0);
        b = meshF(iF,1);
        c = meshF(iF,2);
        patchF.row(FCnt) << reorderMap[a], reorderMap[b], reorderMap[c];
        // std::cout << patchF.row(FCnt) << std::endl;
        FCnt ++;

    }

    // stich path to hole
    // std::vector<bool> isProjPoints;
    // std::vector<size_t> projPointsIdx;
    // std::unordered_map<size_t, size_t> patchVToholeVidx;
    Eigen::MatrixXd copypatchV = patchV;
    std::deque<size_t> stack,twoends;
    Eigen::VectorXd first,second,vector,newPos;
    int firstIdx,secondIdx,patchBoundaryIdx;
    for (int j = 0; j <isProjPoints.size()+1; ++j) {
        stack.push_back(j%isProjPoints.size());
        if (stack.size()>=2 && isProjPoints[stack[0]] && isProjPoints[stack.back()]) {
            firstIdx = HoleboundaryPointsIdx[stack[0]];
            secondIdx = HoleboundaryPointsIdx[stack.back()];

            first = holemeshV.row(firstIdx);
            second = holemeshV.row(secondIdx);

            copypatchV.row(reorderMap[projPointsIdx[stack[0]]]) = first;
            copypatchV.row(reorderMap[projPointsIdx[stack.back()]]) = second;

            printf("first idx: %d  ",firstIdx); std::cout<<first.transpose()<<"\t";
            printf("second idx: %d  ",secondIdx); std::cout<<second.transpose()<<std::endl;
            stack.pop_front();

            vector = second - first;

            int nSeg = stack.size();

            double addCnt = 1;
            while (stack.size()>1 ) {
                size_t curIdx = stack.front();
                stack.pop_front();
                // printf("addCnt: %f << nSeg: %d << curIdx: %lu  \n",addCnt ,nSeg ,curIdx );
                newPos = first+ (addCnt/nSeg) * vector;
                addCnt ++;
                patchBoundaryIdx = projPointsIdx[curIdx];
                copypatchV.row(reorderMap[patchBoundaryIdx]) = newPos;
                printf("update  patchBoundaryIdx -> reorderMap[patchBoundaryIdx]: %d %lu  \n",patchBoundaryIdx,reorderMap[patchBoundaryIdx]);
                std::cout<<first.transpose()<<"\t";
            }
        }
    }


        glm::vec3 insidecolor = glm::vec3(1., 0., 0.);
    // // Load mesh
    polyscope::CurveNetwork *trace = polyscope::registerCurveNetworkLine("cutTrace ", cutTrace);
    trace->setColor(glm::vec3(247./255,5./255,5./255));
    trace->setRadius(0.0003);


    drawMesh("holeMesh", holemeshV, holemeshF);
    drawMesh("poissonMesh", meshV, meshF);
    drawMesh("patchMesh", patchV, patchF);
    drawMesh("copypatchV", copypatchV, patchF);
    polyscope::getSurfaceMesh("poissonMesh")->addFaceColorQuantity("fColor", fColor);

    polyscope::show();
}
