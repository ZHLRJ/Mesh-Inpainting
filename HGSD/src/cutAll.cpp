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
#include <iostream>
#include <format>
#include <unistd.h>
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr;
std::unique_ptr<VertexPositionGeometry> geometry_uptr;
// so we can more easily pass these to different classes while preserving syntax
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;

// 临时解决non manifold
std::unique_ptr<SurfaceMesh> tmpmesh_uptr;
std::unique_ptr<VertexPositionGeometry> tmpgeometry_uptr;
SurfaceMesh* tmpmesh;
VertexPositionGeometry* tmpgeometry;

int main() {
    // Initialize Polyscope
    polyscope::init();
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    // cut patch from given information
    std::string originalMeshFile = "/Users/haoliang/Downloads/ongoing/compareModel/result/2-1-Hole.obj";
    Eigen::MatrixXd patchV,originalMeshV; Eigen::MatrixXi patchF,originalMeshF;
    igl::readOBJ(originalMeshFile, originalMeshV, originalMeshF);
    drawMesh("original",originalMeshV, originalMeshF);
    // int fileIdx = 4;

    size_t newAddidx = originalMeshV.rows();
    Eigen::MatrixXd stitchedMeshV;
    Eigen::MatrixXi stitchedMeshF;
    // std::cout << "originalMeshV.size() = " << originalMeshV.size() << std::endl;
    for (int fileIdx = 0; fileIdx < 1; fileIdx++) {
        printf("Start loading patch %d \n", fileIdx);
        std::unordered_map<size_t, size_t> addPatchreMap;

        std::string inputMeshFile = "/Users/haoliang/Downloads/ongoing/HGSD/2-1Cut/holeCutting-2-1.obj";
        std::string inputMeshCutInfoFile = "/Users/haoliang/Downloads/ongoing/HGSD/2-1Cut/loopinfo-2-1.txt";
        Eigen::MatrixXd inputmeshV; Eigen::MatrixXi inputmeshF;
        igl::readOBJ(inputMeshFile, inputmeshV, inputmeshF);

        //  1.   store inputMeshCutInfoFile as std::vector < std::array<size_t > >
        //  {inputVidx, isProjected, holdMeshVidx}  ex: {9060	1	5908}
        std::ifstream inputFile(inputMeshCutInfoFile);
        std::vector<geometrycentral::Vector3> cutTrace;
        std::vector<std::array<size_t,3>> inputInfo;
        std::unordered_set<size_t> boundary,insideFaces,visitedV,visitedF;
        std::unordered_map<size_t, size_t> orientationcheck;
        std::deque<size_t> queue;
        std::unordered_map<size_t, std::set<int>> meshGraph;
        std::unordered_map<std::pair<size_t, size_t>, size_t, PairHash> edgeAdjFaces;
        size_t a,b,c;
        while (inputFile >> a >> b>> c) {
            inputInfo.push_back({a,b,c});
            geometrycentral::Vector3 point =geometrycentral::Vector3{inputmeshV(a,0),inputmeshV(a,1),inputmeshV(a,2)};
            cutTrace.push_back(point);
            boundary.insert(a);

        }
        cutTrace.push_back(cutTrace.front());
        inputFile.close(); // Close the file stream

        for (int i=0; i<inputInfo.size();i++) {
            orientationcheck[inputInfo[(i+1)%inputInfo.size()][0]] = inputInfo[i][0];
        }

        for (int iF=0; iF<inputmeshF.rows();iF++) {
            //a -> b -> c
            a = inputmeshF(iF,0);b = inputmeshF(iF,1);c = inputmeshF(iF,2);
            meshGraph[a].insert(iF);
            meshGraph[b].insert(iF);
            meshGraph[c].insert(iF);
            edgeAdjFaces[std::pair<int, int>{a,b}]= iF;
            edgeAdjFaces[std::pair<int, int>{b,c}]= iF;
            edgeAdjFaces[std::pair<int, int>{c,a}]= iF;
            if (boundary.count(a) && boundary.count(b) && boundary.count(c)) {
                // keep orinitation
                if (orientationcheck[a]==b || orientationcheck[b]==c || orientationcheck[c]==a) {
                    insideFaces.insert(iF);
                }
            }
            else if (boundary.count(a) && orientationcheck[a]==b) {
                insideFaces.insert(iF);
                queue.push_back(c);
            }
            else if (boundary.count(b) && orientationcheck[b]==c) {
                insideFaces.insert(iF);
                queue.push_back(a);
            }
            else if (boundary.count(c) && orientationcheck[c]==a) {
                insideFaces.insert(iF);
                queue.push_back(b);
            }
        }

        // BFS for inside face

        while (!queue.empty()) {
            size_t vcur =queue.front();
            queue.pop_front();
            // visitedV.insert(vcur);
            if (visitedV.count(vcur) == 0 && boundary.count(vcur)==0) {
                visitedV.insert(vcur);
                for (auto iF:meshGraph[vcur]) {
                    if (insideFaces.count(iF) == 0) {
                        insideFaces.insert(iF);
                        a = inputmeshF(iF,0);b = inputmeshF(iF,1);c = inputmeshF(iF,2);
                        if (visitedV.count(a)==0) {queue.push_back(a);}
                        if (visitedV.count(b)==0) {queue.push_back(b);}
                        if (visitedV.count(c)==0) {queue.push_back(c);}

                    }
                }
            }

        }
        // special case: three points on boundary around by insideFaces
        for (int iF=0; iF<inputmeshF.rows();iF++) {
            //a -> b -> c
            a = inputmeshF(iF,0);b = inputmeshF(iF,1);c = inputmeshF(iF,2);
            size_t f1,f2,f3;
            f1 = edgeAdjFaces[std::pair<size_t, size_t>{c,b}];
            f2 = edgeAdjFaces[std::pair<size_t, size_t>{b,a}];
            f3 = edgeAdjFaces[std::pair<size_t, size_t>{a,c}];
            if (insideFaces.count(iF)==0 && insideFaces.count(f1) && insideFaces.count(f2) && insideFaces.count(f3) ) {
                insideFaces.insert(iF);
            }

        }




        std::vector<std::array<double, 3>> fColor(inputmeshF.rows());
        std::array<double, 3> basecolor = {241./255,128./255,223./255};
        std::fill(fColor.begin(), fColor.end(), basecolor);
        for (int iF=0; iF<inputmeshF.rows();iF++) {
            if (insideFaces.count(iF)) {
                fColor[iF] = std::array<double, 3>{0.,1.,0.};
            }
        }
        // for (auto v:boundary) {
        //     std::cout<< visitedV.count(v)<<std::endl;
        // }

        // cut off patch

        std::unordered_map<size_t, size_t> reorderMap;
        Eigen::MatrixXd patchV(boundary.size()+visitedV.size(),3);
        // printf("Total # of V: %lu  \n",boundary.size()+visitedV.size());
        Eigen::MatrixXi patchF(insideFaces.size(),3);
        size_t vCnt = 0;
        size_t FCnt = 0;
        for (auto iF:insideFaces) {
            //a -> b -> c
            for (int i=0;i<3;++i) {
                size_t idx = inputmeshF(iF,i);
                if (reorderMap.find(idx) == reorderMap.end()) {
                    reorderMap[idx] = vCnt++;
                    // printf("reordered map %lu --> %lu \n",idx,reorderMap[idx]);
                    patchV.row(reorderMap[idx]) = inputmeshV.row(idx);
                }
            }
            a = inputmeshF(iF,0);
            b = inputmeshF(iF,1);
            c = inputmeshF(iF,2);
            patchF.row(FCnt) << reorderMap[a], reorderMap[b], reorderMap[c];
            // std::cout << patchF.row(FCnt) <<std::endl;
            FCnt ++;
        }
        // writeSurfaceMesh(*mesh, *geometry, "patch" + std::to_string(fileIdx) + ".obj");


        // prepare stitch


        Eigen::MatrixXd copypatchV = patchV;
        std::deque<size_t> stack,twoends;
        Eigen::VectorXd first,second,vector,newPos;
        int firstIdx,secondIdx,patchBoundaryIdx;
        for (int j = 0; j <inputInfo.size()+1; ++j) {
            // std::cout<<inputInfo[j][1]<<std::endl;
            if (j<inputInfo.size() && inputInfo[j][1]) {
                // printf("v map info check : %lu, %lu,  %lu \n",j,reorderMap[inputInfo[j][0] ],inputInfo[j][2]);
                // std::cout << reorderMap[inputInfo[j][0] ] << std::endl;
                addPatchreMap[reorderMap[inputInfo[j][0]] ] = inputInfo[j][2];
            }
            stack.push_back(j%inputInfo.size());
            if (stack.size()>=2 && inputInfo[stack[0]][1] && inputInfo[stack.back()][1]) {
                firstIdx = inputInfo[stack[0]][2];
                secondIdx = inputInfo[stack.back()][2];

                first = originalMeshV.row(firstIdx);
                second = originalMeshV.row(secondIdx);

                copypatchV.row(reorderMap[inputInfo[stack[0]][0] ]) = first;
                copypatchV.row(reorderMap[inputInfo[stack.back()][0]]) = second;
                //
                // printf("first idx: %d  ",firstIdx); std::cout<<first.transpose()<<"\t";
                // printf("second idx: %d  ",secondIdx); std::cout<<second.transpose()<<std::endl;
                stack.pop_front();

                vector = second - first;

                int nSeg = stack.size();

                double addCnt = 1;
                while (stack.size()>1 ) {
                    size_t curIdx = stack.front();
                    stack.pop_front();
                    // std::cout << addCnt << nSeg << curIdx << std::endl;
                    // printf("addCnt: %f << nSeg: %d << curIdx: %lu  \n",addCnt ,nSeg ,curIdx );
                    newPos = first+ (addCnt/nSeg) * vector;
                    addCnt ++;
                    patchBoundaryIdx = inputInfo[curIdx][0];
                    copypatchV.row(reorderMap[patchBoundaryIdx]) = newPos;
                    // printf("update  patchBoundaryIdx -> reorderMap[patchBoundaryIdx]: %d %lu  \n",patchBoundaryIdx,reorderMap[patchBoundaryIdx]);
                }
            }
        }
        stitchedMeshV.resize(originalMeshV.rows()+patchV.rows()-addPatchreMap.size(),3);
        stitchedMeshF.resize(originalMeshF.rows()+patchF.rows(),3);

        stitchedMeshV.block(0, 0, originalMeshV.rows(), originalMeshV.cols()) = originalMeshV;
        stitchedMeshF.block(0, 0, originalMeshF.rows(), originalMeshF.cols()) = originalMeshF;
        for (int i=0; i<patchF.rows(); ++i) {
            // std::cout<< newAddidx << std::endl;
            for (int j=0;j<3;++j) {
                size_t patchvidx = patchF(i,j);
                if (addPatchreMap.count(patchvidx) == 0) {
                    stitchedMeshV.row(newAddidx) = copypatchV.row(patchvidx);
                    addPatchreMap[patchvidx] = newAddidx++;
                    // printf("will add to stitch mesh %lu --> %lu \n",patchvidx,addPatchreMap[patchvidx]);
                }
                // printf("will add to stitch mesh %lu --> %lu  %lu \n",patchvidx,addPatchreMap[patchvidx],addPatchreMap.size());
            }
        }

        for (int i=0; i<patchF.rows(); ++i) {
            patchF(i,0);
            stitchedMeshF.row(originalMeshF.rows()+i) <<
                addPatchreMap[patchF(i,0)],addPatchreMap[patchF(i,1)],addPatchreMap[patchF(i,2)];
        }


        // printf("Total # of V: %lu  \n",boundary.size()+visitedV.size());
        originalMeshV = stitchedMeshV;
        originalMeshF = stitchedMeshF;



        // drawMesh("Patch" + std::to_string(fileIdx), copypatchV, patchF,156,245,190);

        // drawMesh("inputMesh", inputmeshV, inputmeshF,156,245,190);
        // polyscope::getSurfaceMesh("inputMesh")->addFaceColorQuantity("fColor", fColor);

        // drawMesh("stitchpatchV"+ std::to_string(fileIdx), copypatchV, patchF);
        // polyscope::CurveNetwork *trace = polyscope::registerCurveNetworkLine("cutTrace"+ std::to_string(fileIdx), cutTrace);
        // trace->setColor(glm::vec3(0./255,0./255,0./255));
        // trace->setRadius(0.001);
    }
    drawMesh("StitchedMesh", stitchedMeshV, stitchedMeshF,156,245,190);

    std::tie(tmpmesh_uptr, tmpgeometry_uptr) = makeSurfaceMeshAndGeometry(stitchedMeshV, stitchedMeshF);
    tmpmesh = tmpmesh_uptr.release();
    tmpgeometry = tmpgeometry_uptr.release();
    writeSurfaceMesh(*tmpmesh, *tmpgeometry, "StitchedMesh2-1.obj");

    polyscope::show();
}
