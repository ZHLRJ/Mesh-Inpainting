//
// Created by haoliang on 9/11/25.
//
#include "utils.h"

// std::pair<std::vector<std::array<size_t,3>>,std::vector<geometrycentral::Vector3>> ReadinputMeshCutInfoFile(std::string inputMeshCutInfoFile) {
//     std::ifstream inputFile(inputMeshCutInfoFile);
//     std::vector<geometrycentral::Vector3> cutTrace;
//     std::vector<std::array<size_t,3>> inputInfo;
//     size_t a,b,c;
//     while (inputFile >> a >> b>> c) {
//         inputInfo.push_back({a,b,c});
//         geometrycentral::Vector3 point =geometrycentral::Vector3{inputmeshV(a,0),inputmeshV(a,1),inputmeshV(a,2)};
//         cutTrace.push_back(point);
//     }
//     // add first for loop
//     cutTrace.push_back(cutTrace.front());
//     inputFile.close(); // Close the file stream
// }


polyscope::SurfaceMesh* drawMesh(std::string name, Eigen::MatrixXd &V,Eigen::MatrixXi &F,
    double R,double G,double B) {
    polyscope::SurfaceMesh * meshHandle =  polyscope::registerSurfaceMesh(name, V,F);
    meshHandle->setSurfaceColor(glm::vec3(R/255.,G/255.,B/255.));
    meshHandle->setEdgeWidth(0.8);
    // return  meshHandle;
}

std::pair<std::vector<std::vector<int>>,std::vector<Eigen::MatrixXd>> MeshHoleFind(ManifoldSurfaceMesh* mesh,VertexPositionGeometry* geometry) {
    std::unordered_map<int, int> halfEdgeMap;
    for (Halfedge e :mesh->exteriorHalfedges()) {
        halfEdgeMap[e.sibling().vertex().getIndex()] = e.vertex().getIndex();
    }
    std::unordered_set<int> isVisited;
    std::vector<std::vector<int>> holesVindex;
    std::vector<int> curLoop;
    std::vector<Eigen::MatrixXd> holeVertices;
    for (std::pair<int,int> v: halfEdgeMap) {
        if (isVisited.find(v.first) == isVisited.end()) {
            curLoop.push_back(v.second);
            while (isVisited.find(halfEdgeMap[curLoop.back()]) == isVisited.end()){

                // isVisited.insert(halfEdgeMap[curLoop.back()]);
                isVisited.insert(curLoop.back());
                // std::cout << curLoop.back() << "-- " << halfEdgeMap[curLoop.back()] << std::endl;
                if (isVisited.find(halfEdgeMap[curLoop.back()]) == isVisited.end())
                    curLoop.push_back(halfEdgeMap[curLoop.back()]);


            }
            // std::cout << curLoop.back() << "-- ";
            Eigen::MatrixXd curLoopPosMatrix(curLoop.size(),3);
            for (int i = 0; i < curLoop.size(); i++) {
                curLoopPosMatrix.row(i) <<geometry->inputVertexPositions[curLoop[i]].x, geometry->inputVertexPositions[curLoop[i]].y,geometry->inputVertexPositions[curLoop[i]].z;
            }
            holeVertices.push_back(curLoopPosMatrix);
            holesVindex.push_back(curLoop);
            // printf("curLoop : %d, curLoopPosMatrix %d \n",curLoop.size(),curLoopPosMatrix.rows());
            curLoop.clear();
            // std::cout << curLoopPosMatrix.row(0) << " -- " << curLoopPosMatrix.row(curLoopPosMatrix.rows()-1)<< std::endl;
        }

    }
    return std::make_pair(holesVindex,holeVertices);
    // return holeVertices;
}
std::pair<std::vector<Eigen::MatrixXd>,std::vector<Eigen::VectorXi>> ProjToMesh(std::vector<Eigen::MatrixXd> holeVerticesPos,Eigen::MatrixXd isomeshV,Eigen::MatrixXi isomeshF) {
    std::vector<Eigen::MatrixXd> closestPoints;
    std::vector<Eigen::VectorXi> closestFaces;
    igl::AABB<Eigen::MatrixXd,3> tree;
    tree.init(isomeshV,isomeshF);
    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;
    for (auto P: holeVerticesPos) {
        tree.squared_distance(isomeshV,isomeshF,P,sqrD,I,C);
        closestPoints.push_back(C);
        closestFaces.push_back(I);
    }
    return std::make_pair(closestPoints,closestFaces);
}

std::pair<Eigen::MatrixXd,Eigen::VectorXi> PointOnMesh(Eigen::MatrixXd P,Eigen::MatrixXd V,Eigen::MatrixXi F) {
    igl::AABB<Eigen::MatrixXd,3> tree;
    tree.init(V,F);
    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;
    tree.squared_distance(V,F,P,sqrD,I,C);

    return std::make_pair(C,I);
}
Eigen::MatrixXd ToBCCoordinates(Eigen::MatrixXd V,Eigen::MatrixXi F,Eigen::MatrixXd points,Eigen::VectorXi onFaces) {
    Eigen::MatrixXd L,Va,Vb,Vc; // barycentric coordinates
    Va = V(F(onFaces,0),Eigen::all);
    Vb = V(F(onFaces,1),Eigen::all);
    Vc = V(F(onFaces,2),Eigen::all);
    // igl::barycentric_coordinates(closestPoints[0],faceVertices.col(0),faceVertices.col(1),faceVertices.col(2),L);
    // std::cout << "Barycentric coordinates : " << Va.row(0) <<Vb.row(0) <<Vc.row(0) << std::endl;
    igl::barycentric_coordinates(points,Va,Vb,Vc,L);
    return L;
}

void bcProcess(Eigen::MatrixXd& bc, float alpha,float beta) {
    for (auto i=0; i<bc.rows(); ++i) {
        Eigen::Index maxidx,minidx;
        bc.row(i).maxCoeff(&maxidx);
        bc.row(i).minCoeff(&minidx);
        if  (bc.row(i).minCoeff() <=alpha) {
            bc(i,maxidx) += bc.row(i).minCoeff();
            bc(i,minidx) = 0.;
        }
        for (int j = 0; j < 3; ++j) {
            if (bc(i,j)>=beta) {
                bc(i,j) = 1.; bc(i,(j+1)%3) = 0.;bc(i,(j+2)%3) = 0.;
                break;
            }
        }
    }
}
std::vector<SurfacePoint> robustPath(std::vector<SurfacePoint> path) {
    std::vector<SurfacePoint> newPath,final;
    for (int pi=0;pi<path.size();pi++) {
        // approx some edge point to vertex
        SurfacePoint p = path[pi];
        if (path[pi].type == SurfacePointType::Edge) {
            if (path[pi].tEdge <=0.07){
                p = path[pi].edge.firstVertex();
            }
            else if (path[pi].tEdge >=0.93) {
                p = path[pi].edge.secondVertex();
            }
        }
        if (path[pi].type == SurfacePointType::Face) {
            std::unordered_map<int,double> vTobc;
            int index = 0;
            for (auto v:path[pi].face.adjacentVertices()) {
                if (path[pi].faceCoords[index]==1) {
                    p = SurfacePoint(v);
                    break;
                }
                if (path[pi].faceCoords[index]!=0) {
                    vTobc[v.getIndex()] = path[pi].faceCoords[index] ;
                }
                index +=1;
            }
            if (vTobc.size()==2) {
                for (auto e:path[pi].face.adjacentEdges()) {
                    if (vTobc.find(e.firstVertex().getIndex()) != vTobc.end() && vTobc.find(e.secondVertex().getIndex()) != vTobc.end()) {
                        p = SurfacePoint(e,vTobc[e.secondVertex().getIndex()]);
                    }
                }
            }
        }
        if (newPath.empty() || p!=newPath.back()) {
            newPath.push_back(p);
        }
    }

    for (auto p : newPath) {
        if (final.empty()) {
            final.push_back(p);
        }
        else if(p.type!=final.back().type){
            final.push_back(p);
        }
        else {
            if (p.type==SurfacePointType::Edge && abs(p.tEdge-final.back().tEdge)>1e-8 || p.edge!=final.back().edge) {
                // printf( "two tedge %f , -- %f",p.tEdge, final.back().tEdge);
                final.push_back(p);
            }
        }
    }
    // for (auto p : final) {std::cout << p << std::endl;}
    return final;
}

std::vector<SurfacePoint> robustPath2(std::vector<SurfacePoint> path,ManifoldSurfaceMesh* mesh,VertexPositionGeometry* geometry) {
    std::vector<SurfacePoint> newPath,final;
    for (int pi=0;pi<path.size();pi++) {
        // approx some edge point to vertex
        SurfacePoint p = path[pi];
        if (path[pi].type == SurfacePointType::Edge) {
            if (path[pi].tEdge <=0.05){
                p = path[pi].edge.firstVertex();
            }
            else if (path[pi].tEdge >=0.95) {
                p = path[pi].edge.secondVertex();
            }
        }

        if (path[pi].type == SurfacePointType::Face) {
            std::unordered_map<int,double> vTobc;
            int index = 0;
            for (auto v:path[pi].face.adjacentVertices()) {
                if (path[pi].faceCoords[index]==1) {
                    p = SurfacePoint(v);
                    break;
                }
                if (path[pi].faceCoords[index]!=0) {
                    vTobc[v.getIndex()] = path[pi].faceCoords[index] ;
                }
                index +=1;
            }
            if (vTobc.size()==2) {
                for (auto e:path[pi].face.adjacentEdges()) {
                    if (vTobc.find(e.firstVertex().getIndex()) != vTobc.end() && vTobc.find(e.secondVertex().getIndex()) != vTobc.end()) {
                        p = SurfacePoint(e,vTobc[e.secondVertex().getIndex()]);
                    }
                }
            }
        }
        newPath.push_back(p);
        // if (newPath.empty() || p!=newPath.back()) {
        //     newPath.push_back(p);
        // }
    }

    for (auto p : newPath) {
        if (final.empty()) {
            final.push_back(p);
        }
        else if(  (toXYZ(p,mesh,geometry) - toXYZ(final.back(),mesh,geometry)).norm2()>1e-8    ) {
            final.push_back(p);
        }
        // else if(p.type!=final.back().type){
        //     final.push_back(p);
        // }
        // else {
        //     if (p.type==SurfacePointType::Edge && abs(p.tEdge-final.back().tEdge)>1e-8 || p.edge!=final.back().edge) {
        //         // printf( "two tedge %f , -- %f",p.tEdge, final.back().tEdge);
        //         final.push_back(p);
        //     }
        // }
    }
    // for (auto p : final) {std::cout << p << std::endl;}
    return final;
}
std::vector<SurfacePoint> robustPath3(std::vector<SurfacePoint> path,ManifoldSurfaceMesh* mesh,VertexPositionGeometry* geometry) {
    // all cross should be edge point
    std::vector<SurfacePoint> approxPath,finalPath;
    for (auto p:path) {
        // approx some edge point to vertex
        if (p.type == SurfacePointType::Edge) {
            if (p.tEdge <=0.1){
                p = p.edge.firstVertex();
            }
            else if (p.tEdge >=0.9) {
                p = p.edge.secondVertex();
            }
        }
        // bc process will eliminate all bc(1,0,0) type, but bc(1,0,0) may show in test
        if (p.type == SurfacePointType::Face) {
            int index = 0;
            for (auto v:p.face.adjacentVertices()) {
                if (p.faceCoords[index]==1) {
                    // std::cout<<"original p"<<p<<std::endl;
                    p = SurfacePoint(v);
                    // printf("update p ");
                    // std::cout<<p<<std::endl;
                }
                index ++;
            }
        }
        approxPath.push_back(p);
    }
    // remove duplicate points
    for (auto p : approxPath) {
        if (finalPath.empty()) {
            finalPath.push_back(p);
        }
        else if(  (toXYZ(p,mesh,geometry) - toXYZ(finalPath.back(),mesh,geometry)).norm2()>1e-8    ) {
            finalPath.push_back(p);
        }
    }
    return finalPath;
}

std::vector<SurfacePoint> robustPath4(std::vector<SurfacePoint> path,ManifoldSurfaceMesh* mesh,VertexPositionGeometry* geometry) {
    // first check if first and last points is represent by bc(a,b,0) , if then transfer to edge format for consistence
    for (int pi = 0; pi < path.size();pi++) {
        SurfacePoint p = path[pi];
        if (path[pi].type == SurfacePointType::Face ) {
            std::unordered_map<int,double> vTobc;
            int index = 0;
            for (auto v:path[pi].face.adjacentVertices()) {
                if (path[pi].faceCoords[index]==1) {
                    p = SurfacePoint(v);
                    break;
                }
                if (path[pi].faceCoords[index]!=0) {
                    vTobc[v.getIndex()] = path[pi].faceCoords[index] ;
                }
                index +=1;
            }
            if (vTobc.size()==2) {
                for (auto e:path[pi].face.adjacentEdges()) {
                    if (vTobc.find(e.firstVertex().getIndex()) != vTobc.end() && vTobc.find(e.secondVertex().getIndex()) != vTobc.end()) {
                        p = SurfacePoint(e,vTobc[e.secondVertex().getIndex()]);
                    }
                }
            }
        }
        if (p.type == SurfacePointType::Edge) {
            if (p.tEdge <=0.08){
                p = p.edge.firstVertex();
            }
            else if (p.tEdge >=0.92) {
                p = p.edge.secondVertex();
            }
        }
        path[pi] = p;

    }
    std::vector<SurfacePoint> newPath,Final;


    for (auto p : path) {

        if (newPath.empty()) {
            newPath.push_back(p);
        }

        else if(  (toXYZ(p,mesh,geometry) - toXYZ(newPath.back(),mesh,geometry)).norm2()>1e-12    ) {
            newPath.push_back(p);

        }
    }
    std::unordered_set<int> badPoints;
    SurfacePoint p,pre,nxt;

    for (int i=0;i<newPath.size();i++) {
        p = newPath[i];
        // e --> v, del e
        if (p.type == SurfacePointType::Vertex) {
            size_t vidx = p.vertex.getIndex();
            if (i-1!=0 && newPath[i-1].type == SurfacePointType::Edge) {
                pre = newPath[i-1];
                size_t ef  = pre.edge.firstVertex().getIndex();
                size_t es = pre.edge.secondVertex().getIndex();
                if (vidx == ef || vidx == es) {
                    badPoints.insert(i-1);
                }
            }
            if (i+1!=newPath.size()-1 && newPath[i+1].type == SurfacePointType::Edge) {
                nxt = newPath[i+1];
                size_t ef  = nxt.edge.firstVertex().getIndex();
                size_t es = nxt.edge.secondVertex().getIndex();
                if (vidx == ef || vidx == es) {
                    badPoints.insert(i+1);
                }
            }
        }
    }

    for (int i=0;i<newPath.size();i++) {
        if (badPoints.find(i) == badPoints.end()) {
            Final.push_back(newPath[i]);
        }
    }
    return Final;
}
bool isVertex(SurfacePoint p) {
    if (p.type == SurfacePointType::Vertex){return true;}
    if (p.type == SurfacePointType::Face && p.faceCoords[0]==1 or p.faceCoords[1]==1 or p.faceCoords[2]==1 ){return true;}
    return false ;
}

std::unordered_set<int> findFace(SurfacePoint p) {
    std::unordered_set<int> fpool ;
    if (p.type == SurfacePointType::Edge) {
        for (auto f:p.edge.adjacentFaces()){fpool.insert(f.getIndex());}
    }
    if (p.type == SurfacePointType::Face) {fpool.insert(p.face.getIndex());}

    if (p.type == SurfacePointType::Vertex) {for (auto f:p.vertex.adjacentFaces()){fpool.insert(f.getIndex());}}
    // for (auto f:fpool){std::cout<<f<<"-- ";}
    // std::cout<<std::endl;
    return fpool;
}
int findCommonFace(std::unordered_set<int> a,std::unordered_set<int> b) {
    int cnt = 0;
    int ans = -1;
    for (auto i:a) {
        if (b.find(i) != b.end()){ cnt++;ans = i; }
    }
    if (cnt ==1) return ans;
    return -1;

}
int findCommonFaceV2(SurfacePoint v1,SurfacePoint v2) {
    std::unordered_set<int> a = findFace(v1);
    std::unordered_set<int> b = findFace(v2);
    int cnt = 0;
    int ans = -1;
    for (auto i:a) {
        if (b.find(i) != b.end()){ cnt++;ans = i; }
    }
    if (cnt ==1) return ans;
    for (auto i:a) {std::cout<<i<<" ";}
    std::cout<<std::endl;
    for (auto i:b) {std::cout<<i<<" ";}
    std::cout<<std::endl;
    printf("No common face exist! \n");
    return -1;
}


std::pair<Eigen::MatrixXd,Eigen::MatrixXi> updateMesh(Eigen::MatrixXd V,Eigen::MatrixXi F,
    std::vector<std::vector<size_t>> addF,std::unordered_set<size_t> removedFaces,std::vector<geometrycentral::Vector3> addV) {
    Eigen::MatrixXd newV = V;
    Eigen::MatrixXi newF(0,3);

    // add mew vertices and faces
    for (auto v :addV) {
        // std::cout << v << std::endl;
        newV.conservativeResize(newV.rows()+1, Eigen::NoChange);
        newV.row(newV.rows()-1)[0] = v.x;
        newV.row(newV.rows()-1)[1] = v.y;
        newV.row(newV.rows()-1)[2] = v.z;

    }
    // printf("total # of vertices in newmesh %lu \n",newprojmeshV.rows());
    for (int i = 0;i<F.rows();i++) {
        if (removedFaces.find(i) == removedFaces.end()) {
            // std::cout << i << std::endl;
            newF.conservativeResize(newF.rows()+1, Eigen::NoChange);
            newF.row(newF.rows() - 1) = F.row(i);
        }
    }
    for (auto f :addF) {
        Eigen::Vector3i new_eigen_vector;
        new_eigen_vector << f[0], f[1], f[2];
        newF.conservativeResize(newF.rows()+1, Eigen::NoChange);
        newF.row(newF.rows()-1) = new_eigen_vector;
    }

    return std::make_pair(newV,newF);

}
std::vector<geometrycentral::Vector3> findGeopath(SurfacePoint source,SurfacePoint target,ManifoldSurfaceMesh* mesh,VertexPositionGeometry* geometry) {
    std::vector<geometrycentral::Vector3> onTrace;
    GeodesicAlgorithmExact mmp(*mesh, *geometry);
    mmp.propagate(target);
    std::vector<SurfacePoint> path = mmp.traceBack(source);

    for (auto p:path) {
        onTrace.push_back(toXYZ(p,mesh,geometry));
    }
return onTrace;
}
remeshInfo MeshModify::reTriangulation(std::vector<SurfacePoint> tmppath) {
    remeshInfo info;
    // std::vector<std::vector<size_t>> addF;
    std::vector<SurfacePoint> path = robustPath(tmppath);
    // size_t cntAddVertex = mesh->nVertices();
    size_t nV = mesh->nVertices();
    // for (int pi=0;pi<path.size();pi++) {
    //     if (!isVertex(path[pi])){
    //         if (info.addV.empty() || (info.addV.back()-toXYZ(path[pi])).norm2()>1e-8  )
    //             info.addV.push_back(toXYZ(path[pi]));
    //             // std::cout << toXYZ(path[pi]) << std::endl;
    //         // std::cout << (info.addV.back()-toXYZ(path[pi])).norm2() << std::endl;
    //             }
    // }
    size_t increment = 0;
    // std::unordered_set<size_t> removedFaces;
    for (int pi=0;pi<path.size()-1;pi++) {
        SurfacePoint v1 = path[pi];
        SurfacePoint v2 = path[(pi+1)%path.size()];
        // std::cout<<info.addV.size()<<std::endl;
        std::unordered_set<int> v1f = findFace(v1);
        std::unordered_set<int> v2f = findFace(v2);
        printf("faces include in v1: ");
        for (auto f:v1f) {std::cout<<f<<" , ";}std::cout<<std::endl;
        printf("faces include in v2: ");
        for (auto f:v2f) {std::cout<<f<<" , ";}std::cout<<std::endl;
        int fcommonidx = findCommonFace(v1f,v2f);

        std::cout << v1  << " -- " << v2 << std::endl;
        std::cout << "find common face in: "<< fcommonidx << std::endl;
        // printf("have # vertices:   %lu,  %lu, pi = %d \n",nV,info.addV.size(),pi );
        // break;
        // 1. v—e  or e-v
        if ((isVertex(v1) && v2.type == SurfacePointType::Edge) | (isVertex(v2) && v1.type == SurfacePointType::Edge)) {
            bool swaped = false;
            if (isVertex(v2)) {
                std::swap(v1,v2);
                swaped = true;
            }
            // three point not in same edge
            if (fcommonidx!=-1) {
                for (int vidx = 0;vidx<3;vidx++) {
                    if (F(fcommonidx,vidx) ==v1.vertex.getIndex()) {
                        size_t a = F(fcommonidx,vidx);
                        size_t b = F(fcommonidx,(vidx+1)%3);
                        size_t c = F(fcommonidx,(vidx+2)%3);
                        info.removedFaces.insert(fcommonidx);
                        info.addF.push_back(std::vector<size_t>{a,b,nV+increment});
                        info.addF.push_back(std::vector<size_t>{nV+increment,c,a});
                        printf("add triangle type v-e/e-v - %lu, %lu,%lu; and %lu, %lu,%lu. \n\n",a,b,nV+increment,nV+increment,c,a);
                        if (swaped){increment++;}

                        break;
                    }
                }
            }
        }
        // 2. v-f or f-v

        if ((isVertex(v1) && v2.type == SurfacePointType::Face) | (isVertex(v2) && v1.type == SurfacePointType::Face)) {
            bool swped=false;
            if (isVertex(v2)) {
                std::swap(v1,v2);
                swped=true;
            }

            // three point not in same edge
            if (fcommonidx!=-1) {
                size_t a = F(fcommonidx,0);
                size_t b = F(fcommonidx,1);
                size_t c = F(fcommonidx,2);
                info.removedFaces.insert(fcommonidx);
                info.addF.push_back(std::vector<size_t>{a,nV+increment,c});
                info.addF.push_back(std::vector<size_t>{c,nV+increment,b});
                info.addF.push_back(std::vector<size_t>{a,b,nV+increment});

                printf("add triangle type v-f/f-v - %lu, %lu,%lu; and %lu, %lu,%lu  and %lu, %lu,%lu. \n\n ",a,nV+increment,
                    c,c,nV+increment,b,a,b,nV+increment);
                if (swped) {increment++;}
                // if (info.addV.empty() || info.addV.back()!=toXYZ(v2)){info.addV.push_back(toXYZ(v2));}
            }
        }
        // 2. f-e
        if (v1.type == SurfacePointType::Face && v2.type == SurfacePointType::Edge) {

            for (int vidx = 0;vidx<3;vidx++) {
                if ( F(fcommonidx,vidx) != v2.edge.firstVertex().getIndex() &&
                     F(fcommonidx,vidx) != v2.edge.secondVertex().getIndex()  ) {
                    size_t a = F(fcommonidx,vidx);  // not on v1 edge
                    size_t b = F(fcommonidx,(vidx+1)%3);
                    size_t c = F(fcommonidx,(vidx+2)%3);

                    info.removedFaces.insert(fcommonidx);
                    info.addF.push_back(std::vector<size_t>{a,nV+increment,c});
                    info.addF.push_back(std::vector<size_t>{a,b,nV+increment});
                    printf("add triangle type f-e first part  - %lu, %lu,%lu; --  %lu, %lu,%lu. \n ",
                        a,nV+increment,c,a,b,nV+increment);

                    increment++;

                    info.addF.push_back(std::vector<size_t>{nV+increment-1,nV+increment,c});
                    info.addF.push_back(std::vector<size_t>{nV+increment-1,b,nV+increment});

                    printf("add triangle type f-e second part - %lu, %lu,%lu; --  %lu, %lu,%lu. \n\n ",
                        nV+increment-1,nV+increment,c,nV+increment,b,nV+increment);

                    break;
                }
            }

        }
        // 2. e-f
        if (v1.type == SurfacePointType::Edge && v2.type == SurfacePointType::Face) {

            for (int vidx = 0;vidx<3;vidx++) {
                if ( F(fcommonidx,vidx) != v1.edge.firstVertex().getIndex() &&
                     F(fcommonidx,vidx) != v1.edge.secondVertex().getIndex()  ) {
                    size_t a = F(fcommonidx,vidx);  // not on v1 edge
                    size_t b = F(fcommonidx,(vidx+1)%3);
                    size_t c = F(fcommonidx,(vidx+2)%3);
                    increment++;
                    info.removedFaces.insert(fcommonidx);
                    info.addF.push_back(std::vector<size_t>{a,nV+increment,c});
                    info.addF.push_back(std::vector<size_t>{a,b,nV+increment});
                    printf("add triangle type e-f first part  - %lu, %lu,%lu; --  %lu, %lu,%lu. \n ",
                        a,nV+increment,c,a,b,nV+increment);


                    info.addF.push_back(std::vector<size_t>{nV+increment-1,c,nV+increment});
                    info.addF.push_back(std::vector<size_t>{nV+increment,b,nV+increment-1});
                    // cntAddVertex ++;
                    printf("add triangle type e-f second part - %lu, %lu,%lu; --  %lu, %lu,%lu. \n\n ",
                        nV+increment-1,nV+increment,c,nV+increment-1,b,nV+increment);

                    break;
                }
            }

        }
        // 4 e-e
        if (v1.type == SurfacePointType::Edge && v2.type == SurfacePointType::Edge) {
            // three point not in same edge
            if (fcommonidx!=-1) {

                for (int vidx = 0;vidx<3;vidx++) {
                    if ( F(fcommonidx,vidx) != v1.edge.firstVertex().getIndex() &&
                     F(fcommonidx,vidx) != v1.edge.secondVertex().getIndex()  ) {
                        size_t a = F(fcommonidx,vidx);
                        size_t b = F(fcommonidx,(vidx+1)%3);
                        size_t c = F(fcommonidx,(vidx+2)%3);
                        info.removedFaces.insert(fcommonidx);

                        if (b== v2.edge.firstVertex().getIndex() || b== v2.edge.secondVertex().getIndex()) {
                            info.addF.push_back(std::vector<size_t>{a,nV+increment,c});
                            printf("add triangle type e-e-1 add first - %lu, %lu,%lu;\n ",a,nV+increment,c);
                            increment++;

                            info.addF.push_back(std::vector<size_t>{a,nV+increment,nV+increment-1});
                            info.addF.push_back(std::vector<size_t>{nV+increment,b,nV+increment-1});
                            printf("add triangle type e-e-1 add second part %lu, %lu,%lu, --  %lu, %lu,%lu \n\n ", a,nV+increment,nV+increment-1,nV+increment,b,nV+increment-1);
                            // increment++;
                        }
                        else {
                            info.addF.push_back(std::vector<size_t>{a,b,nV+increment});
                            printf("add triangle type e-e-2 add first - %lu, %lu,%lu;\n ",a,b,nV+increment);

                            increment++;
                            info.addF.push_back(std::vector<size_t>{a,nV+increment-1,nV+increment});
                            info.addF.push_back(std::vector<size_t>{nV+increment,nV+increment-1,c});
                            printf("add triangle type e-e-2 add second part %lu, %lu,%lu, --  %lu, %lu,%lu \n\n ", a,nV+increment-1,nV+increment,nV+increment,nV+increment-1,c);
                            // increment++;
                        }

                        break;
                    }
                }
            }
        }


    }
    for (auto f:info.removedFaces) {
        std::cout << f << std::endl;
        for (auto v:mesh->face(f).adjacentVertices()){std::cout<<v<<"- ";}
    }
    return info;
}
geometrycentral::Vector3 toXYZ(SurfacePoint p,ManifoldSurfaceMesh* mesh,VertexPositionGeometry* geometry) {
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
// ConnectupdateInfo updateConnection(ConnectupdateInfo info){
//     geometrycentral::Vector3 s = info.s;
//     geometrycentral::Vector3 t = info.t;
//     Eigen::MatrixXd V= info.V;
//     Eigen::MatrixXi F = info.F;
//     Eigen::MatrixXd newV;
//     Eigen::MatrixXi newF;
//     size_t nV = V.rows();
//     std::vector<std::vector<size_t>> addF;
//     std::unordered_set<size_t> removedFaces;
//     std::vector<geometrycentral::Vector3> addV;
//     ConnectupdateInfo nextinfo;
//     // return, finished
//     if ((s-t).norm2()<1e-6) {
//         // Eigen::MatrixXd newV = V;
//         // // std::cout << v << std::endl;
//         // newV.conservativeResize(newV.rows()+1, Eigen::NoChange);
//         // newV.row(newV.rows()-1)[0] = t.x;
//         // newV.row(newV.rows()-1)[1] = t.y;
//         // newV.row(newV.rows()-1)[2] = t.z;
//         // info.V = newV;
//         return info;
//     }
//     // do operation ============================
//
//     // point projection
//     Eigen::MatrixXd st(2,3);
//     st<<s.x,s.y,s.z,t.x,t.y,t.z;
//     Eigen::VectorXi I;
//     Eigen::MatrixXd C;
//     printf("next loop ====");
//
//     std::tie(C,I) = PointOnMesh(st,V,F);
//
//     Eigen::MatrixXd bc = ToBCCoordinates(V,F,C,I);
//     bcProcess(bc);
//
//     // use geometricial class
//     ManifoldSurfaceMesh* mesh;
//     VertexPositionGeometry* geometry;
//     std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr;
//     std::unique_ptr<VertexPositionGeometry> geometry_uptr;
//     std::tie(mesh_uptr, geometry_uptr) = makeManifoldSurfaceMeshAndGeometry(V, F);
//     mesh = mesh_uptr.release();
//     geometry = geometry_uptr.release();
//
//     GeodesicAlgorithmExact mmp(*mesh, *geometry);
//     SurfacePoint source = SurfacePoint(mesh->face(I[0]),geometrycentral::Vector3{bc(0,0),bc(0,1),bc(0,2)});
//     SurfacePoint target = SurfacePoint(mesh->face(I[1]),geometrycentral::Vector3{bc(1,0),bc(1,1),bc(1,2)});
//     mmp.propagate(target);
//     std::cout<<"Source: "<<source<< "target: "<<target <<std::endl;
//     std::vector<SurfacePoint> tmppath = mmp.traceBack(source);
//         // the path is s -> t, if s is vertex delete duplicate
//     if (source.faceCoords[0]==1 or source.faceCoords[1]==1 or source.faceCoords[2]==1) {tmppath.erase(tmppath.begin()+1);}
//     std::vector<SurfacePoint> path = robustPath2(tmppath,mesh,geometry);
//
//     SurfacePoint first = path[0];
//     SurfacePoint second = path[1];
//
//     // first is V
//     std::cout<<"first: "<<first<< "second: "<<second<<std::endl;
//     if (isVertex(first)) {
//         if (isVertex(second)) {
//             nextinfo.s =toXYZ(second,mesh,geometry) ;
//             nextinfo.t =t;
//             nextinfo.V = V;
//             nextinfo.F = F;
//             return updateConnection(nextinfo);
//         }
//         if (second.type == SurfacePointType::Edge) {
//             if (first.vertex.getIndex() == second.edge.firstVertex().getIndex() || first.vertex.getIndex() == second.edge.secondVertex().getIndex()) {
//                 nextinfo.s =toXYZ(second,mesh,geometry) ;
//                 nextinfo.t =t;
//                 nextinfo.V = V;
//                 nextinfo.F = F;
//                 return updateConnection(nextinfo);
//             }
//             else {
//                 int fcommonidx = findCommonFaceV2(first,second);
//                 printf("find common face: %d ",fcommonidx);
//                 for (int vidx = 0;vidx<3;vidx++) {
//                     if (F(fcommonidx,vidx) ==first.vertex.getIndex()) {
//                         size_t a = F(fcommonidx,vidx);
//                         size_t b = F(fcommonidx,(vidx+1)%3);
//                         size_t c = F(fcommonidx,(vidx+2)%3);
//                         removedFaces.insert(fcommonidx);
//                         addF.push_back(std::vector<size_t>{a,b,nV});
//                         addF.push_back(std::vector<size_t>{nV,c,a});
//                         addV.push_back(toXYZ(second,mesh,geometry));
//                         printf("add triangle type v-e :  %lu, %lu,%lu; and %lu, %lu,%lu. \n\n",a,b,nV,nV,c,a);
//                         break;
//                     }
//                 }
//                 std::tie(newV,newF) = updateMesh(V,F,addF,removedFaces,addV);
//                 nextinfo.s =toXYZ(second,mesh,geometry) ;
//                 nextinfo.t =t;
//                 nextinfo.V = newV;
//                 nextinfo.F = newF;
//                 // for (auto v:addV)   {std::cout<<v<<std::endl;}
//
//                 return updateConnection(nextinfo);
//                 }
//             }
//         if (second.type == SurfacePointType::Face) {
//             int fcommonidx = findCommonFaceV2(first,second);
//             for (int vidx = 0;vidx<3;vidx++) {
//                 size_t a = F(fcommonidx,vidx);
//                 size_t b = F(fcommonidx,(vidx+1)%3);
//                 size_t c = F(fcommonidx,(vidx+2)%3);
//                 removedFaces.insert(fcommonidx);
//                 addF.push_back(std::vector<size_t>{a,b,nV});
//                 addF.push_back(std::vector<size_t>{nV,b,c});
//                 addF.push_back(std::vector<size_t>{a,nV,c});
//                 // addV.push_back(t);
//                 printf("add triangle type v-f :  %lu, %lu,%lu; -- %lu, %lu,%lu;-- %lu, %lu,%lu \n\n",a,b,nV,nV,b,c,a,nV,c);
//                 std::tie(newV,newF) = updateMesh(V,F,addF,removedFaces,addV);
//                 nextinfo.s =toXYZ(second,mesh,geometry) ;
//                 nextinfo.t =t;
//                 nextinfo.V = newV;
//                 nextinfo.F = newF;
//                 // for (auto v:addV)   {std::cout<<v<<std::endl;}
//                 return updateConnection(nextinfo);
//             }
//
//         }
//
//         }
//     if (first.type == SurfacePointType::Edge) {
//
//         // addV.push_back(toXYZ(first,mesh,geometry));
//         if (second.type == SurfacePointType::Vertex) {
//             if (second.vertex.getIndex() == first.edge.firstVertex().getIndex() || second.vertex.getIndex() == first.edge.secondVertex().getIndex()) {
//
//                 std::tie(newV,newF) = updateMesh(V,F,addF,removedFaces,addV);
//                 nextinfo.s =toXYZ(second,mesh,geometry) ;
//                 nextinfo.t =t;
//                 nextinfo.V = V;
//                 nextinfo.F = F;
//                 return updateConnection(nextinfo);
//             }
//             else {
//                 int fcommonidx = findCommonFaceV2(first,second);
//                 printf("find common face e-v: %d ",fcommonidx);
//                 for (int vidx = 0;vidx<3;vidx++) {
//                     if (F(fcommonidx,vidx) ==second.vertex.getIndex()) {
//                         size_t a = F(fcommonidx,vidx);
//                         size_t b = F(fcommonidx,(vidx+1)%3);
//                         size_t c = F(fcommonidx,(vidx+2)%3);
//                         removedFaces.insert(fcommonidx);
//                         addF.push_back(std::vector<size_t>{a,b,nV});
//                         addF.push_back(std::vector<size_t>{nV,c,a});
//                         printf("add triangle type v-e :  %lu, %lu,%lu; and %lu, %lu,%lu. \n\n",a,b,nV,nV,c,a);
//                         break;
//                     }
//                 }
//                 std::tie(newV,newF) = updateMesh(V,F,addF,removedFaces,addV);
//                 nextinfo.s =toXYZ(second,mesh,geometry) ;
//                 nextinfo.t =t;
//                 nextinfo.V = newV;
//                 nextinfo.F = newF;
//                 // for (auto v:addV)   {std::cout<<v<<std::endl;}
//                 return updateConnection(nextinfo);
//             }
//         }
//         if (second.type == SurfacePointType::Face) {
//             int fcommonidx = findCommonFaceV2(first,second);
//             for (int vidx = 0;vidx<3;vidx++) {
//                 if ( F(fcommonidx,vidx) != first.edge.firstVertex().getIndex() &&
//                      F(fcommonidx,vidx) != first.edge.secondVertex().getIndex()  ) {
//                     size_t a = F(fcommonidx,vidx);  // not on v1 edge
//                     size_t b = F(fcommonidx,(vidx+1)%3);
//                     size_t c = F(fcommonidx,(vidx+2)%3);
//
//                     removedFaces.insert(fcommonidx);
//                     addF.push_back(std::vector<size_t>{a,b,nV+1});
//                     addF.push_back(std::vector<size_t>{a,nV+1,c});
//                     addF.push_back(std::vector<size_t>{nV+1,b,nV});
//                     addF.push_back(std::vector<size_t>{nV+1,nV,c});
//                     std::tie(newV,newF) = updateMesh(V,F,addF,removedFaces,addV);
//                     nextinfo.s =toXYZ(second,mesh,geometry) ;
//                     nextinfo.t =t;
//                     nextinfo.V = newV;
//                     nextinfo.F = newF;
//                     return updateConnection(nextinfo);
//                 }
//             }
//
//         }
//         if (second.type == SurfacePointType::Edge) {
//             int fcommonidx = findCommonFaceV2(first,second);
//             for (int vidx = 0;vidx<3;vidx++) {
//                     if ( F(fcommonidx,vidx) != first.edge.firstVertex().getIndex() &&
//                      F(fcommonidx,vidx) != first.edge.secondVertex().getIndex()  ) {
//                         size_t a = F(fcommonidx,vidx);
//                         size_t b = F(fcommonidx,(vidx+1)%3);
//                         size_t c = F(fcommonidx,(vidx+2)%3);
//                         removedFaces.insert(fcommonidx);
//                         if (b== second.edge.firstVertex().getIndex() || b== second.edge.secondVertex().getIndex()) {
//                             addF.push_back(std::vector<size_t>{a,nV,c});
//                             addF.push_back(std::vector<size_t>{a,nV+1,nV});
//                             addF.push_back(std::vector<size_t>{nV+1,b,nV});
//                         }
//                         else {
//                             addF.push_back(std::vector<size_t>{a,b,nV});
//                             addF.push_back(std::vector<size_t>{a,nV,nV+1});
//                             addF.push_back(std::vector<size_t>{nV+1,nV,c});
//                         }
//                         std::tie(newV,newF) = updateMesh(V,F,addF,removedFaces,addV);
//                         nextinfo.s =toXYZ(second,mesh,geometry) ;
//                         nextinfo.t =t;
//                         nextinfo.V = newV;
//                         nextinfo.F = newF;
//                         return updateConnection(nextinfo);
//
//                     }
//                 }
//
//         }
//
//     }
//
//
// }}
// void EasyDebugPrint3(int fi,size_t v1,size_t v2,size_t v3,size_t t1,size_t t2,size_t t3,
//     size_t t4,size_t t5,size_t t6,size_t t7,size_t t8,size_t t9) {
//     printf("remove face: %d (%lu,%lu,%lu), add triangle(f-v) :  %lu, %lu,%lu -- %lu, %lu,%lu-- %lu, %lu,%lu; \n",
//                 fi,v1,v2,v3,t1,t2,t3,t4,t5,t6,t7,t8,t9);
// }
// void EasyDebugPrint2(int fi,size_t v1,size_t v2,size_t v3,size_t t1,size_t t2,size_t t3,
//     size_t t4,size_t t5,size_t t6,size_t t7,size_t t8,size_t t9) {
//     printf("remove face: %d (%lu,%lu,%lu), add triangle(f-v) :  %lu, %lu,%lu -- %lu, %lu,%lu-- %lu, %lu,%lu; \n",
//                 fi,v1,v2,v3,t1,t2,t3,t4,t5,t6,t7,t8,t9);
// }
PathInfo updateConnection(PathInfo info){
    Eigen::MatrixXd V= info.V;
    Eigen::MatrixXi F = info.F;
    std::vector<SurfacePoint> inputPath = info.inputPath;

    Eigen::MatrixXd newV;
    Eigen::MatrixXi newF;
    size_t nV = V.rows();
    PathInfo afterUpdate;
    std::vector<std::vector<size_t>> addF;
    std::unordered_set<size_t> removedFaces;
    std::vector<geometrycentral::Vector3> addV;
    SurfacePoint first,second;

    // printf(" # of V: %lu \n", V.rows());

    // do operation ============================
    // use geometricial class
    ManifoldSurfaceMesh* mesh = info.mesh;
    VertexPositionGeometry* geometry = info.geometry;

    std::vector<SurfacePoint> path = robustPath4(inputPath,mesh,geometry);

    // std::vector<SurfacePoint> path = inputPath;

    // printf("============ path use for update;============== \n");
    // for (auto p:path) {
    //     std::cout<<p<<std::endl;
    //     if (p.type == SurfacePointType::Edge) {std::cout<<p.edge.firstVertex()<<"--"<<p.edge.secondVertex()<<std::endl;}
    //     else if (p.type == SurfacePointType::Face) {
    //         printf("remove face: %d (%lu,%lu,%lu) ",p.face.getIndex(),
    //             F(p.face.getIndex(),0),F(p.face.getIndex(),1),F(p.face.getIndex(),2));
    //     }
    // }
    // printf("============ path use for update;============== \n");
    // printf("path length： %lu \n",path.size());

    std::unordered_map<size_t,size_t> addVindexMap;
    int increment = 0;
    for (int i = 0; i < path.size(); i++) {
        SurfacePoint p = path[i];
        if (p.type == SurfacePointType::Edge || p.type == SurfacePointType::Face) {
            addV.push_back(toXYZ(p,mesh,geometry));
            addVindexMap[i] = nV+increment;
            increment++;
            afterUpdate.vLink.push_back(addVindexMap[i]);
        }
        else {
            afterUpdate.vLink.push_back(p.vertex.getIndex());
        }
    }

    for (int pathidx = 0; pathidx < path.size()-1; pathidx++) {
        first = path[pathidx];
        second = path[pathidx+1];
        // 1.  f--v   只会出现在开始
        if (first.type == SurfacePointType::Face && second.type == SurfacePointType::Vertex) {
            int opfaceIdx = first.face.getIndex();
            size_t a,b,c;
            a = F(opfaceIdx,0);
            b = F(opfaceIdx,1);
            c = F(opfaceIdx,2);
            addF.push_back(std::vector<size_t>{a,b,addVindexMap[pathidx]});
            addF.push_back(std::vector<size_t>{addVindexMap[pathidx],b,c});
            addF.push_back(std::vector<size_t>{a,addVindexMap[pathidx],c});
            // printf("remove face: %d (%lu,%lu,%lu), add triangle(f-v) :  %lu, %lu,%lu -- %lu, %lu,%lu-- %lu, %lu,%lu; \n",
            //     opfaceIdx,a,b,c,a,b,addVindexMap[pathidx],addVindexMap[pathidx],b,c,a,addVindexMap[pathidx],c);
            // EasyDebugPrint3(opfaceIdx,a,b,c,a,b,addVindexMap[pathidx],addVindexMap[pathidx],b,c,a,addVindexMap[pathidx],c);
            removedFaces.insert(opfaceIdx);

        }
        // 2.  f --> e   只会出现在开始
        if (first.type == SurfacePointType::Face && second.type == SurfacePointType::Edge) {
            int opfaceIdx = first.face.getIndex();
            size_t a,b,c;
            size_t ef  = second.edge.firstVertex().getIndex();
            size_t es = second.edge.secondVertex().getIndex();
            for (int vidx = 0;vidx<3;vidx++) {
                if ( F(opfaceIdx,vidx) != ef  && F(opfaceIdx,vidx) != es) {
                    a = F(opfaceIdx,vidx);
                    b = F(opfaceIdx,(vidx+1)%3);
                    c = F(opfaceIdx,(vidx+2)%3);
                    break;
                }
            }
            addF.push_back(std::vector<size_t>{a,b,addVindexMap[pathidx]});
            addF.push_back(std::vector<size_t>{a,addVindexMap[pathidx],c});
            // printf("remove face: %d (%lu,%lu,%lu), add triangle(f-e-first): %lu, %lu,%lu; -- %lu, %lu,%lu; \n",
            //     opfaceIdx,a,b,c,a,b,addVindexMap[pathidx], a,addVindexMap[pathidx],c);

            addF.push_back(std::vector<size_t>{addVindexMap[pathidx],b,addVindexMap[pathidx+1]});
            addF.push_back(std::vector<size_t>{addVindexMap[pathidx],addVindexMap[pathidx+1],c});
            // printf("remove face: %d (%lu,%lu,%lu), add triangle(f-e-second): %lu, %lu,%lu; -- %lu, %lu,%lu; \n",
            //     opfaceIdx,a,b,c,addVindexMap[pathidx],b,addVindexMap[pathidx+1],addVindexMap[pathidx],addVindexMap[pathidx+1],c);
            removedFaces.insert(opfaceIdx);

        }
        if (first.type == SurfacePointType::Vertex) {
            size_t firstvidx = first.vertex.getIndex();
            // 3. v-v do nothing
            if (second.type == SurfacePointType::Vertex) {}
            // 4. v --> e
            if (second.type == SurfacePointType::Edge) {
                size_t a,b,c;
                size_t ef  = second.edge.firstVertex().getIndex();
                size_t es = second.edge.secondVertex().getIndex();
                if (firstvidx == ef || firstvidx == es) {
                    printf("v--e happen, and v on e case.\n ");
                    // may do some operation
                    for (auto adface:second.edge.adjacentFaces()) {
                        int opfaceIdx = adface.getIndex();
                        removedFaces.insert(opfaceIdx);
                        for (int vidx = 0;vidx<3;vidx++) {
                            if ( F(opfaceIdx,vidx) != ef  && F(opfaceIdx,vidx) != es) {
                                a = F(opfaceIdx,vidx);
                                b = F(opfaceIdx,(vidx+1)%3);
                                c = F(opfaceIdx,(vidx+2)%3);
                                break;
                            }
                        }

                        addF.push_back(std::vector<size_t>{a,b,addVindexMap[pathidx+1]});
                        addF.push_back(std::vector<size_t>{a,addVindexMap[pathidx+1],c});
                        // printf("remove face: %d (%lu,%lu,%lu), add triangle: %lu, %lu,%lu; -- %lu, %lu,%lu;\n ",opfaceIdx,a,b,c,a,b,addVindexMap[pathidx+1],a,addVindexMap[pathidx+1],c);

                    }
                }
                else {
                    int opfaceIdx = findCommonFaceV2(first, second);
                    removedFaces.insert(opfaceIdx);
                    for (int vidx = 0;vidx<3;vidx++) {
                        if ( F(opfaceIdx,vidx) != ef  && F(opfaceIdx,vidx) != es){
                            a = F(opfaceIdx, vidx);
                            b = F(opfaceIdx, (vidx + 1) % 3);
                            c = F(opfaceIdx, (vidx + 2) % 3);
                            break;
                            }
                    }
                    // printf("traget triangle v to e: %lu, %lu, %lu \n",a,b,c);

                    addF.push_back(std::vector<size_t>{a,b,addVindexMap[pathidx+1]});
                    addF.push_back(std::vector<size_t>{a,addVindexMap[pathidx+1],c});
                    // printf("remove face %d (%lu,%lu,%lu), add triangle(v-e): %lu, %lu,%lu; %lu, %lu,%lu.\n",
                    //     opfaceIdx,a,b,c,a,b,addVindexMap[pathidx+1],a,addVindexMap[pathidx+1],c);

                }
            }
            // 5. v--> f 只可能是最后一个
            if (second.type == SurfacePointType::Face) {
                int opfaceIdx = second.face.getIndex();
                size_t a,b,c;
                a = F(opfaceIdx,0);
                b = F(opfaceIdx,1);
                c = F(opfaceIdx,2);
                addF.push_back(std::vector<size_t>{a,b,addVindexMap[pathidx+1]});
                addF.push_back(std::vector<size_t>{addVindexMap[pathidx+1],b,c});
                addF.push_back(std::vector<size_t>{a,addVindexMap[pathidx+1],c});
                // printf("remove face %d (%lu,%lu,%lu), add triangle type v-->f  :  %lu, %lu,%lu -- %lu, %lu,%lu-- %lu, %lu,%lu; \n",
                //     opfaceIdx,a,b,c,a,b,addVindexMap[pathidx+1],addVindexMap[pathidx+1],b,c,a,addVindexMap[pathidx+1],c);
                removedFaces.insert(opfaceIdx);
            }
        } // end of first.type == Vertex

        if (first.type == SurfacePointType::Edge) {
            // 6 e --> v

            if (second.type == SurfacePointType::Vertex) {
                size_t a,b,c;
                size_t secondvidx = second.vertex.getIndex();
                size_t ef  = first.edge.firstVertex().getIndex();
                size_t es = first.edge.secondVertex().getIndex();
                if (secondvidx == ef || secondvidx == es) {
                    // may do some operation
                    // printf("========= e-->v situation happen, v is on e ! =========\n");
                    for (auto adface:first.edge.adjacentFaces()) {
                        int opfaceIdx = adface.getIndex();
                        if (removedFaces.find(opfaceIdx) != removedFaces.end()) {continue;}
                        removedFaces.insert(opfaceIdx);
                        for (int vidx = 0;vidx<3;vidx++) {
                            if ( F(opfaceIdx,vidx) != ef  && F(opfaceIdx,vidx) != es) {
                                a = F(opfaceIdx,vidx);
                                b = F(opfaceIdx,(vidx+1)%3);
                                c = F(opfaceIdx,(vidx+2)%3);
                                break;
                            }
                        }

                        addF.push_back(std::vector<size_t>{a,b,addVindexMap[pathidx]});
                        addF.push_back(std::vector<size_t>{a,addVindexMap[pathidx],c});
                        // printf("remove face: %d (%lu,%lu,%lu), add triangle: %lu, %lu,%lu; -- %lu, %lu,%lu; \n",opfaceIdx,a,b,c,a,b,addVindexMap[pathidx],a,addVindexMap[pathidx],c);

                    }

                }
                else {
                    int opfaceIdx = findCommonFaceV2(first, second);
                    removedFaces.insert(opfaceIdx);
                    for (int vidx = 0;vidx<3;vidx++) {
                        if ( F(opfaceIdx,vidx) != ef  && F(opfaceIdx,vidx) != es) {
                            a = F(opfaceIdx,vidx);
                            b = F(opfaceIdx,(vidx+1)%3);
                            c = F(opfaceIdx,(vidx+2)%3);
                            break;
                        }
                    }

                    // printf("corresponding e to v: %lu, %lu, %lu \n",a,b,c);
                    addF.push_back(std::vector<size_t>{a,b,addVindexMap[pathidx]});
                    addF.push_back(std::vector<size_t>{a,addVindexMap[pathidx],c});
                    // printf("remove face %d (%lu,%lu,%lu), add triangle(e-v): %lu, %lu,%lu; and %lu, %lu,%lu. \n",
                    //     opfaceIdx,a,b,c,a,b,addVindexMap[pathidx],a,addVindexMap[pathidx],c);

                }
            }
            // 7 e --> e
            if (second.type == SurfacePointType::Edge) {

                int opfaceIdx = findCommonFaceV2(first, second);
                // std::cout << opfaceIdx << std::endl;

                size_t a,b,c;

                size_t firstef  = first.edge.firstVertex().getIndex();
                size_t firstes = first.edge.secondVertex().getIndex();

                size_t secondef  = second.edge.firstVertex().getIndex();
                size_t secondes = second.edge.secondVertex().getIndex();
                // printf("arrive here \n, common face %lu",opfaceIdx);
                for (int vidx = 0;vidx<3;vidx++) {
                    if ( F(opfaceIdx,vidx) != firstef  && F(opfaceIdx,vidx) != firstes) {
                        a = F(opfaceIdx,vidx);
                        b = F(opfaceIdx,(vidx+1)%3);
                        c = F(opfaceIdx,(vidx+2)%3);
                        break;
                    }
                }

                removedFaces.insert(opfaceIdx);

                // printf("remove face %d (%lu,%lu,%lu). ",opfaceIdx,a,b,c);
                if (b== secondes | b== secondef) {

                    addF.push_back(std::vector<size_t>{a,addVindexMap[pathidx],c});
                    // printf("add triangle(e-e-1 first): %lu, %lu,%lu;\n ",a,addVindexMap[pathidx],c);

                    addF.push_back(std::vector<size_t>{a,addVindexMap[pathidx+1],addVindexMap[pathidx]});
                    addF.push_back(std::vector<size_t>{addVindexMap[pathidx+1],b,addVindexMap[pathidx]});
                    // printf("add triangle(e-e-1 second) %lu, %lu,%lu, --  %lu, %lu,%lu \n\n ", a,addVindexMap[pathidx+1],addVindexMap[pathidx],
                    //     addVindexMap[pathidx+1],b,addVindexMap[pathidx]);
                }
                else {
                    // printf("corresponding v to edge-2: %lu,other two %lu, %lu \n",a,b,c);
                    addF.push_back(std::vector<size_t>{a,b,addVindexMap[pathidx]});
                    // printf("add triangle(e-e-2 first) - %lu, %lu,%lu;\n ",a,b,addVindexMap[pathidx]);

                    addF.push_back(std::vector<size_t>{a,addVindexMap[pathidx],addVindexMap[pathidx+1]});
                    addF.push_back(std::vector<size_t>{addVindexMap[pathidx+1],addVindexMap[pathidx],c});
                    // printf("add triangle(e-e-2 second) %lu, %lu,%lu, --  %lu, %lu,%lu \n\n ", a,addVindexMap[pathidx],addVindexMap[pathidx+1],
                    //     addVindexMap[pathidx+1],addVindexMap[pathidx],c);
                }
            }
            // 8. e--f // 只会在最后出现
            if (second.type == SurfacePointType::Face) {
                int opfaceIdx = second.face.getIndex();
                size_t a,b,c;
                size_t ef  = first.edge.firstVertex().getIndex();
                size_t es = first.edge.secondVertex().getIndex();
                for (int vidx = 0;vidx<3;vidx++) {
                    if ( F(opfaceIdx,vidx) != ef  && F(opfaceIdx,vidx) != es) {
                        a = F(opfaceIdx,vidx);
                        b = F(opfaceIdx,(vidx+1)%3);
                        c = F(opfaceIdx,(vidx+2)%3);
                        break;
                    }
                }
                addF.push_back(std::vector<size_t>{a,b,addVindexMap[pathidx+1]});
                addF.push_back(std::vector<size_t>{a,addVindexMap[pathidx+1],c});
                // printf("remove face: %d (%lu,%lu,%lu), add triangle(e-f-first):  %lu, %lu,%lu; -- %lu, %lu,%lu; \n",
                //     opfaceIdx,a,b,c,a,b,addVindexMap[pathidx+1], a,addVindexMap[pathidx+1],c);

                addF.push_back(std::vector<size_t>{addVindexMap[pathidx+1],b,addVindexMap[pathidx]});
                addF.push_back(std::vector<size_t>{addVindexMap[pathidx+1],addVindexMap[pathidx],c});
                // printf("remove face: %d (%lu,%lu,%lu), add triangle(e-f-second) :  %lu, %lu,%lu; -- %lu, %lu,%lu; \n",
                //     opfaceIdx,a,b,c,addVindexMap[pathidx+1],b,addVindexMap[pathidx],addVindexMap[pathidx+1],addVindexMap[pathidx],c);
                removedFaces.insert(opfaceIdx);
            }
        }

        }
    // if last point on edge
    SurfacePoint lastPoint = path.back();
    if (lastPoint.type == SurfacePointType::Edge) {
        size_t ef  = lastPoint.edge.firstVertex().getIndex();
        size_t es = lastPoint.edge.secondVertex().getIndex();
        size_t a,b,c;
        for (auto adface:lastPoint.edge.adjacentFaces()) {
            int opfaceIdx = adface.getIndex();
            if (removedFaces.find(opfaceIdx) == removedFaces.end()) {
                removedFaces.insert(opfaceIdx);
                for (int vidx = 0;vidx<3;vidx++) {
                    if ( F(opfaceIdx,vidx) != ef  && F(opfaceIdx,vidx) != es) {
                        a = F(opfaceIdx,vidx);
                        b = F(opfaceIdx,(vidx+1)%3);
                        c = F(opfaceIdx,(vidx+2)%3);
                        break;
                    }
                }

                addF.push_back(std::vector<size_t>{a,b,addVindexMap[path.size()-1]});
                addF.push_back(std::vector<size_t>{a,addVindexMap[path.size()-1],c});
                // printf("Final remove face: %d (%lu,%lu,%lu), add triangle: %lu, %lu,%lu; -- %lu, %lu,%lu;\n ",opfaceIdx,a,b,c,a,b,
                //     addVindexMap[path.size()-1],a,addVindexMap[path.size()-1],c);

            }
        }

    }

    std::tie(newV,newF) = updateMesh(V,F,addF,removedFaces,addV);

    afterUpdate.V = newV;
    afterUpdate.F = newF;
    afterUpdate.inputPath = path;
    printf("the update finish! \n");

    if (isVertex(path.back())) {
        afterUpdate.endPoint = path.back().vertex.getIndex();
    }
    else {
        afterUpdate.endPoint = newV.rows()-1;
    }

    if (isVertex(path[0])) {
        afterUpdate.start = path[0].vertex.getIndex();
    }
    else {
        afterUpdate.start = nV;
    }
    return afterUpdate;
    }


