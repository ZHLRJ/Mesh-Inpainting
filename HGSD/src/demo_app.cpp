
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

int main(int argc, char **argv) {
    std::string filepath = "/Users/haoliang/Downloads/ongoing/crumpleddevelopable/CrumpledDevelopablepart.obj";
    std::string projMeshfile = "/Users/haoliang/Downloads/ongoing/result/poissonMesh.obj";
    Eigen::MatrixXd projmeshV;
    Eigen::MatrixXi projmeshF;
    igl::readOBJ(filepath, projmeshV, projmeshF);

    // use geometricial class
    std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr;
    std::unique_ptr<VertexPositionGeometry> geometry_uptr;
    std::tie(mesh_uptr, geometry_uptr) = makeManifoldSurfaceMeshAndGeometry(projmeshV, projmeshF);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();

    polyscope::init();

    std::vector<geometrycentral::Vector3> onTrace, MapPts;

    drawMesh("demoTestMesh", projmeshV, projmeshF);
    // test loop
    SurfacePoint p1,p2,p3,p4,p5,v1,v2,px,py;
    MeshModify testMesh(projmeshV,projmeshF);
    // Source: [SurfacePoint: type=Face, face= f_18115 faceCoords= <0, 0, 1>]
    // target: [SurfacePoint: type=Face, face= f_12846 faceCoords= <0, 0.367898, 0.632102>]

    // px = SurfacePoint(mesh->face(15887),geometrycentral::Vector3{0.350274, 0.277123, 0.372602});
    // py = SurfacePoint(mesh->face(15905),geometrycentral::Vector3{0.140207, 0, 0.859793});

    px = SurfacePoint(mesh->face(0),geometrycentral::Vector3{0.350274, 0.277123, 0.372602});
    py = SurfacePoint(mesh->face(10),geometrycentral::Vector3{0.2, 0.3, 0.5});
    // px = SurfacePoint(mesh->vertex(9045));
    // py = SurfacePoint(mesh->vertex(9020));
    p1 = SurfacePoint(testMesh.mesh->face(5),geometrycentral::Vector3{1.,0.,0.});
    p2 = SurfacePoint(testMesh.mesh->face(4),geometrycentral::Vector3{1.,0.,0.});
    p3 = SurfacePoint(testMesh.mesh->face(10),geometrycentral::Vector3{0.3,0.5,0.2});
    p4 = SurfacePoint(testMesh.mesh->face(16),geometrycentral::Vector3{1,0.,0.});
    // p4 = SurfacePoint(testMesh.mesh->vertex(3));
    p5 = SurfacePoint(testMesh.mesh->face(15),geometrycentral::Vector3{0.,0.5,0.5});
    std::vector<SurfacePoint> arr{px,px,p1,p2,p3,p4,p5};
    remeshInfo InfoAll;
    std::vector<SurfacePoint> Allpath;

    for (int i = 0; i < arr.size(); i++) {
        MapPts.push_back(toXYZ(arr[i],mesh,geometry));
    }
    Eigen::MatrixXd newV;Eigen::MatrixXi newF;
    PathInfo ans,input;
    input.V = projmeshV;
    input.F = projmeshF;
    input.mesh = mesh;
    input.geometry = geometry;
    for (auto i=0; i < 1; ++i) {

        printf(" hole : %d \n",i);


        GeodesicAlgorithmExact mmp(*mesh, *geometry);
        mmp.propagate(py);
        std::vector<SurfacePoint> path = mmp.traceBack(px);
        input.inputPath = path;

        ans =  updateConnection(input);
        printf("Finished . ");
        for (auto p:path) {
            // std::cout << p<<std::endl;
            onTrace.push_back(toXYZ(p,mesh,geometry));

        }
        printf("path finished");
    }

    // polyscope::CurveNetwork* trace = polyscope::registerCurveNetworkLine("onTrace", onTrace);
    polyscope::CurveNetwork* trace = polyscope::registerCurveNetworkLine("onTrace", onTrace);
    trace->setColor(glm::vec3(247./255,5./255,5./255));
    drawMesh("ansMesh", ans.V, ans.F,81,212,131);
    // drawMesh("newMesh", newV,newF,81,212,131);
    polyscope::show();


}
