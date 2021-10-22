#include "heat-method.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
HeatMethod::HeatMethod(ManifoldSurfaceMesh* surfaceMesh, VertexPositionGeometry* geo) {

    this->mesh = surfaceMesh;
    this->geometry = geo;

    geo->requireVertexPositions();

    this->A = geo->laplaceMatrix();

    SparseMatrix<double> M = geometry->massMatrix();
    double h = geo->meanEdgeLength();
    double dt = h*h;
    this->F = M + dt*this->A;
}

/*
 * Computes the vector field X = -∇u / |∇u|.
 *
 * Input: <u>, a dense vector representing the heat that is allowed to diffuse on the input mesh for a brief period of
 * time.
 * Returns: A MeshData container that stores a Vector3 per face.
 */
FaceData<Vector3> HeatMethod::computeVectorField(const Vector<double>& u) const {

    const auto u0 = u;
    SparseMatrix<double> FL = F;
    Vector<double> ut = geometrycentral::solvePositiveDefinite(FL, u0);

    auto faceData = FaceData<Vector3>(*mesh, {0, 0, 0});

    for(const auto& f : mesh->faces()){

        const auto he = f.halfedge();

        size_t v1 = he.tipVertex().getIndex();
        size_t v2 = he.tailVertex().getIndex();
        size_t v3 = he.next().tipVertex().getIndex();

        const auto e1 = geometry->vertexPositions[v3]-geometry->vertexPositions[v2];
        const auto e2 = geometry->vertexPositions[v1]-geometry->vertexPositions[v3];
        const auto e3 = geometry->vertexPositions[v2]-geometry->vertexPositions[v1];

        auto N = geometry->faceNormal(f);

        auto sum =
            u[v1]*(cross(N,e1)) +
            u[v2]*(cross(N,e2)) +
            u[v3]*(cross(N,e3));

        double denominator = 2 * geometry->faceArea(f);
        assert(0 != denominator);

        const auto preX = sum / denominator;

        if(preX.norm() > 0){

            const auto X = preX.normalize();
            faceData[f] = X;
        }
    }

    return faceData;
}

/*
 * Computes the integrated divergence ∇.X.
 *
 * Input: <X>, the vector field -∇u / |∇u| represented as a FaceData container
 * Returns: A dense vector
 */
Vector<double> HeatMethod::computeDivergence(const FaceData<Vector3>& X) const {

    Vector<double> result(mesh->nVertices());

    for(const auto v : mesh->vertices()){

        double integratedDivAtV = 0.0;

        for(const auto he1 : v.outgoingHalfedges()){

            const auto Xf = X[he1.face()];
            const auto he2 = he1.next().next();

            const auto e1 = geometry->vertexPositions[he1.tipVertex().getIndex()]-
                            geometry->vertexPositions[v.getIndex()];

            const auto e2 = geometry->vertexPositions[he2.tailVertex().getIndex()]-
                            geometry->vertexPositions[v.getIndex()];

            integratedDivAtV += 0.5 *
                            (geometry->cotan(he1) * dot(e1,Xf) +
                             geometry->cotan(he2) * dot(e2,Xf));
        }

        result[v.getIndex()] = integratedDivAtV;
    }

    return result;
}

/*
 * Computes the geodesic distances φ using the heat method.
 *
 * Input: <delta>, a dense vector representing the heat sources, i.e., u0 = δ(x). Returns: A dense vector containing the
 * geodesic distances per vertex.
 */
Vector<double> HeatMethod::compute(const Vector<double>& delta) const {

    // TODO
    Vector<double> phi = Vector<double>::Zero(delta.rows());

    // Since φ is unique up to an additive constant, it should be shifted such that the smallest distance is zero
    this->subtractMinimumDistance(phi);

    return phi;
}