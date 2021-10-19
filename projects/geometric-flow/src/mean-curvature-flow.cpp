// Implement member functions for MeanCurvatureFlow class.
#include "mean-curvature-flow.h"
#include "geometrycentral/numerical/linear_solvers.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
MeanCurvatureFlow::MeanCurvatureFlow(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    // Build member variables: mesh, geometry
    mesh = inputMesh;
    geometry = inputGeo;
}

/*
 * Build the mean curvature flow operator.
 *
 * Input: The mass matrix <M> of the mesh, and the timestep <h>.
 * Returns: A sparse matrix representing the mean curvature flow operator.
 */
SparseMatrix<double> MeanCurvatureFlow::buildFlowOperator(const SparseMatrix<double>& M, double h) const {

    SparseMatrix<double> L = geometry->laplaceMatrix();
    return M + h*L;

}

/*
 * Performs mean curvature flow.
 *
 * Input: The timestep <h>.
 * Returns:
 */
void MeanCurvatureFlow::integrate(double h) {

    SparseMatrix<double> M = geometry->massMatrix();
    SparseMatrix<double> A = buildFlowOperator(M, h);
    DenseMatrix<double> poistions(mesh->nVertices(),3);

    int i=0;
    for(auto v : mesh->vertices()){
        auto pos = geometry->inputVertexPositions[v];
        poistions.row(i) << pos.x , pos.y, pos.z;
        ++i;
    }

    Vector<double> x_0 = M*poistions.col(0);
    Vector<double> x_h = geometrycentral::solvePositiveDefinite(A, x_0);

    Vector<double> y_0 = M*poistions.col(1);
    Vector<double> y_h = geometrycentral::solvePositiveDefinite(A, y_0);

    Vector<double> z_0 = M*poistions.col(2);
    Vector<double> z_h = geometrycentral::solvePositiveDefinite(A, z_0);

    i=0;
    for (Vertex v : mesh->vertices()) {
        Vector3 pos_h{x_h[i], y_h[i], z_h[i]};
        geometry->inputVertexPositions[v] = pos_h;
        ++i;
    }
}