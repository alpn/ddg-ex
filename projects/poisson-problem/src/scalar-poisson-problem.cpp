// Implement member functions for ScalarPoissonProblem class.
#include "scalar-poisson-problem.h"
#include "geometrycentral/numerical/linear_solvers.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
ScalarPoissonProblem::ScalarPoissonProblem(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;

    this->A = geometry->laplaceMatrix();
    this->M = geometry->massMatrix();
    this->totalArea = geometry->totalArea();
}

/*
 * Computes the solution of the poisson problem Ax = -M(rho - rhoBar), where A is the POSITIVE DEFINITE Laplace matrix
 * and M is the mass matrix.
 *
 * Input: <rho>, the density of vertices in the mesh.
 * Returns: The solution vector.
 */
Vector<double> ScalarPoissonProblem::solve(const Vector<double>& rho) const {

    SparseMatrix<double> matrix = this->A;

    assert(0!=totalArea);
    double rhoBarVal = ((M*rho) / totalArea).sum();
    Vector<double> rhoBar = Vector<double>::Constant(rho.size(), rhoBarVal);

    Vector<double> rhs = -M*(rho - rhoBar);
    return geometrycentral::solvePositiveDefinite(matrix, rhs);

}