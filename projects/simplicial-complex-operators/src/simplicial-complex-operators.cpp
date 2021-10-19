// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        geometry->vertexIndices[v] = idx;
        ++idx;
    }

    idx = 0;
    for(auto e : mesh->edges()){
        geometry->edgeIndices[e] = idx;
        ++idx;
    }

    idx = 0;
    for(auto f : mesh->faces()){
        geometry->faceIndices[f] = idx;
        ++idx;
    }
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    // SparseMatrix elements are i,j,value triplets
    typedef Eigen::Triplet<size_t> T;
    std::vector<T> triplets;

    for(auto e : mesh->edges()){
        triplets.push_back(T(e.getIndex(), e.firstVertex().getIndex(), 1));
        triplets.push_back(T(e.getIndex(), e.secondVertex().getIndex(), 1));
    }

    Eigen::SparseMatrix<size_t> sparseMatrix(mesh->nEdges(), mesh->nVertices());
    sparseMatrix.setFromTriplets(triplets.begin(), triplets.end());

    
    return sparseMatrix;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    // SparseMatrix elements are i,j,value triplets
    typedef Eigen::Triplet<size_t> T;
    std::vector<T> triplets;

    for(auto f : mesh->faces()){
        for(auto e : f.adjacentEdges()){
            triplets.push_back(T(f.getIndex(), e.getIndex(), 1));
        }
    }

    Eigen::SparseMatrix<size_t> sparseMatrix(mesh->nFaces(), mesh->nEdges());
    sparseMatrix.setFromTriplets(triplets.begin(), triplets.end());

    return sparseMatrix;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    Vector<size_t> vertices (mesh->nVertices());
    for(size_t i=0; i < mesh->nVertices(); i++){
        if(subset.vertices.count(i)){
            vertices[i]=1;
        }
        else{
            vertices[i]=0;
        }
    }
    return vertices;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    Vector<size_t> edges (mesh->nEdges());
    for(size_t i=0; i < mesh->nEdges(); i++){
        if(subset.edges.count(i)){
            edges[i]=1;
        }
        else{
            edges[i]=0;
        }
    }
    return edges;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    Vector<size_t> faces (mesh->nFaces());
    for(size_t i=0; i < mesh->nFaces(); i++){
        if(subset.edges.count(i)){
            faces[i]=1;
        }
        else{
            faces[i]=0;
        }
    }
    return faces;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    auto star(subset);

    for(auto v : subset.vertices){
        for(auto e : mesh->vertex(v).adjacentEdges()){
            star.edges.insert(e.getIndex());
        }
        for(auto f : mesh->vertex(v).adjacentFaces()){
            star.faces.insert(f.getIndex());
        }
    }

    for(auto e : subset.edges){
        for(auto f : mesh->edge(e).adjacentFaces()){
            star.faces.insert(f.getIndex());
        }
    }

    return star;
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    auto closure(subset);

    for(auto e : subset.edges){
        for(auto v : mesh->edge(e).adjacentVertices()){
            closure.vertices.insert(v.getIndex());
        }
    }

    for(auto f : subset.faces){
        for(auto v : mesh->face(f).adjacentVertices()){
            closure.vertices.insert(v.getIndex());
        }
        for(auto e : mesh->face(f).adjacentEdges()){
            closure.edges.insert(e.getIndex());
        }
    }

    return closure;
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    auto st = star(subset);
    auto cl = closure(subset);

    auto stcl = star(cl);
    auto clst = closure(st);

    auto sub = clst.deepCopy();
    sub.deleteSubset(stcl);

    return sub;
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    for(auto e : subset.edges){
        for(auto v : mesh->edge(e).adjacentVertices()){
            if(0 == subset.vertices.count(v.getIndex())){
                return false;
            }
        }
    }

    for(auto f : subset.faces){

        for(auto v : mesh->face(f).adjacentVertices()){
            if(0 == subset.vertices.count(v.getIndex())){
                return false;
            }
        }

        for(auto e : mesh->face(f).adjacentEdges()){
            if(0 == subset.edges.count(e.getIndex())){
                return false;
            }
        }

    }

    return true;
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    if(false == isComplex(subset)){
        return -1;
    }

    if(0 == subset.faces.size() && 0 == subset.edges.size()){
        return 0;
    }

    for(auto v : subset.vertices){

        bool contained = false;

        for(auto e : mesh->vertex(v).adjacentEdges()){
            if(subset.edges.count(e.getIndex())){
                contained = true;
                break;
            }
        }

        if(!contained){
            return -1;
        }
    }

    if(0 == subset.faces.size()){
        return 1;
    }

    for(auto e : subset.edges){

        bool contained = false;

        for(auto f : mesh->edge(e).adjacentFaces()){
            if(subset.faces.count(f.getIndex())){
                contained = true;
                break;
            }
        }

        if(!contained){
            return -1;
        }
    }

    return 2;
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    // TODO
    return subset; // placeholder
}