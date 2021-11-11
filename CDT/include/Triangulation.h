/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

/**
 * @file
 * Triangulation class
 */

#ifndef CDT_vW1vZ0lO8rS4gY4uI4fB
#define CDT_vW1vZ0lO8rS4gY4uI4fB

#include "CDTUtils.h"
#include "LocatorKDTree.h"

#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <memory>
#include <stack>
#include <vector>

namespace CDT
{

/**
 * Enum of strategies specifying order in which a range of vertices is inserted
 * @note @ref VertexInsertionOrder::Randomized will only randomize order of
 * inserting in triangulation, vertex indices will be preserved as they were
 * specified in the final triangulation
 */
struct CDT_EXPORT VertexInsertionOrder
{
    /**
     * The Enum itself
     * @note needed to pre c++11 compilers that don't support 'class enum'
     */
    enum Enum
    {
        Randomized, ///< vertices will be inserted in random order
        AsProvided, ///< vertices will be inserted in the same order as provided
    };
};

/// Enum of what type of geometry used to embed triangulation into
struct CDT_EXPORT SuperGeometryType
{
    /**
     * The Enum itself
     * @note needed to pre c++11 compilers that don't support 'class enum'
     */
    enum Enum
    {
        SuperTriangle, ///< conventional super-triangle
        Custom,        ///< user-specified custom geometry (e.g., grid)
    };
};

/**
 * Data structure representing a 2D constrained Delaunay triangulation
 *
 * @tparam T type of vertex coordinates (e.g., float, double)
 * @tparam TNearPointLocator class providing locating near point for efficiently
 * inserting new points. Provides methods: 'addPoint(vPos, iV)' and
 * 'nearPoint(vPos) -> iV'
 */
template <typename T, typename TNearPointLocator = LocatorKDTree<T> >
class CDT_EXPORT Triangulation
{
public:
    typedef std::vector<V2d<T> > V2dVec;              ///< Vertices vector
    typedef std::vector<TriIndVec> VerticesTriangles; ///< Triangles by vertex
    V2dVec vertices;            ///< triangulation's vertices
    TriangleVec triangles;      ///< triangulation's triangles
    EdgeUSet fixedEdges;        ///<  triangulation's constraints (fixed edges)
    VerticesTriangles vertTris; ///< triangles adjacent to each vertex

    /** Stores count of overlapping boundaries for a fixed edge. If no entry is
     * present for an edge: no boundaries overlap.
     * @note map only has entries for fixed for edges that represent overlapping
     * boundaries
     * @note needed for handling depth calculations and hole-removel in case of
     * overlapping boundaries
     */
    unordered_map<Edge, BoundaryOverlapCount> overlapCount;

    /*____ API _____*/
    /// Default constructor
    Triangulation();
    /**
     * Constructor
     * @param vertexInsertionOrder strategy used for ordering vertex insertions
     */
    Triangulation(VertexInsertionOrder::Enum vertexInsertionOrder);
    /**
     * Constructor
     * @param vertexInsertionOrder strategy used for ordering vertex insertions
     * @param nearPtLocator class providing locating near point for efficiently
     * inserting new points
     */
    Triangulation(
        VertexInsertionOrder::Enum vertexInsertionOrder,
        const TNearPointLocator& nearPtLocator);
    /**
     * Insert custom point-types specified by iterator range and X/Y-getters
     * @tparam TVertexIter iterator that dereferences to custom point type
     * @tparam TGetVertexCoordX function object getting x coordinate from
     * vertex. Getter signature: const TVertexIter::value_type& -> T
     * @tparam TGetVertexCoordY function object getting y coordinate from
     * vertex. Getter signature: const TVertexIter::value_type& -> T
     * @param first beginning of the range of vertices to add
     * @param last end of the range of vertices to add
     * @param getX getter of X-coordinate
     * @param getY getter of Y-coordinate
     */
    template <
        typename TVertexIter,
        typename TGetVertexCoordX,
        typename TGetVertexCoordY>
    void insertVertices(
        TVertexIter first,
        TVertexIter last,
        TGetVertexCoordX getX,
        TGetVertexCoordY getY);
    /// Insert vertices into triangulation
    void insertVertices(const std::vector<V2d<T> >& vertices);
    /**
     * Insert constraints (custom-type fixed edges) into triangulation
     * @note If some edge appears more than once in @ref insertEdges input this
     * means that multiple boundaries overlap at the edge and impacts
     * how hole detection algorithm of @ref eraseOuterTrianglesAndHoles works.
     * <b>Make sure there are no erroneous duplicates.</b>
     * @tparam TEdgeIter iterator that dereferences to custom edge type
     * @tparam TGetEdgeVertexStart function object getting start vertex index
     * from an edge.
     * Getter signature: const TEdgeIter::value_type& -> CDT::VertInd
     * @tparam TGetEdgeVertexEnd function object getting end vertex index from
     * an edge. Getter signature: const TEdgeIter::value_type& -> CDT::VertInd
     * @param first beginning of the range of edges to add
     * @param last end of the range of edges to add
     * @param getStart getter of edge start vertex index
     * @param getEnd getter of edge end vertex index
     */
    template <
        typename TEdgeIter,
        typename TGetEdgeVertexStart,
        typename TGetEdgeVertexEnd>
    void insertEdges(
        TEdgeIter first,
        TEdgeIter last,
        TGetEdgeVertexStart getStart,
        TGetEdgeVertexEnd getEnd);
    /**
     * Insert constraints (fixed edges) into triangulation
     * @note If some edge appears more than once in @ref insertEdges input this
     * means that multiple boundaries overlap at the edge and impacts
     * how hole detection algorithm of @ref eraseOuterTrianglesAndHoles works.
     * <b>Make sure there are no erroneous duplicates.</b>
     */
    void insertEdges(const std::vector<Edge>& edges);
    /**
     * Erase triangles adjacent to super triangle
     *
     * @note does nothing if custom geometry is used
     */
    void eraseSuperTriangle();
    /// Erase triangles outside of constrained boundary using growing
    void eraseOuterTriangles();
    /**
     * Erase triangles outside of constrained boundary and auto-detected holes
     *
     * @note detecting holes relies on layer peeling based on layer depth
     * @note supports overlapping or touching boundaries
     */
    void eraseOuterTrianglesAndHoles();
    /**
     * Call this method after directly setting custom super-geometry via
     * vertices and triangles members
     */
    void initializedWithCustomSuperGeometry();

private:
    /*____ Detail __*/
    void addSuperTriangle(const Box2d<T>& box);
    void addNewVertex(const V2d<T>& pos, const TriIndVec& tris);
    void insertVertex(const VertInd iVert);
    void insertEdge(Edge edge);
    tuple<TriInd, VertInd, VertInd> intersectedTriangle(
        const VertInd iA,
        const std::vector<TriInd>& candidates,
        const V2d<T>& a,
        const V2d<T>& b) const;
    /// Returns indices of three resulting triangles
    std::stack<TriInd> insertPointInTriangle(const VertInd v, const TriInd iT);
    /// Returns indices of four resulting triangles
    std::stack<TriInd>
    insertPointOnEdge(const VertInd v, const TriInd iT1, const TriInd iT2);
    array<TriInd, 2> trianglesAt(const V2d<T>& pos) const;
    array<TriInd, 2> walkingSearchTrianglesAt(const V2d<T>& pos) const;
    TriInd walkTriangles(const VertInd startVertex, const V2d<T>& pos) const;
    bool isFlipNeeded(
        const V2d<T>& pos,
        const TriInd iT,
        const TriInd iTopo,
        const VertInd iVert) const;
    void flipEdge(const TriInd iT, const TriInd iTopo);
    void changeNeighbor(
        const TriInd iT,
        const TriInd oldNeighbor,
        const TriInd newNeighbor);
    void changeNeighbor(
        const TriInd iT,
        const VertInd iVedge1,
        const VertInd iVedge2,
        const TriInd newNeighbor);
    void addAdjacentTriangle(const VertInd iVertex, const TriInd iTriangle);
    void addAdjacentTriangles(
        const VertInd iVertex,
        const TriInd iT1,
        const TriInd iT2,
        const TriInd iT3);
    void addAdjacentTriangles(
        const VertInd iVertex,
        const TriInd iT1,
        const TriInd iT2,
        const TriInd iT3,
        const TriInd iT4);
    void removeAdjacentTriangle(const VertInd iVertex, const TriInd iTriangle);
    TriInd triangulatePseudopolygon(
        const VertInd ia,
        const VertInd ib,
        const std::vector<VertInd>& points);
    VertInd findDelaunayPoint(
        const VertInd ia,
        const VertInd ib,
        const std::vector<VertInd>& points) const;
    TriInd pseudopolyOuterTriangle(const VertInd ia, const VertInd ib) const;
    TriInd addTriangle(const Triangle& t); // note: invalidates iterators!
    TriInd addTriangle(); // note: invalidates triangle iterators!
    void makeDummy(const TriInd iT);
    void eraseDummies();
    void eraseSuperTriangleVertices(); // no effect if custom geometry is used
    template <typename TriIndexIter>
    void eraseTrianglesAtIndices(TriIndexIter first, TriIndexIter last);
    TriIndUSet growToBoundary(std::stack<TriInd> seeds) const;
    void fixEdge(const Edge& edge);

    std::vector<TriInd> m_dummyTris;
    TNearPointLocator m_nearPtLocator;
    std::size_t m_nTargetVerts;
    SuperGeometryType::Enum m_superGeomType;
    VertexInsertionOrder::Enum m_vertexInsertionOrder;
};

} // namespace CDT

namespace CDT
{

inline size_t randomCDT(const size_t i)
{
    static mt19937 g(9001);
    return g() % i;
}

//-----------------------
// Triangulation methods
//-----------------------
template <typename T, typename TNearPointLocator>
template <
    typename TVertexIter,
    typename TGetVertexCoordX,
    typename TGetVertexCoordY>
void Triangulation<T, TNearPointLocator>::insertVertices(
    const TVertexIter first,
    const TVertexIter last,
    TGetVertexCoordX getX,
    TGetVertexCoordY getY)
{
    if(vertices.empty())
    {
        addSuperTriangle(envelopBox<T>(first, last, getX, getY));
    }

    const std::size_t nExistingVerts = vertices.size();

    vertices.reserve(nExistingVerts + std::distance(first, last));
    for(TVertexIter it = first; it != last; ++it)
        addNewVertex(V2d<T>::make(getX(*it), getY(*it)), TriIndVec());

    switch(m_vertexInsertionOrder)
    {
    case VertexInsertionOrder::AsProvided:
        for(TVertexIter it = first; it != last; ++it)
            insertVertex(VertInd(nExistingVerts + std::distance(first, it)));
        break;
    case VertexInsertionOrder::Randomized:
        std::vector<VertInd> ii(std::distance(first, last));
        typedef std::vector<VertInd>::iterator Iter;
        VertInd value = nExistingVerts;
        for(Iter it = ii.begin(); it != ii.end(); ++it, ++value)
            *it = value;
        std::random_shuffle(ii.begin(), ii.end(), randomCDT);
        for(Iter it = ii.begin(); it != ii.end(); ++it)
            insertVertex(*it);
        break;
    }
}

template <typename T, typename TNearPointLocator>
template <
    typename TEdgeIter,
    typename TGetEdgeVertexStart,
    typename TGetEdgeVertexEnd>
void Triangulation<T, TNearPointLocator>::insertEdges(
    TEdgeIter first,
    const TEdgeIter last,
    TGetEdgeVertexStart getStart,
    TGetEdgeVertexEnd getEnd)
{
    for(; first != last; ++first)
    {
        // +3 to account for super-triangle vertices
        insertEdge(Edge(
            VertInd(getStart(*first) + m_nTargetVerts),
            VertInd(getEnd(*first) + m_nTargetVerts)));
    }
    eraseDummies();
}

} // namespace CDT

#ifndef CDT_USE_AS_COMPILED_LIBRARY
#include "Triangulation.hpp"
#endif

#endif // header-guard
