/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

/**
 * @file
 * Public API
 */

#ifndef CDT_lNrmUayWQaIR5fxnsg9B
#define CDT_lNrmUayWQaIR5fxnsg9B

#include "CDTUtils.h"
#include "Triangulation.h"
#include "remove_at.hpp"

#include <algorithm>
#include <iterator>
#include <stack>
#include <utility>
#include <vector>

namespace CDT
{

/**
 * Information about removed duplicated vertices.
 *
 * Contains mapping information and removed duplicates indices.
 * @note vertices {0,1,2,3,4} where 0 and 3 are the same will produce mapping
 *       {0,1,2,0,3} (to new vertices {0,1,2,3}) and duplicates {3}
 */
struct CDT_EXPORT DuplicatesInfo
{
    std::vector<std::size_t> mapping;    ///< vertex index mapping
    std::vector<std::size_t> duplicates; ///< duplicates' indices
};

/**
 * Find duplicates in given custom point-type range
 * @note duplicates are points with exactly same X and Y coordinates
 * @tparam TVertexIter iterator that dereferences to custom point type
 * @tparam TGetVertexCoordX function object getting x coordinate from vertex.
 * Getter signature: const TVertexIter::value_type& -> T
 * @tparam TGetVertexCoordY function object getting y coordinate from vertex.
 * Getter signature: const TVertexIter::value_type& -> T
 * @param first beginning of the range of vertices
 * @param last end of the range of vertices
 * @param getX getter of X-coordinate
 * @param getY getter of Y-coordinate
 * @returns information about vertex duplicates
 */
template <
    typename T,
    typename TVertexIter,
    typename TGetVertexCoordX,
    typename TGetVertexCoordY>
CDT_EXPORT DuplicatesInfo FindDuplicates(
    TVertexIter first,
    TVertexIter last,
    TGetVertexCoordX getX,
    TGetVertexCoordY getY);

/**
 * Remove duplicates in-place from vector of custom points
 * @tparam TVertexIter iterator that dereferences to custom point type
 * @tparam TAllocator allocator used by input vector of vertices
 * @param vertices vertices to remove duplicates from
 * @param duplicates information about duplicates
 */
template <typename TVertex, typename TAllocator>
CDT_EXPORT void RemoveDuplicates(
    std::vector<TVertex, TAllocator>& vertices,
    const std::vector<std::size_t>& duplicates);

/**
 * Remove duplicated points in-place
 *
 * @tparam T type of vertex coordinates (e.g., float, double)
 * @param[in, out] vertices collection of vertices to remove duplicates from
 * @returns information about duplicated vertices that were removed.
 */
template <typename T>
CDT_EXPORT DuplicatesInfo RemoveDuplicates(std::vector<V2d<T> >& vertices);

/**
 * Remap vertex indices in edges (in-place) using given vertex-index mapping.
 *
 * @note Mapping can be a result of RemoveDuplicates function
 * @param[in,out] edges collection of edges to remap
 * @param mapping vertex-index mapping
 */
CDT_EXPORT void
RemapEdges(std::vector<Edge>& edges, const std::vector<std::size_t>& mapping);

/**
 * Find point duplicates, remove them from vector (in-place) and remap edges
 * (in-place)
 * @note Same as a chained call of @ref FindDuplicates, @ref RemoveDuplicates,
 * and @ref RemapEdges
 * @tparam TVertexIter iterator that dereferences to custom point type
 * @tparam TGetVertexCoordX function object getting x coordinate from vertex.
 * Getter signature: const TVertexIter::value_type& -> T
 * @tparam TGetVertexCoordY function object getting y coordinate from vertex.
 * Getter signature: const TVertexIter::value_type& -> T
 * @param[in, out] vertices vertices to remove duplicates from
 * @param[in, out] edges collection of edges connecting vertices
 * @param getX getter of X-coordinate
 * @param getY getter of Y-coordinate
 * @returns information about vertex duplicates
 */
template <
    typename T,
    typename TVertex,
    typename TGetVertexCoordX,
    typename TGetVertexCoordY,
    typename TVertexAllocator,
    typename TEdgeAllocator>
CDT_EXPORT DuplicatesInfo RemoveDuplicatesAndRemapEdges(
    std::vector<TVertex, TVertexAllocator>& vertices,
    std::vector<Edge, TEdgeAllocator>& edges,
    TGetVertexCoordX getX,
    TGetVertexCoordY getY);

/**
 * Same as a chained call of @ref RemoveDuplicates + @ref RemapEdges
 *
 * @tparam T type of vertex coordinates (e.g., float, double)
 * @param[in, out] vertices collection of vertices to remove duplicates from
 * @param[in,out] edges collection of edges to remap
 */
template <typename T>
CDT_EXPORT DuplicatesInfo RemoveDuplicatesAndRemapEdges(
    std::vector<V2d<T> >& vertices,
    std::vector<Edge>& edges);

/**
 * Calculate depth of each triangle in constraint triangulation.
 *
 * Perform depth peeling from super triangle to outermost boundary,
 * then to next boundary and so on until all triangles are traversed.
 * For example depth is:
 *  - 0 for triangles outside outermost boundary
 *  - 1 for triangles inside boundary but outside hole
 *  - 2 for triangles in hole
 *  - 3 for triangles in island and so on...
 *
 * @param seed seed triangle to begin depth-peeling from
 * @param triangles triangles of triangulation
 * @param fixedEdges constraint edges of triangulation
 * @return vector where element at index i stores depth of i-th triangle
 */
CDT_EXPORT std::vector<LayerDepth> CalculateTriangleDepths(
    const TriInd seed,
    const TriangleVec& triangles,
    const EdgeUSet& fixedEdges);

/**
 * Calculate depth of each triangle in constraint triangulation. Supports
 * overlapping boundaries.
 *
 * Perform depth peeling from super triangle to outermost boundary,
 * then to next boundary and so on until all triangles are traversed.@n
 * For example depth is:
 *  - 0 for triangles outside outermost boundary
 *  - 1 for triangles inside boundary but outside hole
 *  - 2 for triangles in hole
 *  - 3 for triangles in island and so on...
 *
 * @param seed seed triangle to begin depth-peeling from
 * @param triangles triangles of triangulation
 * @param fixedEdges constraint edges of triangulation
 * @param overlapCount counts of boundary overlaps at fixed edges (map has
 * entries only for edges that represent overlapping boundaries)
 * @return vector where element at index i stores depth of i-th triangle
 */
CDT_EXPORT std::vector<LayerDepth> CalculateTriangleDepths(
    const TriInd seed,
    const TriangleVec& triangles,
    const EdgeUSet& fixedEdges,
    const unordered_map<Edge, BoundaryOverlapCount>& overlapCount);

/**
 * Depth-peel a layer in triangulation, used when calculating triangle depths
 *
 * It takes starting seed triangles, traverses neighboring triangles, and
 * assigns given layer depth to the traversed triangles. Traversal is blocked by
 * constraint edges. Triangles behind constraint edges are recorded as seeds of
 * next layer and returned from the function.
 *
 * @tparam T type of vertex coordinates (e.g., float, double)
 * @param seeds indices of seed triangles
 * @param triangles triangles of triangulation
 * @param fixedEdges constraint edges of triangulation
 * @param layerDepth current layer's depth to mark triangles with
 * @param[in, out] triDepths depths of triangles
 * @return triangles of the next (deeper) layer that are adjacent to the peeled
 * layer. To be used as seeds when peeling the next layer.
 */
CDT_EXPORT
TriIndUSet PeelLayer(
    std::stack<TriInd> seeds,
    const TriangleVec& triangles,
    const EdgeUSet& fixedEdges,
    const LayerDepth layerDepth,
    std::vector<LayerDepth>& triDepths);

/**
 * Depth-peel a layer in triangulation, used when calculating triangle depths
 *
 * It takes starting seed triangles, traverses neighboring triangles, and
 * assigns given layer depth to the traversed triangles. Traversal is blocked by
 * constraint edges. Triangles behind constraint edges are recorded as seeds of
 * next layer and returned from the function.
 *
 * @tparam T type of vertex coordinates (e.g., float, double)
 * @param seeds indices of seed triangles
 * @param triangles triangles of triangulation
 * @param fixedEdges constraint edges of triangulation
 * @param overlapCount counts of boundary overlaps at fixed edges (map has
 * entries only for edges that represent overlapping boundaries)
 * @param layerDepth current layer's depth to mark triangles with
 * @param[in, out] triDepths depths of triangles
 * @return triangles of the deeper layers that are adjacent to the peeled layer.
 *         To be used as seeds when peeling deeper layers.
 */
CDT_EXPORT
unordered_map<TriInd, LayerDepth> PeelLayer(
    std::stack<TriInd> seeds,
    const TriangleVec& triangles,
    const EdgeUSet& fixedEdges,
    const unordered_map<Edge, BoundaryOverlapCount>& overlapCount,
    const LayerDepth layerDepth,
    std::vector<LayerDepth>& triDepths);

/**
 * Extract all edges of triangles
 *
 * @param triangles triangles used to extract edges
 * @return an unordered set of all edges of triangulation
 */
CDT_EXPORT EdgeUSet extractEdgesFromTriangles(const TriangleVec& triangles);

} // namespace CDT

//*****************************************************************************
// Implementations of template functionlity
//*****************************************************************************
// hash for CDT::V2d<T>
#ifdef CDT_CXX11_IS_SUPPORTED
namespace std
#else
namespace boost
#endif
{
template <typename T>
struct hash<CDT::V2d<T> >
{
    size_t operator()(const CDT::V2d<T>& xy) const
    {
#ifdef CDT_CXX11_IS_SUPPORTED
        typedef std::hash<T> Hasher;
#else
        typedef boost::hash<T> Hasher;
#endif
        return Hasher()(xy.x) ^ Hasher()(xy.y);
    }
};
} // namespace std

namespace CDT
{

template <
    typename T,
    typename TVertexIter,
    typename TGetVertexCoordX,
    typename TGetVertexCoordY>
DuplicatesInfo FindDuplicates(
    TVertexIter first,
    TVertexIter last,
    TGetVertexCoordX getX,
    TGetVertexCoordY getY)
{
    typedef unordered_map<V2d<T>, std::size_t> PosToIndex;
    PosToIndex uniqueVerts;
    const std::size_t verticesSize = std::distance(first, last);
    DuplicatesInfo di = {
        std::vector<std::size_t>(verticesSize), std::vector<std::size_t>()};
    for(std::size_t iIn = 0, iOut = iIn; iIn < verticesSize; ++iIn, ++first)
    {
        typename PosToIndex::const_iterator it;
        bool isUnique;
        tie(it, isUnique) = uniqueVerts.insert(
            std::make_pair(V2d<T>::make(getX(*first), getY(*first)), iOut));
        if(isUnique)
        {
            di.mapping[iIn] = iOut++;
            continue;
        }
        di.mapping[iIn] = it->second; // found a duplicate
        di.duplicates.push_back(iIn);
    }
    return di;
}

template <typename TVertex, typename TAllocator>
void RemoveDuplicates(
    std::vector<TVertex, TAllocator>& vertices,
    const std::vector<std::size_t>& duplicates)
{
    vertices.erase(
        remove_at(
            vertices.begin(),
            vertices.end(),
            duplicates.begin(),
            duplicates.end()),
        vertices.end());
}

template <
    typename T,
    typename TVertex,
    typename TGetVertexCoordX,
    typename TGetVertexCoordY,
    typename TVertexAllocator,
    typename TEdgeAllocator>
DuplicatesInfo RemoveDuplicatesAndRemapEdges(
    std::vector<TVertex, TVertexAllocator>& vertices,
    std::vector<Edge, TEdgeAllocator>& edges,
    TGetVertexCoordX getX,
    TGetVertexCoordY getY)
{
    const DuplicatesInfo di =
        FindDuplicates<T>(vertices.begin(), vertices.end(), getX, getY);
    RemoveDuplicates(vertices, di.duplicates);
    RemapEdges(edges, di.mapping);
    return di;
}

} // namespace CDT

#ifndef CDT_USE_AS_COMPILED_LIBRARY
#include "CDT.hpp"
#endif

#endif // header-guard
