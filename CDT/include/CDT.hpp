/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

/**
 * @file
 * Public API - implementation
 */

#include "CDT.h"

#include "CDTUtils.h"
#include <algorithm>
#include <deque>
#include <stdexcept>

namespace CDT
{

typedef std::deque<TriInd> TriDeque;


template <typename T>
DuplicatesInfo RemoveDuplicates(std::vector<V2d<T> >& vertices)
{
    const DuplicatesInfo di = FindDuplicates<T>(
        vertices.begin(), vertices.end(), getX_V2d<T>, getY_V2d<T>);
    RemoveDuplicates(vertices, di.duplicates);
    return di;
}

CDT_INLINE_IF_HEADER_ONLY void
RemapEdges(std::vector<Edge>& edges, const std::vector<std::size_t>& mapping)
{
    for(std::vector<Edge>::iterator it = edges.begin(); it != edges.end(); ++it)
    {
        *it = Edge(mapping[it->v1()], mapping[it->v2()]); // remap
    }
}

template <typename T>
DuplicatesInfo RemoveDuplicatesAndRemapEdges(
    std::vector<V2d<T> >& vertices,
    std::vector<Edge>& edges)
{
    return RemoveDuplicatesAndRemapEdges<T>(
        vertices, edges, getX_V2d<T>, getY_V2d<T>);
}

CDT_INLINE_IF_HEADER_ONLY
unordered_map<TriInd, LayerDepth> PeelLayer(
    std::stack<TriInd> seeds,
    const TriangleVec& triangles,
    const EdgeUSet& fixedEdges,
    const unordered_map<Edge, BoundaryOverlapCount>& overlapCount,
    const LayerDepth layerDepth,
    std::vector<LayerDepth>& triDepths)
{
    unordered_map<TriInd, LayerDepth> behindBoundary;
    while(!seeds.empty())
    {
        const TriInd iT = seeds.top();
        seeds.pop();
        triDepths[iT] = layerDepth;
        behindBoundary.erase(iT);
        const Triangle& t = triangles[iT];
        for(Index i(0); i < Index(3); ++i)
        {
            const Edge opEdge(t.vertices[ccw(i)], t.vertices[cw(i)]);
            const TriInd iN = t.neighbors[opoNbr(i)];
            if(iN == noNeighbor || triDepths[iN] <= layerDepth)
                continue;
            if(fixedEdges.count(opEdge))
            {
                const unordered_map<Edge, LayerDepth>::const_iterator cit =
                    overlapCount.find(opEdge);
                const LayerDepth triDepth = cit == overlapCount.end()
                                                ? layerDepth + 1
                                                : layerDepth + cit->second + 1;
                behindBoundary[iN] = triDepth;
                continue;
            }
            seeds.push(iN);
        }
    }
    return behindBoundary;
}

CDT_INLINE_IF_HEADER_ONLY
TriIndUSet PeelLayer(
    std::stack<TriInd> seeds,
    const TriangleVec& triangles,
    const EdgeUSet& fixedEdges,
    const LayerDepth layerDepth,
    std::vector<LayerDepth>& triDepths)
{
    TriIndUSet behindBoundary;
    while(!seeds.empty())
    {
        const TriInd iT = seeds.top();
        seeds.pop();
        triDepths[iT] = layerDepth;
        behindBoundary.erase(iT);
        const Triangle& t = triangles[iT];
        for(Index i(0); i < Index(3); ++i)
        {
            const Edge opEdge(t.vertices[ccw(i)], t.vertices[cw(i)]);
            const TriInd iN = t.neighbors[opoNbr(i)];
            if(iN == noNeighbor || triDepths[iN] <= layerDepth)
                continue;
            if(fixedEdges.count(opEdge))
            {
                behindBoundary.insert(iN);
                continue;
            }
            seeds.push(iN);
        }
    }
    return behindBoundary;
}

std::vector<LayerDepth> CalculateTriangleDepths(
    const TriInd seed,
    const TriangleVec& triangles,
    const EdgeUSet& fixedEdges,
    const unordered_map<Edge, BoundaryOverlapCount>& overlapCount)
{
    std::vector<LayerDepth> triDepths(
        triangles.size(), std::numeric_limits<LayerDepth>::max());
    std::stack<TriInd> seeds(TriDeque(1, seed));
    LayerDepth layerDepth = 0;
    LayerDepth deepestSeedDepth = 0;

    unordered_map<LayerDepth, TriIndUSet> seedsByDepth;
    do
    {
        const unordered_map<TriInd, LayerDepth>& newSeeds = PeelLayer(
            seeds, triangles, fixedEdges, overlapCount, layerDepth, triDepths);

        seedsByDepth.erase(layerDepth);
        typedef unordered_map<TriInd, LayerDepth>::const_iterator Iter;
        for(Iter it = newSeeds.begin(); it != newSeeds.end(); ++it)
        {
            deepestSeedDepth = std::max(deepestSeedDepth, it->second);
            seedsByDepth[it->second].insert(it->first);
        }
        const TriIndUSet& nextLayerSeeds = seedsByDepth[layerDepth + 1];
        seeds = std::stack<TriInd>(
            TriDeque(nextLayerSeeds.begin(), nextLayerSeeds.end()));
        ++layerDepth;
    } while(!seeds.empty() || deepestSeedDepth > layerDepth);

    return triDepths;
}

std::vector<LayerDepth> CalculateTriangleDepths(
    const TriInd seed,
    const TriangleVec& triangles,
    const EdgeUSet& fixedEdges)
{
    std::vector<LayerDepth> triDepths(
        triangles.size(), std::numeric_limits<LayerDepth>::max());
    std::stack<TriInd> seeds(TriDeque(1, seed));
    LayerDepth layerDepth = 0;

    do
    {
        const TriIndUSet& newSeeds =
            PeelLayer(seeds, triangles, fixedEdges, layerDepth++, triDepths);
        seeds = std::stack<TriInd>(TriDeque(newSeeds.begin(), newSeeds.end()));
    } while(!seeds.empty());

    return triDepths;
}

CDT_INLINE_IF_HEADER_ONLY EdgeUSet
extractEdgesFromTriangles(const TriangleVec& triangles)
{
    EdgeUSet edges;
    typedef TriangleVec::const_iterator CIt;
    for(CIt t = triangles.begin(); t != triangles.end(); ++t)
    {
        edges.insert(Edge(VertInd(t->vertices[0]), VertInd(t->vertices[1])));
        edges.insert(Edge(VertInd(t->vertices[1]), VertInd(t->vertices[2])));
        edges.insert(Edge(VertInd(t->vertices[2]), VertInd(t->vertices[0])));
    }
    return edges;
}

} // namespace CDT
