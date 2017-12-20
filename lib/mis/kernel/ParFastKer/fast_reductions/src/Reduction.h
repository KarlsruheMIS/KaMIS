 /******************************************************************************
 * Copyright (C) 2015-2017 Darren Strash <strash@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *****************************************************************************/

#ifndef REDUCTION_H
#define REDUCTION_H

#include <vector>
#include <utility>

enum ReductionType {ISOLATED_VERTEX, FOLDED_VERTEX, DOMINATED_VERTEX, REMOVED_VERTEX, REMOVED_VERTEX_AND_NEIGHBORS, FOLDED_TWINS};

class Reduction
{

public:
    Reduction(ReductionType const type)
    : m_iVertex(-1)
    , m_vNeighbors()
    , m_vRemovedEdges()
    , reductionType(type) {
    }

    void SetVertex(int const vertex)
    {
        m_iVertex = vertex;
    }

    void SetKeptVertex(int const vertex) {
        keptVertex = vertex;
    }

    void SetTwin(int const twin)
    {
        m_iTwin = twin;
    }

    int GetKeptVertex() const {return keptVertex; }

    int GetVertex() const { return m_iVertex; }

    int GetTwin() const { return m_iTwin; }

    void AddNeighbor(int const neighbor) {
        m_vNeighbors.push_back(neighbor);
    }

    std::vector<int> const &GetNeighbors() const
    {
        return m_vNeighbors;
    }

    void AddRemovedEdge(int const v1, int const v2)
    {
        m_vRemovedEdges.push_back(std::make_pair(v1, v2));
    }

    std::vector<std::pair<int,int>> const &GetRemovedEdges() const
    {
        return m_vRemovedEdges;
    }

    ReductionType GetType() const { return reductionType; }

private:
    int m_iVertex;
    int m_iTwin;
    std::vector<int>                m_vNeighbors;
    std::vector<std::pair<int,int>> m_vRemovedEdges;
    ReductionType                   reductionType;
    int keptVertex;
};

#endif // REDUCTION_H
