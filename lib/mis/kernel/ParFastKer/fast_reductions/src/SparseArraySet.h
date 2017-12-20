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

#ifndef SPARSE_ARRAY_SET_H
#define SPARSE_ARRAY_SET_H

#include <set>
#include <vector>
#include <iostream>
#include <cassert>
#include <utility>

class SparseArraySet
{
public:
    SparseArraySet(size_t const size) : m_Elements(size, -1), m_iBegin(0), m_iEnd(-1), m_States(0), m_bInserted(false), m_bRemoved(false)
    {
    }

    SparseArraySet() : m_Elements(), m_iBegin(0), m_iEnd(-1), m_States(0), m_bInserted(false), m_bRemoved(false)
    {
    }

    ~SparseArraySet() {}

    void Resize(size_t const size)
    {
        assert(m_Elements.size() >= Size());
        m_Elements.resize(size, -1);
        assert(m_Elements.size() >= Size());
    }

    void InitializeFromAdjacencyArray(std::vector<std::vector<int>> const &adjacencyArray, int const vertex)
    {
        Resize(adjacencyArray[vertex].size());
        for (int const neighbor : adjacencyArray[vertex]) {
            assert(!m_bRemoved); // not allowed to insert and remove when saving states
            if (!m_States.empty()) m_bInserted = true;
            m_iEnd++;
            m_Elements[m_iEnd] = neighbor;
        }
        assert(m_Elements.size() >= Size());
    }

    void PrintSummary() const
    {
        std::cout << "Array[" << m_iBegin << ":" << m_iEnd << "] : ";
        for (int const p : *this) {
            std::cout << p << (Contains(p)? "":"<-inconsistency") << " ";
        }
        std::cout << std::endl;
    }

    bool Contains(int const x) const {
        for (int const value : *this) {
            if (value == x) return true;
        }
        return false;
    }

    // Inserts are not allowed after saving state, as it is currently not supported.
    void Insert(int const x) {
        if (Contains(x)) return;
        assert(m_Elements.size() > Size());
        assert(!m_bRemoved); // not allowed to insert and remove when saving states
        if (!m_States.empty()) m_bInserted = true;
        m_iEnd++;
        m_Elements[m_iEnd] = x;
    }

    void Remove(int const x) {
        assert(m_Elements.size() >= Size());
        for (int index = m_iBegin; index <= m_iEnd; index++) {
            if (m_Elements[index] == x) {
                if (!m_States.empty()) m_bRemoved = true;
                assert(!m_bInserted); // not allowed to insert and remove when saving states
                m_Elements[index] = m_Elements[m_iEnd];
                m_Elements[m_iEnd] = x;
                m_iEnd--;
                assert(m_Elements.size() > Size());
                return;
            }
        }
        return;
    }

    void MoveTo(int const x, SparseArraySet &other) {
        Remove(x);
        other.Insert(x);
    }

    size_t Size()  const { return m_iEnd - m_iBegin + 1; }
    bool   Empty() const { return m_iEnd < m_iBegin; }

    std::vector<int>::iterator begin() { return m_Elements.begin() + m_iBegin;   }
    std::vector<int>::iterator end()   { return m_Elements.begin() + m_iEnd + 1; }

    std::vector<int>::const_iterator begin() const { return m_Elements.begin() + m_iBegin;   }
    std::vector<int>::const_iterator end()   const { return m_Elements.begin() + m_iEnd + 1; }

    int At(size_t const index) const
    {
        return m_Elements[index];
    }

    int operator[](size_t const index) const { return At(index); }


    void SaveState()
    {
////        std::cout << "Saving State" << std::endl << std::flush;
        m_States.push_back(std::make_pair(m_iBegin, m_iEnd));
    }

    void RestoreState()
    {
////        std::cout << "Restoring State" << std::endl << std::flush;
        std::pair<int,int> const &range(m_States.back());
        m_iBegin = range.first;
        m_iEnd   = range.second;
        m_States.pop_back();
////        std::cout << "#States = " << m_States.size() << std::endl << std::flush;
        if (m_States.empty()) {
            m_bRemoved  = false;
            m_bInserted = false;
        }
    }

    void Clear()
    {
        m_iBegin = 0;
        m_iEnd = -1;
    }

    bool Equals(SparseArraySet &other) {
        for(int element: other) {
            if(!(this->Contains(element)))
                return false;
        }
        for(int element: *this) {
            if(!other.Contains(element))
                return false;
        }
        return true;
    }

    static bool Test() {
        std::cout << "SparseArraySet: ";
        SparseArraySet testSet(3);
        testSet.Insert(0);
        testSet.Insert(1);
        testSet.Insert(2);

        if (!testSet.Contains(0) || !testSet.Contains(1) || !testSet.Contains(2) || testSet.Size() != 3) {
            std::cout << "FAILED: SparseArraySet failed ContainsTest" << std::endl;
            return false;
        }

        testSet.Remove(0);
        if (testSet.Contains(0) || testSet.Size() != 2) {
            std::cout << "FAILED: SparseArraySet failed RemoveFirst Test" << std::endl;
            return false;
        }

        testSet.Insert(0);
        if (!testSet.Contains(0) || !testSet.Contains(1) || !testSet.Contains(2) || testSet.Size() != 3) {
            std::cout << "FAILED: SparseArraySet failed RemoveAdd Test" << std::endl;
            return false;
        }

        testSet.Remove(0);
        if (testSet.Contains(0) || testSet.Size() != 2) {
            std::cout << "FAILED: SparseArraySet failed RemoveLast Test" << std::endl;
            return false;
        }

        testSet.Insert(0);
        testSet.Remove(1);
        if (testSet.Contains(1) || testSet.Size() != 2) {
            std::cout << "FAILED: SparseArraySet failed RemoveMiddle Test" << std::endl;
            return false;
        }

        testSet.Insert(1);
        testSet.Remove(1);
        testSet.Remove(2);
        testSet.Remove(0);
        if (testSet.Contains(0) || testSet.Contains(1) || testSet.Contains(2) || testSet.Size() != 0 || !testSet.Empty()) {
            std::cout << "FAILED: SparseArraySet failed RemoveAll Test" << std::endl;
            return false;
        }

        testSet.Insert(0);
        testSet.Insert(1);
        testSet.Insert(2);
        if (!testSet.Contains(0) || !testSet.Contains(1) || !testSet.Contains(2) || testSet.Size() != 3) {
            std::cout << "FAILED: SparseArraySet failed InsertAll Test" << std::endl;
            return false;
        }

        testSet.SaveState();
        if (!testSet.Contains(0) || !testSet.Contains(1) || !testSet.Contains(2) || testSet.Size() != 3) {
            std::cout << "FAILED: SparseArraySet failed SaveState Test" << std::endl;
            return false;
        }

        testSet.Remove(0);
        if (testSet.Contains(0) || !testSet.Contains(1) || !testSet.Contains(2) || testSet.Size() != 2) {
            std::cout << "FAILED: SparseArraySet failed RemoveAfterSave Test" << std::endl;
            return false;
        }

        testSet.RestoreState();
        if (!testSet.Contains(0) || !testSet.Contains(1) || !testSet.Contains(2) || testSet.Size() != 3) {
            std::cout << "FAILED: SparseArraySet failed RestoreState Test" << std::endl;
            return false;
        }

        std::cout << "PASSED!" << std::endl << std::flush;
        return true;
    }

private:
    std::vector<int> m_Elements;
    int m_iBegin;
    int m_iEnd;
    std::vector<std::pair<int,int>> m_States;
    bool m_bInserted;
    bool m_bRemoved;
};

#endif // SPARSE_ARRAY_SET_H
