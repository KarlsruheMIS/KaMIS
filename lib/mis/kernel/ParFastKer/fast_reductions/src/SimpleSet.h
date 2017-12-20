 /******************************************************************************
 * Copyright (C) 2019 Demian Hespe <hespe@kit.edu>
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

#ifndef SIMPLE_SET_H
#define SIMPLE_SET_H

#include <vector>

class SimpleSet
{
public:
    SimpleSet(size_t const size, bool const init) : elements(size, init)
    {
    }

    SimpleSet() : elements()
    {
    }

    ~SimpleSet() {}

    void Resize(size_t const size)
    {
        elements.resize(size, false);
    }

    bool Contains(int const x) const {
        assert(x < elements.size());
        return elements[x];
    }

    void Insert(int const x) {
        assert(x < elements.size());
        elements[x] = true;
    }

    void Remove(int const x) {
        assert(x < elements.size());
        elements[x] = false;
    }

    size_t Size()  const {
        int count = 0;
        for(auto element: elements) {
            if (element) {
                count++;
            }
        }
        return count; 
    }

private:
    std::vector<char> elements;
};

#endif // SIMPLE_SET_H
