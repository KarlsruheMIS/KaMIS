 /******************************************************************************
 * modified.h
 *
 *****************************************************************************/

#ifndef MODIFIED_H
#define MODIFIED_H

#include <vector>
#include <cassert>

class branch_and_reduce_algorithm;

class modified {

public:
    int add;
    std::vector<int> removed;
    std::vector<int> vs;
    std::vector<std::vector<int>> oldAdj;
    branch_and_reduce_algorithm *pAlg;

public:
    modified(int const add, std::vector<int> &removed, std::vector<int> &vs, std::vector<std::vector<int>> &newAdj, branch_and_reduce_algorithm *_pAlg);

    modified(std::vector<int> &removed, std::vector<int> &vs, branch_and_reduce_algorithm *_pAlg);

    virtual ~modified() {};

    void restore();

    virtual void reverse(std::vector<int> &x) = 0;
};

class fold : public modified
{

public:
    fold(int const add, std::vector<int> &removed, std::vector<int> &vs, std::vector<std::vector<int>> &newAdj, branch_and_reduce_algorithm *_pAlg)
    : modified(add, removed, vs, newAdj, _pAlg)
    { }

    fold(std::vector<int> &removed, std::vector<int> &vs, branch_and_reduce_algorithm *_pAlg)
    : modified(removed, vs, _pAlg)
    { }

    virtual ~fold() {}

    virtual void reverse(std::vector<int> &x) {
        int k = removed.size() / 2;
        if (x[vs[0]] == 0) {
            for (int i = 0; i < k; i++) x[removed[i]] = 1;
            for (int i = 0; i < k; i++) x[removed[k + i]] = 0;
        } else if (x[vs[0]] == 1) {
            for (int i = 0; i < k; i++) x[removed[i]] = 0;
            for (int i = 0; i < k; i++) x[removed[k + i]] = 1;
        }
    }
};

class alternative : public modified {

public:
    int k;

    alternative(int const add, std::vector<int> &removed, std::vector<int> &vs, std::vector<std::vector<int>> &newAdj, branch_and_reduce_algorithm *_pAlg, int k)
    : modified(add, removed, vs, newAdj, _pAlg)
    {
        this->k = k;
    }

    alternative(std::vector<int> &removed, std::vector<int> &vs, branch_and_reduce_algorithm *_pAlg, int k)
    : modified(removed, vs, _pAlg)
    {
        this->k = k;
    }

    virtual ~alternative() {}

    void reverse(std::vector<int> &x) {
        bool A0 = false, A1 = true;
        bool B0 = false, B1 = true;
        for (int i = 0; i < k; i++) {
            if (x[vs[i]] == 0) A0 = true;
            if (x[vs[i]] != 1) A1 = false;
        }
        for (int i = k; i < static_cast<int>(vs.size()); i++) {
            if (x[vs[i]] == 0) B0 = true;
            if (x[vs[i]] != 1) B1 = false;
        }
        if (A1 || B0) {
            for (int i = 0; i < static_cast<int>(removed.size()/ 2); i++) x[removed[i]] = 0;
            for (int i = static_cast<int>(removed.size() / 2); i < static_cast<int>(removed.size()); i++) x[removed[i]] = 1;
        } else if (B1 || A0) {
            for (int i = 0; i < static_cast<int>(removed.size() / 2); i++) x[removed[i]] = 1;
            for (int i = static_cast<int>(removed.size() / 2); i < static_cast<int>(removed.size()); i++) x[removed[i]] = 0;
        }
    }
};

#endif // MODIFIED_H
