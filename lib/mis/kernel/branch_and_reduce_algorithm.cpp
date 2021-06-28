 
#include "branch_and_reduce_algorithm.h"
#include "mis_config.h"
#include "fast_set.h"
#include "modified.h"
#include "data_structure/graph_access.h"
#include "data_structure/priority_queues/maxNodeHeap.h"

#include <vector>
#include <set>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <algorithm>  // sort()
#include <deque>
#include <chrono>

////#define debug(x) {fprintf(2, x);}

using namespace std;

int  branch_and_reduce_algorithm::REDUCTION   = 3;
int  branch_and_reduce_algorithm::LOWER_BOUND = 4;
int  branch_and_reduce_algorithm::BRANCHING   = 2;
bool branch_and_reduce_algorithm::outputLP    = false;
long branch_and_reduce_algorithm::nBranchings = 0;
int  branch_and_reduce_algorithm::debug       = 0;

branch_and_reduce_algorithm::branch_and_reduce_algorithm(vector<vector<int>> &_adj, int const _N)
: adj() 
, n(_adj.size())
, used(n*2)
{
////    srand(4327897);
    SHRINK = 0.5;
    depth = 0;
    maxDepth = 10;
    rootDepth = -1; // invalid value

    n = _adj.size();
    adj.swap(_adj);

    N = _N;
    opt = n;
    y.resize(N, 0);
    for (int i = 0; i < n; i++) y[i] = 1;
    for (int i = n; i < N; i++) y[i] = 2;
    crt = 0;
    x.resize(N, 0);
    for (int i = 0; i < n; i++) x[i] = -1;
    for (int i = n; i < N; i++) x[i] = 2;
    rn = n;
    in.resize(n, -1);
    out.resize(n, -1);
    lb = -1; // invalid value

    vRestore.resize(n, 0);

    que.resize(n * 2, 0);
    level.resize(n * 2, 0);
////    cout << "level.size=" << level.size() << endl << flush;
    iter.resize(n * 2, 0);

    modTmp.resize(n, 0);

    modifiedN = 0;
    modifieds.resize(N, shared_ptr<modified>());
////    packing.reserve(N);
}

int branch_and_reduce_algorithm::deg(int v) {
    assert(x[v] < 0);
    int deg = 0;
    for (int u : adj[v]) if (x[u] < 0) deg++;
    return deg;
}

void branch_and_reduce_algorithm::set(int v, int a)
{
    assert(x[v] < 0);
    crt += a;
    x[v] = a;
    vRestore[--rn] = v;
    if (a == 0) {
        for (int u : adj[v]) if (x[u] < 0) {
            x[u] = 1;
            crt++;
            vRestore[--rn] = u;
        }
    }
}

// methods that modify the graph

void branch_and_reduce_algorithm::compute_fold(vector<int> const &S, vector<int> const &NS) {
    assert(NS.size() == S.size() + 1);
    vector<int> removed(S.size() * 2);
    for (unsigned int i = 0; i < S.size(); i++) removed[i] = S[i];
    for (unsigned int i = 0; i < S.size(); i++) removed[S.size() + i] = NS[1 + i];
    int s = NS[0];
    used.clear();
    for (int v : S) used.add(v);
    vector<int> &tmp = modTmp;
    int p = 0;
    for (int v : NS) {
        assert(!used.get(v));
        for (int u : adj[v]) if (x[u] < 0 && used.add(u)) {
            tmp[p++] = u;
        }
    }
    vector<vector<int>> newAdj(p + 1);
    {
        vector<int> copyOfTmp(tmp.begin(), tmp.begin() + p);
        newAdj[0].swap(copyOfTmp);
    }
    std::sort(newAdj[0].begin(), newAdj[0].end());
    vector<int> vs(p + 1);
    vs[0] = s;
    used.clear();
    for (int v : S) used.add(v);
    for (int v : NS) used.add(v);
    for (unsigned int i = 0; i < newAdj[0].size(); i++) {
        int v = newAdj[0][i];
        p = 0;
        bool add = false;
        for (int u : adj[v]) if (x[u] < 0 && !used.get(u)) {
            if (!add && s < u) {
                tmp[p++] = s;
                add = true;
            }
            tmp[p++] = u;
        }
        if (!add) tmp[p++] = s;
        vs[1 + i] = v;

        {
            vector<int> copyOfTmp(tmp.begin(), tmp.begin() + p);
            newAdj[i+1].swap(copyOfTmp);
        }
    }
    modifieds[modifiedN++] = make_shared<fold>(fold(S.size(), removed, vs, newAdj, this));
////    cout << __LINE__ << ", " << this << ", " << depth << ": Setting modifieds[" << modifiedN-1 << "]=" << modifieds[modifiedN-1] << endl << flush;
}

void branch_and_reduce_algorithm::compute_alternative(vector<int> const &A, vector<int> const &B) {
    assert(A.size() == B.size());
    used.clear();
    for (int b : B) for (int u : adj[b]) if (x[u] < 0) used.add(u);
    for (int a : A) for (int u : adj[a]) if (x[u] < 0 && used.get(u)) set(u, 1);
    NodeID p = 0, q = 0;
    vector<int> &tmp = modTmp;
    used.clear();
    for (int b : B) used.add(b);
    for (int a : A) for (int u : adj[a]) if (x[u] < 0 && used.add(u)) tmp[p++] = u;
    vector<int> A2(tmp.begin(), tmp.begin() + p);
    std::sort(A2.begin(), A2.end());
    p = 0;
    used.clear();
    for (int a : A) used.add(a);
    for (int b : B) for (int u : adj[b]) if (x[u] < 0 && used.add(u)) tmp[p++] = u;
    vector<int> B2(tmp.begin(), tmp.begin() + p);
    std::sort(B2.begin(), B2.end());
    vector<int> removed(A.size() + B.size());
    for (unsigned int i = 0; i < A.size(); i++) removed[i] = A[i];
    for (unsigned int i = 0; i < B.size(); i++) removed[A.size() + i] = B[i];
    vector<int> vs(A2.size() + B2.size());
    for (unsigned int i = 0; i < A2.size(); i++) vs[i] = A2[i];
    for (unsigned int i = 0; i < B2.size(); i++) vs[A2.size() + i] = B2[i];
    vector<vector<int>> newAdj(vs.size());
    used.clear();
    for (int a : A) used.add(a);
    for (int b : B) used.add(b);
    for (unsigned int i = 0; i < vs.size(); i++) {
        unsigned int v = (i < A2.size()) ? A2[i] : B2[i - A2.size()];
        vector<int> const &C = (i < A2.size()) ? B2 : A2;
        p = q = 0;
        for (int u : adj[v]) if (x[u] < 0 && !used.get(u)) {
            while (q < C.size() && C[q] <= u) {
                if (used.get(C[q])) q++;
                else tmp[p++] = C[q++];
            }
            if (p == 0 || tmp[p - 1] != u) tmp[p++] = u;
        }
        while (q < C.size()) {
            if (used.get(C[q])) q++;
            else tmp[p++] = C[q++];
        }
        {
            vector<int> copyOfTmp(tmp.begin(), tmp.begin() + p);
            newAdj[i].swap(copyOfTmp);
        }
    }
    modifieds[modifiedN++] = make_shared<alternative>(alternative(removed.size()/ 2, removed, vs, newAdj, this, A2.size()));
////    cout << __LINE__ << ", " << this << ", " << depth << ": Setting modifieds[" << modifiedN-1 << "]=" << modifieds[modifiedN-1] << endl << flush;
}

void branch_and_reduce_algorithm::restore(int n) {
    while (rn < n) {
        int v = vRestore[rn];
        if (v >= 0) {
            crt -= x[v];
            x[v] = -1;
            rn++;
        } else {
            modifieds[--modifiedN]->restore();
////            cout << __LINE__ << ", " << this << ", " << depth << ": Deleting " << modifieds[modifiedN] << endl << flush;
            modifieds[modifiedN] = shared_ptr<modified>();
////            delete modifieds[modifiedN];
////            modifieds[modifiedN] = nullptr;
        }
    }
}

void branch_and_reduce_algorithm::reverse() {
    for (int i = modifiedN - 1; i >= 0; i--) {
////        cout << __LINE__ << ", " << this << ": modifieds[" << i << "] is " << modifieds[i] << endl << flush;
        modifieds[i]->reverse(y);
    }
}


// lower bounds

int branch_and_reduce_algorithm::lpLowerBound() {
    ////    try (Stat stat = new Stat("lb_LP")) {
    return crt + (rn + 1) / 2;
    ////    }
}

int branch_and_reduce_algorithm::cycleLowerBound() {
    ////    try (Stat stat = new Stat("lb_cycle")) {
    int lb = crt;
    vector<int> &id = iter;
    for (int i = 0; i < n; i++) id[i] = -1;
    vector<int> &pos = que;
    vector<int> &S  = level;
    vector<int> &S2 = modTmp;
    for (int i = 0; i < n; i++) if (x[i] < 0 && id[i] < 0) {
        int v = i;
        int size = 0;
        do {
            assert(id[v] < 0);
            id[v] = i;
            v = out[v];
            pos[v] = size;
            S[size++] = v;
        } while (v != i);
        bool clique = true;
        for (int j = 0; j < size; j++) {
            v = S[j];
            int num = 0;
            for (int u : adj[v]) if (x[u] < 0 && id[u] == id[v]) num++;
            if (num != size - 1) {
                clique = false;
                break;
            }
        }
        if (clique) {
            lb += size - 1;
        } else {
            while (size >= 6) {
                int minSize = size, s = 0, t = size;
                for (int j = 0; j < size; j++) {
                    used.clear();
                    v = S[j];
                    for (int u : adj[v]) if (x[u] < 0 && id[u] == id[v]) {
                        used.add(u);
                    }
                    v = S[(j + 1) % size];
                    for (int u : adj[v]) if (x[u] < 0 && id[u] == id[v]) {
                        if (used.get(S[(pos[u] + 1) % size])) {
                            int size2 = (pos[u] - j + size) % size;
                            if (minSize > size2 && size2 % 2 != 0) {
                                minSize = size2;
                                s = (j + 1) % size;
                                t = (pos[u] + 1) % size;
                            }
                        }
                    }
                }
                if (minSize == size) break;
                int p = 0;
                for (int j = t; j != s; j = (j + 1) % size) {
                    S2[p++] = S[j];
                }
                for (int j = s; j != t; j = (j + 1) % size) {
                    id[S[j]] = n;
                }
                
////                cout << "Before level.size=" << level.size() << endl << flush;
////                vector<int> &S3 = S; S = S2; S2 = S3;
////                cout << "After  level.size=" << level.size() << endl << flush;
                S.swap(S2);

                size -= minSize;
                assert(size == p);
                assert(minSize > 1);
                lb += (minSize + 1) / 2;
                for (int j = 0; j < size; j++) pos[S[j]] = j;
            }
            assert(size > 1);
            lb += (size + 1) / 2;
        }
    }

    if (static_cast<int>(level.size()) != n*2) {
        level.swap(modTmp);
    }
    return lb;
    ////    }
}

int branch_and_reduce_algorithm::cliqueLowerBound() {
    ////    try (Stat stat = new Stat("lb_clique")) {
    int need = crt;
    vector<long long> ls(rn, 0);
    int k = 0;
    for (int i = 0; i < n; i++) if (x[i] < 0) ls[k++] = ((long long)deg(i)) << 32 | i;
    std::sort(ls.begin(), ls.end());
    vector<int> &clique = que;
    vector<int> &size = level;
    vector<int> &tmp = iter;
    used.clear();
    for (int i = 0; i < rn; i++) {
        int v = (int)ls[i];
        int to = v, max = 0;
        for (int u : adj[v]) if (x[u] < 0 && used.get(u)) tmp[clique[u]] = 0;
        for (int u : adj[v]) if (x[u] < 0 && used.get(u)) {
            int c = clique[u];
            tmp[c]++;
            if (tmp[c] == size[c] && max < size[c]) {
                to = c;
                max = size[c];
            }
        }
        clique[v] = to;
        if (to != v) {
            size[to]++;
            need++;
        } else {
            size[v] = 1;
        }
        used.add(v);
    }
    return need;
    ////}
}

int branch_and_reduce_algorithm::lowerBound() {
    int type = 0, tmp;
    if (lb < crt) {
        lb = crt;
        type = 1;
    }
    if (LOWER_BOUND == 1 || LOWER_BOUND == 4) {
        tmp = cliqueLowerBound();
        if (lb < tmp) {
            lb = tmp;
            type = 4;
        }
    }
    if (LOWER_BOUND == 2 || LOWER_BOUND == 4) {
        tmp = lpLowerBound();
        if (lb < tmp) {
            lb = tmp;
            type = 2;
        }
    }
    if (LOWER_BOUND == 3 || LOWER_BOUND == 4) {
        tmp = cycleLowerBound();
        if (lb < tmp) {
            lb = tmp;
            type = 3;
        }
    }
    if (debug >= 2 && depth <= maxDepth) fprintf(stderr, "%slb: %d (%d), %d\n", debugString().c_str(), lb, type, opt);
    return lb;
}

// helper for lpReduction
bool branch_and_reduce_algorithm::dinicDFS(int v) {
    while (iter[v] >= 0) {
        int u = adj[v][iter[v]--], w = in[u];
        if (x[u] >= 0) continue;
        if (w < 0 || (level[v] < level[w] && iter[w] >= 0 && dinicDFS(w))) {
            in[u] = v;
            out[v] = u;
            return true;
        }
    }
    return false;
}

// helper for lpReduction
void branch_and_reduce_algorithm::updateLP() {
    ////    try (Stat stat = new Stat("updateLP")) {
#if 1
    for (int v = 0; v < n; v++) if (out[v] >= 0 && ((x[v] < 0) ^ (x[out[v]] < 0))) {
        in[out[v]] = -1;
        out[v] = -1;
    }
    for (;;) {
        used.clear();
        int qs = 0, qt = 0;
        for (int v = 0; v < n; v++) if (x[v] < 0 && out[v] < 0) {
            level[v] = 0;
            used.add(v);
            que[qt++] = v;
        }
        bool ok = false;
        while (qs < qt) {
            int v = que[qs++];
            iter[v] = adj[v].size() - 1;
            for (int u : adj[v]) if (x[u] < 0 && used.add(n + u)) {
                int w = in[u];
                if (w < 0) ok = true;
                else {
                    level[w] = level[v] + 1;
                    used.add(w);
                    que[qt++] = w;
                }
            }
        }
        if (!ok) break;
        for (int v = n - 1; v >= 0; v--) if (x[v] < 0 && out[v] < 0) {
            dinicDFS(v);
        }
    }

    ////    if (debug >= 3 && depth <= maxDepth) fprintf(stderr, "%sUpdateLP", debugString().c_str());
#endif // 0
    ////    }
}

// reductions

bool branch_and_reduce_algorithm::lpReduction() {
    int oldn = rn;
    ////    if (debug >= 3 && depth <= maxDepth) fprintf(stderr, "%sLP:start\n", debugString().c_str());
    updateLP();

////    if (modifiedN >= 5)
////        cout << __LINE__ << ", " << this << ", " << depth << ": Setting modifieds[" << 5 << "]=" << modifieds[5] << endl << flush;
    ////    if (debug >= 3 && depth <= maxDepth) fprintf(stderr, "%sLP:afterUpdate\n", debugString().c_str());
#if 1
    ////    try (Stat stat = new Stat("reduce_LP")) {
    for (int v = 0; v < n; v++) {
        if (x[v] < 0 && used.get(v) && !used.get(n + v)) set(v, 0);
    }
    ////    if (debug >= 3 && depth <= maxDepth) fprintf(stderr, "%sLP:afterSets\n", debugString().c_str());
    used.clear();
    ////    if (debug >= 3 && depth <= maxDepth) fprintf(stderr, "%sLP:afterClear\n", debugString().c_str());
    int p = 0;
    iter.assign(iter.size(), 0);
    ////    if (debug >= 3 && depth <= maxDepth) fprintf(stderr, "%sLP:startFor\n", debugString().c_str());

////    if (modifiedN >= 5)
////        cout << __LINE__ << ", " << this << ", " << depth << ": Setting modifieds[" << 5 << "]=" << modifieds[5] << endl << flush;
    for (int s = 0; s < n; s++) if (x[s] < 0 && used.add(s)) {
        ////        cout << s << "/" << n << ":" << endl;
        int qt = 0;
        que[qt] = s;
        while (qt >= 0) {
            int v = que[qt], u = -1;
////            if (modifiedN >= 5)
////                cout << __LINE__ << ", " << this << ", " << depth << ", s=" << s << ": Setting modifieds[" << 5 << "]=" << modifieds[5] << endl << flush;
            if (v < n) {
                while (iter[v] < static_cast<int>(adj[v].size())) {
                    u = n + adj[v][iter[v]++];
                    ////                    cout << ",u=" << u;
                    if (x[u - n] < 0 && used.add(u)) {
                        ////                        cout << ",break" << endl;
                        break;
                    }
                    u = -1;
                }
////                if (modifiedN >= 5)
////                    cout << __LINE__ << ", " << this << ", " << depth << ", s=" << s << ": Setting modifieds[" << 5 << "]=" << modifieds[5] << endl << flush;
            } else if (used.add(in[v - n])) {
                u = in[v - n];
////                if (modifiedN >= 5)
////                    cout << __LINE__ << ", " << this << ", " << depth << ", s=" << s << ": Setting modifieds[" << 5 << "]=" << modifieds[5] << endl << flush;
            }
            if (u >= 0) {
                que[++qt] = u;
////                if (modifiedN >= 5)
////                    cout << __LINE__ << ", " << this << ", " << depth << ", s=" << s << ": Setting modifieds[" << 5 << "]=" << modifieds[5] << endl << flush;
                ////                cout << ", q[" << qt << "]=" << u << endl;
            } else {
////                if (modifiedN >= 5)
////                    cout << __LINE__ << ", " << this << ", " << depth << ", s=" << s << ": Setting modifieds[" << 5 << "]=" << modifieds[5] << endl << flush;
////                cout << ", level[" << p << "/" << level.size() << "]=" << v << " qt=" << qt << endl;
                level[p++] = v;
                qt--;
////                if (modifiedN >= 5)
////                    cout << __LINE__ << ", " << this << ", " << depth << ", s=" << s << ": Setting modifieds[" << 5 << "]=" << modifieds[5] << endl << flush;
////                cout << ", level[" << p-1 << "/" << level.size() << "]=" << level[p-1] << ", qt=" << qt << endl;
            }
        }

////        if (modifiedN >= 5)
////                cout << __LINE__ << ", " << this << ", " << depth << ", s=" << s << ": Setting modifieds[" << 5 << "]=" << modifieds[5] << endl << flush;
    }

////        if (modifiedN >= 5)
////                cout << __LINE__ << ", " << this << ", " << depth << ": Setting modifieds[" << 5 << "]=" << modifieds[5] << endl << flush;
////    if (debug >= 3 && depth <= maxDepth) fprintf(stderr, "%sLP:firsthalf\n", debugString().c_str());
    used.clear();
    for (int i = p - 1; i >= 0; i--) if (used.add(level[i])) {
        int v = level[i];
        int qs = 0, qt = 0;
        que[qt++] = v;
        bool ok = true;
        while (qs < qt) {
            v = que[qs++];
            if (used.get(v >= n ? (v - n) : (v + n))) ok = false;
            if (v >= n) {
                for (int u : adj[v - n]) if (x[u] < 0 && used.add(u)) {
                    que[qt++] = u;
                }
            } else if (used.add(n + out[v])) {
                que[qt++] = n + out[v];
            }
        }
        ok = false;
        if (ok) {
            for (int j = 0; j < qt; j++) {
                v = que[j];
                if (v >= n) set(v - n, 0);
            }
        }
    }

////        if (modifiedN >= 5)
////                cout << __LINE__ << ", " << this << ", " << depth << ": Setting modifieds[" << 5 << "]=" << modifieds[5] << endl << flush;
    ////    }
#endif // 0
    if (debug >= 3 && depth <= maxDepth && oldn != rn) fprintf(stderr, "%sLP: %d -> %d\n", debugString().c_str(), oldn, rn);
////    else if (debug >= 3 && depth <= maxDepth) fprintf(stderr, "%sLP: %d -> %d\n", debugString().c_str(), oldn, rn);
    ////    if (oldn != rn) Stat.count("reduceN_LP", oldn - rn);
    return oldn != rn;
}


bool branch_and_reduce_algorithm::deg1Reduction() {
    ////    try (Stat stat = new Stat("reduce_deg1")) {
    int oldn = rn;
    vector<int> &deg = iter;
    int qt = 0;
    used.clear();
    for (int v = 0; v < n; v++) if (x[v] < 0) {
        deg[v] = n == rn ? adj[v].size(): this->deg(v);
        if (deg[v] <= 1) {
            que[qt++] = v;
            used.add(v);
        }
    }
    while (qt > 0) {
        int v = que[--qt];
        if (x[v] >= 0) continue;
        assert(deg[v] <= 1);
        for (int u : adj[v]) if (x[u] < 0) {
            for (int w : adj[u]) if (x[w] < 0) {
                deg[w]--;
                if (deg[w] <= 1 && used.add(w)) que[qt++] = w;
            }
        }
////        cout << "Setting " << v << " to zero" << endl << flush;
        set(v, 0);
    }
    if (debug >= 3 && depth <= maxDepth && oldn != rn) fprintf(stderr, "%sdeg1: %d -> %d\n", debugString().c_str(), oldn, rn);
////    else if (debug >= 3 && depth <= maxDepth) fprintf(stderr, "%sdeg1: %d -> %d\n", debugString().c_str(), oldn, rn);
    ////    if (oldn != rn) Stat.count("reduceN_deg1", oldn - rn);
    return oldn != rn;
    ////    }
}

bool branch_and_reduce_algorithm::dominateReduction() {
    ////    try (Stat stat = new Stat("reduce_dominate")) {
    int oldn = rn;
#if 1
    for (int v = 0; v < n; v++) if (x[v] < 0) {
        used.clear();
        used.add(v);
        for (int u : adj[v]) if (x[u] < 0) used.add(u);
        for (int u : adj[v]) if (x[u] < 0) {
            for (int w : adj[u]) {
                if (x[w] < 0 && !used.get(w)) goto loop;
            }
            set(v, 1);
            break;
 loop:       ;
        }
    }
    if (debug >= 3 && depth <= maxDepth && oldn != rn) fprintf(stderr, "%sdominate: %d -> %d\n", debugString().c_str(), oldn, rn);
    ////        if (oldn != rn) Stat.count("reduceN_dominate", oldn - rn);
#endif // 0
    return oldn != rn;
    ////    }
}

bool branch_and_reduce_algorithm::fold2Reduction() {
    ////    try (Stat stat = new Stat("reduce_fold2")) {
    int oldn = rn;
    vector<int> &tmp = level;
    for (int v = 0; v < n; v++) if (x[v] < 0) {
        int p = 0;
        for (int u : adj[v]) if (x[u] < 0) {
            tmp[p++] = u;
            if (p > 2) goto loop;
        }
        if (p < 2) continue;
        for (int u : adj[tmp[0]]) if (u == tmp[1]) {
            set(v, 0);
            goto loop;
        }
        {
        vector<int> copyOfTmp(tmp.begin(), tmp.begin() + 2);
        compute_fold(vector<int>{v}, copyOfTmp);
        }
loop:   ;
    }
    if (debug >= 3 && depth <= maxDepth && oldn != rn) fprintf(stderr, "%sfold2: %d -> %d\n", debugString().c_str(), oldn, rn);
    ////       if (oldn != rn) Stat.count("reduceN_fold2", oldn - rn);
    return oldn != rn;
    ////    }
}

bool branch_and_reduce_algorithm::twinReduction() {
    ////    try (Stat stat = new Stat("reduce_twin")) {
    int oldn = rn;
    vector<int> &vUsed = iter;
    int uid = 0;
    vector<int> NS(3, 0);
    for (int i = 0; i < n; i++) vUsed[i] = 0;
    for (int v = 0; v < n; v++) if (x[v] < 0 && deg(v) == 3) {
        int p = 0;
        for (int u : adj[v]) if (x[u] < 0) {
            NS[p++] = u;
            uid++;
            for (int w : adj[u]) if (x[w] < 0 && w != v) {
                if (p == 1) vUsed[w] = uid;
                else if (vUsed[w] == uid - 1) {
                    vUsed[w]++;
                    if (p == 3 && deg(w) == 3) {
                        uid++;
                        for (int z : NS) vUsed[z] = uid;
                        bool ind = true;
                        for (int z : NS) for (int a : adj[z]) if (x[a] < 0 && vUsed[a] == uid) ind = false;
                        if (ind) {
                            compute_fold(vector<int>{v, w}, NS);
                        } else {
                            set(v, 0);
                            set(w, 0);
                        }
                        goto loop;
                    }
                }
            }
        }
loop:      ;
    }
    if (debug >= 3 && depth <= maxDepth && oldn != rn) fprintf(stderr, "%stwin: %d -> %d\n", debugString().c_str(), oldn, rn);
    ////       if (oldn != rn) Stat.count("reduceN_twin", oldn - rn);
    return oldn != rn;
    ////    }
}

bool branch_and_reduce_algorithm::funnelReduction() {
    ////    try (Stat stat = new Stat("reduce_alternative")) {
    int oldn = rn;
    for (int v = 0; v < n; v++) if (x[v] < 0) {
        used.clear();
        vector<int> &tmp = level;
        int p = 0;
        for (int u : adj[v]) if (x[u] < 0 && used.add(u)) {
            tmp[p++] = u;
        }
        if (p <= 1) {
            set(v, 0);
            continue;
        }
        int u1 = -1;
        for (int i = 0; i < p; i++) {
            int d = 0;
            for (int u : adj[tmp[i]]) if (x[u] < 0 && used.get(u)) d++;
            if (d + 1 < p) {
                u1 = tmp[i];
                break;
            }
        }
        if (u1 < 0) {
            set(v, 0);
            continue;
        } else {
            vector<int> &id = iter;
            for (int i = 0; i < p; i++) id[tmp[i]] = -1;
            for (int u : adj[u1]) if (x[u] < 0) id[u] = 0;
            int u2 = -1;
            for (int i = 0; i < p; i++) if (tmp[i] != u1 && id[tmp[i]] < 0) {
                u2 = tmp[i];
                break;
            }
            assert(u2 >= 0);
            used.remove(u1);
            used.remove(u2);
            int d1 = 0, d2 = 0;
            for (int w : adj[u1]) if (x[w] < 0 && used.get(w)) d1++;
            for (int w : adj[u2]) if (x[w] < 0 && used.get(w)) d2++;
            if (d1 < p - 2 && d2 < p - 2) continue;
            for (int i = 0; i < p; i++) {
                int u = tmp[i];
                if (u == u1 || u == u2) continue;
                int d = 0;
                for (int w : adj[u]) if (x[w] < 0 && used.get(w)) d++;
                if (d < p - 3) {
                    goto loop;
                }
            }
            int u = (d1 == p - 2) ? u2 : u1;
            vector<int> const v1{v};
            vector<int> const v2{u};
            compute_alternative(v1, v2);
        }
loop:      ;
    }
    if (debug >= 3 && depth <= maxDepth && oldn != rn) fprintf(stderr, "%sfunnel: %d -> %d\n", debugString().c_str(), oldn, rn);
    ////       if (oldn != rn) Stat.count("reduceN_funnel", oldn - rn);
    return oldn != rn;
    ////    }
}

bool branch_and_reduce_algorithm::deskReduction() {
    ////    try (Stat stat = new Stat("reduce_desk")) {
    int oldn = rn;
#if 1
    vector<int> &tmp = level;
    vector<int> &nv = iter;
    for (int i = 0; i < n; i++) nv[i] = -1;
    for (int v = 0; v < n; v++) if (x[v] < 0) {
        int d = 0;
        for (int u : adj[v]) if (x[u] < 0) {
            tmp[d++] = u;
            nv[u] = v;
            if (d > 4) break;
        }
        if (d == 3 || d == 4) {
            int d2 = 0;
            for (int i = 0; i < d; i++) {
                int a = deg(tmp[i]);
                if (a == 3 || a == 4) tmp[d2++] = tmp[i];
            }
            for (int i = 0; i < d2; i++) {
                int u1 = tmp[i];
                int sB1 = 0;
                used.clear();
                for (int w : adj[u1]) if (x[w] < 0 && w != v) {
                    used.add(w);
                    sB1++;
                }
                for (int j = i + 1; j < d2; j++) {
                    int u2 = tmp[j];
                    if (used.get(u2)) continue;
                    int sB2 = 0;
                    for (int w : adj[u2]) if (x[w] < 0 && w != v && !used.get(w)) sB2++;
                    if (sB1 + sB2 <= 3) {
                        for (int w : adj[u2]) if (x[w] < 0 && used.get(w) && nv[w] != v) {
                            int d3 = deg(w);
                            if (d3 == 3 || d3 == 4) {
                                int sA = d - 2;
                                for (int z : adj[w]) if (x[z] < 0 && z != u1 && z != u2 && nv[z] != v) {
                                    sA++;
                                }
                                if (sA <= 2) {
                                    compute_alternative(vector<int>{v, w}, vector<int>{u1, u2});
                                    goto loop;
                                }
                            }
                        }
                    } // if
                }
            }
        }
loop:      ;
    }
    if (debug >= 3 && depth <= maxDepth && oldn != rn) fprintf(stderr, "%sdesk: %d -> %d\n", debugString().c_str(), oldn, rn);
////       if (oldn != rn) Stat.count("reduceN_desk", oldn - rn);
#endif // 0
    return oldn != rn;
////    }
}

bool branch_and_reduce_algorithm::unconfinedReduction() {
    ////    try (Stat stat = new Stat("reduce_unconfined")) {
    int oldn = rn;
#if 1
    vector<int> &NS = level;
    vector<int> &deg = iter;
////    vector<int> qs(que);
////    vector<int> vs(que);
    for (int v = 0; v < n; v++) if (x[v] < 0) {
////        cout << "Evaluating v=" << v << endl;
        used.clear();
        used.add(v);
        int p = 1, size = 0;
////        cout << "neighbors of " << v << " : ";
        for (int u : adj[v]) if (x[u] < 0) {
            used.add(u);
            NS[size++] = u;
            deg[u] = 1;
////            cout << u << " ";
        }
////        cout << endl << flush;
        bool ok = false;

        while (!ok) {
            ok = true;
            for (int i = 0; i < size; i++) {
////                cout << i << "/" << size << " : ";
                int const u = NS[i];
                if (deg[u] != 1)  {
////                    cout << " continue..." << endl;
                    continue;
                }
                int z = -1;
                for (int const w : adj[u]) if (x[w] < 0 && !used.get(w)) {
                    if (z >= 0) {
                        z = -2;
////                        cout << " break...." << endl;
                        break;
                    }
                    z = w;
////                    cout << ", z=" << w;
                }
                if (z == -1) {
                    if (REDUCTION >= 3) {
                        vector<int> &qs = que;
                        int q = 0;
                        qs[q++] = 1;
                        for (int w : adj[v]) if (x[w] < 0) qs[q++] = w;
                        vector<int> copyOfqs(qs.begin(), qs.begin() + q);
////                        cout << __LINE__ << ": packing vector of size " << copyOfqs.size() << endl << flush;
                        packing.emplace_back(std::move(copyOfqs));
////                        cout << __LINE__ << ", " << this << ", " << depth << ": packing.size=" << packing.size() << endl << flush;
                    }
                    set(v, 1);
                    goto whileloopend;
                } else if (z >= 0) {
                    ok = false;
                    used.add(z);
                    p++;
                    for (int w : adj[z]) if (x[w] < 0) {
                        if (used.add(w)) {
                            NS[size++] = w;
                            deg[w] = 1;
////                            cout << ", add " << w << " to position " << size-1;
                        } else {
                            deg[w]++;
////                            cout << ", increment d[" << w << "]";
                        }
                    }
                }
            }
        }
whileloopend:
        if (x[v] < 0 && p >= 2) {
            used.clear();
            for (int i = 0; i < size; i++) used.add(NS[i]);
            vector<int> &vs = que;
            for (int i = 0; i < size; i++) {
                vs[i] = vs[n + i] = -1;
                int u = NS[i];
                if (deg[u] != 2) continue;
                int v1 = -1, v2 = -1;
                for (int w : adj[u]) if (x[w] < 0 && !used.get(w)) {
                    if (v1 < 0) v1 = w;
                    else if (v2 < 0) v2 = w;
                    else {
                        v1 = v2 = -1;
                        break;
                    }
                }
                if (v1 > v2) {
                    int t = v1;
                    v1 = v2;
                    v2 = t;
                }
                vs[i] = v1;
                vs[n + i] = v2;
            }
            for (int i = 0; i < size; i++) if (vs[i] >= 0 && vs[n + i] >= 0) {
////               cout << debugIndex++ << " ";
                int u = NS[i];
                used.clear();
                for (int w : adj[u]) if (x[w] < 0) used.add(w);
                for (int j = i + 1; j < size; j++) if (vs[i] == vs[j] && vs[n + i] == vs[n + j] && !used.get(NS[j])) {
                    if (REDUCTION >= 3) {
                        vector<int> &qs = que;
                        int q = 0;
                        qs[q++] = 1;
                        for (int w : adj[v]) if (x[w] < 0) qs[q++] = w;
                        vector<int> copyOfqs(qs.begin(), qs.begin() + q);
////                        cout << "packing : ";
////                        for (int const inqs : copyOfqs) {
////                            cout << inqs << " ";
////                        }
////                        cout << endl << flush;
////                        cout << __LINE__ << ": packing vector of size " << copyOfqs.size() << endl << flush;
                        packing.emplace_back(std::move(copyOfqs));
////                        cout << __LINE__ << ", " << this << ", " << depth << ": packing.size=" << packing.size() << endl << flush;
                    }
                    set(v, 1);
                    ////               Stat.count("reduceN_diamond");
                    goto forloopend;
                }
            }
forloopend: ;
////        cout << endl << flush;
        }
    }
    if (debug >= 3 && depth <= maxDepth && oldn != rn) fprintf(stderr, "%sunconfined: %d -> %d\n", debugString().c_str(), oldn, rn);
////    else if (debug >= 3 && depth <= maxDepth) fprintf(stderr, "%sunconfined: %d -> %d\n", debugString().c_str(), oldn, rn);
////    if (oldn != rn) Stat.count("reduceN_unconfined", oldn - rn);
#endif //0
    return oldn != rn;
    ////    }
}

int branch_and_reduce_algorithm::packingReduction() {
    ////    try (Stat stat = new Stat("reduce_packing")) {
    int oldn = rn;
#if 1
    vector<int> x2(x);
    int a = -1;
////    cout << __LINE__ << ", " << this << ", " << depth << ": packing.size=" << packing.size() << endl << flush;
    for (unsigned int pi=0; pi < packing.size(); ++pi) {
        vector<int> &ps = packing[pi];
        if (a != rn) {
            for (int j = 0; j < N; j++) x2[j] = x[j];
            for (int j = modifiedN - 1; j >= 0; j--) {
////                cout << __LINE__ << ", " << this << ", " << depth << ": modifieds[" << j << "]=" << modifieds[j] << endl;
                modifieds[j]->reverse(x2);
            }
////            cout << "Done with loop..." << endl << flush;
            a = rn;
        }
////        cout << __LINE__ << ": " << endl << flush;
////        cout << "ps.size=" << ps.size() << endl << flush;
        int max = ps.size() - 1 - ps[0], sum = 0, size = 0;
        vector<int> &S = level;
////        cout << __LINE__ << ": " << endl << flush;
        for (unsigned int j = 1; j < ps.size(); j++) {
            int v = ps[j];
            if (x2[v] < 0) S[size++] = v;
            if (x2[v] == 1) sum++;
        }
        if (sum > max) {
////            Stat.count("reduceN_packingR");
            return -1;
        } else if (sum == max && size > 0) {
            vector<int> &count = iter;
            used.clear();
////        cout << __LINE__ << ": " << endl << flush;
            for (int j = 0; j < size; j++) {
                used.add(S[j]);
                count[S[j]] = -1;
            }
////        cout << __LINE__ << ": " << endl << flush;
            for (int j = 0; j < size; j++) {
                for (int u : adj[S[j]]) if (x[u] < 0) {
                    if (used.add(u)) {
                        count[u] = 1;
                    } else if (count[u] < 0) {
                        return -1;
                    } else {
                        count[u]++;
                    }
                }
            }
////        cout << __LINE__ << ": " << endl << flush;
            for (int j = 0; j < size; j++) {
                for (int u : adj[S[j]]) if (x[u] < 0 && count[u] == 1) {
                    vector<int> &tmp = que;
                    int p = 0;
                    tmp[p++] = 1;
                    for (int w : adj[u]) if (x[w] < 0 && !used.get(w)) {
                        tmp[p++] = w;
                    }
                    vector<int> copyOfTmp(tmp.begin(), tmp.begin() + p);
////                    cout << __LINE__ << ": packing vector of size " << copyOfTmp.size() << endl << flush;
                    packing.emplace_back(std::move(copyOfTmp));
////            cout << __LINE__ << ", " << this << ", " << depth << ": packing.size=" << packing.size() << endl << flush;
                }
            }
////        cout << __LINE__ << ": " << endl << flush;
            for (int j = 0; j < size; j++) {
                if (S[j] == 1) return -1;
                assert(x[S[j]] < 0);
                set(S[j], 0);
            }
////        cout << __LINE__ << ": " << endl << flush;
        } else if (sum + size > max) {
////        cout << __LINE__ << ": " << endl << flush;
            assert(size >= 2);
            used.clear();
            for (int j = 0; j < size; j++) used.add(S[j]);
            for (int v : adj[S[0]]) if (x[v] < 0 && !used.get(v)) {
                int p = 0;
                for (int u : adj[v]) if (used.get(u)) p++;
                if (sum + p > max) {
                    vector<int> &qs = que;
                    int q = 0;
                    qs[q++] = 2;
                    for (int u : adj[v]) if (x[u] < 0) qs[q++] = u;
                    vector<int> copyOfqs(qs.begin(), qs.begin() + q);
////                    cout << __LINE__ << ": packing vector of size " << copyOfqs.size() << endl << flush;
                    packing.emplace_back(std::move(copyOfqs));
////            cout << __LINE__ << ", " << this << ", " << depth << ": packing.size=" << packing.size() << endl << flush;
                    set(v, 1);
                    break;
                }
            }
////        cout << __LINE__ << ": " << endl << flush;
        }
    }
    if (debug >= 3 && depth <= maxDepth && oldn != rn) fprintf(stderr, "%spacking: %d -> %d\n", debugString().c_str(), oldn, rn);
////    else if (debug >= 3 && depth <= maxDepth) fprintf(stderr, "%spacking: %d -> %d\n", debugString().c_str(), oldn, rn);
    ////        if (oldn != rn) Stat.count("reduceN_packing", oldn - rn);
#endif // 0
    return oldn != rn ? 1 : 0;
    ////    }
}

void branch_and_reduce_algorithm::branching(timer & t, double time_limit) {
    int oldLB = lb;
    int v = -1, dv = 0;
    vector<int> &mirrors = que;
    int mirrorN = 0;
    ////		try (Stat stat = new Stat("branching")) {
    if (BRANCHING == 0) {
        int p = rand()%rn;
        for (int i = 0; i < n; i++) if (x[i] < 0 && p-- == 0) v = i;
        dv = deg(v);
    } else if (BRANCHING == 1) {
        dv = n + 1;
        for (int u = 0; u < n; u++) if (x[u] < 0) {
            int degree = deg(u);
            if (dv > degree) {
                v = u;
                dv = degree;
            }
        }
    } else if (BRANCHING == 2) {
        dv = -1;
        long long minE = 0;
        for (int u = 0; u < n; u++) if (x[u] < 0) {
            int degree = deg(u);
            if (dv > degree) continue;
            long long e = 0;
            used.clear();
            for (int w : adj[u]) if (x[w] < 0) used.add(w);
            for (int w : adj[u]) if (x[w] < 0) {
                for (int w2 : adj[w]) if (x[w2] < 0 && used.get(w2)) e++;
            }
            if (dv < degree || (dv == degree && minE > e)) {
                dv = degree;
                minE = e;
                v = u;
            }
        }
    }
    vector<int> &ps = iter;
    for (int i = 0; i < n; i++) ps[i] = -2;
    used.clear();
    used.add(v);
    for (int u : adj[v]) if (x[u] < 0) {
        used.add(u);
        ps[u] = -1;
    }
    for (int u : adj[v]) if (x[u] < 0) {
        for (int w : adj[u]) if (x[w] < 0 && used.add(w)) {
            int c1 = dv;
            for (int z : adj[w]) if (x[z] < 0 && ps[z] != -2) {
                ps[z] = w;
                c1--;
            }
            bool ok = true;
            for (int u2 : adj[v]) if (x[u2] < 0 && ps[u2] != w) {
                int c2 = 0;
                for (int w2 : adj[u2]) if (x[w2] < 0 && ps[w2] == w) c2++;
                if (c2 != c1 - 1) {
                    ok = false;
                    break;
                }
            }
            if (ok) mirrors[mirrorN++] = w;
        }
    }
    ////		}
    int pn = rn;
    unsigned int oldP = packing.size();
    if (REDUCTION >= 3) {
        vector<int> &tmp = level;
        int p = 0;
        tmp[p++] = mirrorN > 0 ? 2 : 1;
        for (int u : adj[v]) if (x[u] < 0) tmp[p++] = u;
        vector<int> copyOfTmp(tmp.begin(), tmp.begin() + p);
////        cout << __LINE__ << ": packing vector of size " << copyOfTmp.size() << endl << flush;
        packing.emplace_back(std::move(copyOfTmp));
////            cout << __LINE__ << ", " << this << ", " << depth << ": packing.size=" << packing.size() << endl << flush;
    }
    set(v, 1);
    for (int i = 0; i < mirrorN; i++) set(mirrors[i], 1);
    if (debug >= 2 && depth <= maxDepth) {
        if (mirrorN > 0) fprintf(stderr, "%sbranchMirror (%d, %d): 1\n", debugString().c_str(), dv, mirrorN);
        else fprintf(stderr, "%sbranch (%d): 1\n", debugString().c_str(), dv);
    }
    depth++;
    rec(t, time_limit);
    while (packing.size() > oldP) packing.pop_back();
////    cout << __LINE__ << ", " << this << ", " << depth << ": packing.size=" << packing.size() << endl << flush;
    lb = oldLB;
    depth--;
    restore(pn);
    if (lb >= opt) return;
    nBranchings++;
    if (mirrorN == 0) {
        used.clear();
        used.add(v);
        for (int u : adj[v]) if (x[u] < 0) used.add(u);
        if (REDUCTION >= 3) {
            vector<int> ws(n, -1);
            for (int u : adj[v]) if (x[u] < 0) {
                vector<int> &tmp = level;
                int p = 0;
                tmp[p++] = 1;
                for (int w : adj[u]) if (x[w] < 0 && !used.get(w)) {
                    tmp[p++] = w;
                    ws[w] = u;
                }
                assert(p >= 2);
                for (int u2 : adj[tmp[1]]) if (x[u2] < 0 && used.get(u2) && u2 != u) {
                    int c = 0;
                    for (int w : adj[u2]) if (x[w] < 0) {
                        if (ws[w] == u) c++;
                        else if (w == u || !used.get(w)) {
                            c = -1;
                            break;
                        }
                    }
                    if (c == p - 1) {
                        tmp[0] = 2;
                        break;
                    }
                }
                vector<int> copyOfTmp(tmp.begin(), tmp.begin() + p);
////                    cout << __LINE__ << ": packing vector of size " << copyOfTmp.size() << endl << flush;
                packing.emplace_back(std::move(copyOfTmp));
////            cout << __LINE__ << ", " << this << ", " << depth << ": packing.size=" << packing.size() << endl << flush;
            }
        }
    }
    set(v, 0);
    if (debug >= 2 && depth <= maxDepth) fprintf(stderr, "%sbranch (%d): 0\n", debugString().c_str(), dv);
    depth++;
    rec(t, time_limit);
    while (packing.size() > oldP) packing.pop_back();
////    cout << __LINE__ << ", " << this << ", " << depth << ": packing.size=" << packing.size() << endl << flush;
    lb = oldLB;
    depth--;
    restore(pn);
}

bool branch_and_reduce_algorithm::decompose(timer & t, double time_limit) {
    vector<vector<int>> vss;
////    try (Stat stat = new Stat("decompose")) {
        {
        vector<int> &id = level;
        vector<int> &size = iter;
        int nC = 0;
        {
            for (int i = 0; i < n; i++) id[i] = -1;
            for (int s = 0; s < n; s++) if (x[s] < 0 && id[s] < 0) {
                nC++;
                int qs = 0, qt = 0;
                que[qt++] = s;
                id[s] = s;
                while (qs < qt) {
                    int v = que[qs++];
                    for (int u : adj[v]) if (x[u] < 0 && id[u] < 0) {
                        id[u] = s;
                        que[qt++] = u;
                    }
                }
                size[s] = qt;
            }
        }
        if (nC <= 1 && (n <= 100 || n * SHRINK < rn)) return false;
        vector<long long> cs(nC, 0);
        {
            int p = 0;
            for (int i = 0; i < n; i++) if (x[i] < 0 && id[i] == i) {
                cs[p++] = ((long long)(size[i])) << 32 | i;
            }
            std::sort(cs.begin(), cs.end());
        }
        vss.resize(nC);
        vector<int> qs(n, 0);
        {
            for (int i = 0; i < nC; i++) {
                vss[i].resize(size[(int)cs[i]]);
                qs[(int)cs[i]] = i;
            }
            vector<int> ps(nC);
            for (int i = 0; i < n; i++) if (x[i] < 0) {
                int j = qs[id[i]];
                vss[j][ps[j]++] = i;
            }
        }
        for (int i = 0; i < n; i++) id[i] = -1;
        for (unsigned int i = 0; i < vss.size(); i++) {
            vector<int> &vs = vss[i];
            vector<long long> ls(vs.size());
            for (unsigned int j = 0; j < vs.size(); j++) ls[j] = ((long long)(n - deg(vs[j]))) << 32 | vs[j];
            std::sort(ls.begin(), ls.end());
            for (unsigned int j = 0; j < vs.size(); j++) vs[j] = (int)ls[j];
        }
    }
    vector<int> x2(x);
    for (int i = modifiedN - 1; i >= 0; i--) modifieds[i]->reverse(x2);
    vector<int> size(vss.size());
    for (unsigned int i = 0; i < vss.size(); i++) size[i] = vss[i].size();
    vector<int> pos1(N, -1);
    vector<int> pos2(N, 0);
    vector<vector<int>> packingB;
    {
        for (unsigned int i = 0; i < vss.size(); i++) {
            for (unsigned int j = 0; j < vss[i].size(); j++) {
                pos1[vss[i][j]] = i;
                pos2[vss[i][j]] = j;
            }
        }
        vector<bool> need(N, false);
        for (vector<int> const &ps : packing) {
            
            int max = ps.size() - 1 - ps[0], sum = 0, count = 0;
            for (unsigned int j = 1; j < ps.size(); j++) {
                int v = ps[j];
                if (x2[v] < 0 || x2[v] == 2) {
                    count++;
                }
                if (x2[v] == 1) sum++;
            }
            if (sum > max) return true;
            if (sum + count > max) {
                packingB.push_back(ps);
                for (unsigned int k = 1; k < ps.size(); k++) {
                    if (x2[ps[k]] == 2) need[ps[k]] = true;
                }
            }
        }
        for (int i = 0; i < modifiedN; i++) {
            bool b = false;
            shared_ptr<modified> mod = modifieds[i];
            for (int v : mod->removed) if (need[v]) b = true;
            if (b) {
                if (dynamic_cast<fold*>(mod.get()) != nullptr) {
                    if (x2[mod->vs[0]] == 2) need[mod->vs[0]] = true;
                } else {
                    for (int v : mod->vs) if (x2[v] == 2) {
                        need[v] = true;
                    }
                }
            }
        }
        for (int i = modifiedN - 1; i >= 0; i--) {
            shared_ptr<modified> mod = modifieds[i];
            bool b = false;
            for (int v : mod->removed) if (need[v]) b = true;
            if (b) {
                if (dynamic_cast<fold*>(mod.get()) != nullptr) {
                    for (int v : mod->removed) {
                        assert(pos1[v] == -1);
                        pos1[v] = pos1[mod->vs[0]];
                        assert(pos1[v] >= 0);
                        pos2[v] = size[pos1[v]]++;
                    }
                } else {
                    int max = -1;
                    for (int v : mod->vs) if (max < pos1[v]) max = pos1[v];
                    assert(max >= 0);
                    for (int v : mod->removed) {
                        assert(pos1[v] == -1);
                        pos1[v] = max;
                        pos2[v] = size[pos1[v]]++;
                    }
                }
            }
        }
        for (int i = 0; i < n; i++) {
            if ((x2[i] == 0 || x2[i] == 1) && pos1[i] >= 0) {
////                Debug.print(i, n, x[i]);
                assert(false);
            }
        }
    }
    vector<branch_and_reduce_algorithm*> vcs(vss.size(), nullptr);
    {
        for (int i = 0; i < static_cast<int>(vss.size()); i++) {
            vector<int> &vs = vss[i];
            size[i] += 2;
            vector<vector<int>> adj2(vs.size());
            for (int j = 0; j < static_cast<int>(vs.size()); j++) {
                adj2[j].resize(deg(vs[j]), 0);
                int p = 0;
                for (int u : adj[vs[j]]) if (x[u] < 0) adj2[j][p++] = pos2[u];
                assert(p == adj2[j].size());
                std::sort(adj2[j].begin(), adj2[j].end());
            }

            vcs[i] = new branch_and_reduce_algorithm(adj2, size[i]);
            for (unsigned int j = 0; j < vs.size(); j++) {
                if (in[vs[j]] >= 0 && pos1[in[vs[j]]] == i && pos2[in[vs[j]]] < static_cast<int>(vs.size())) {
                    vcs[i]->in[j] = pos2[in[vs[j]]];
                }
                if (out[vs[j]] >= 0 && pos1[out[vs[j]]] == i && pos2[out[vs[j]]] < static_cast<int>(vs.size())) {
                    vcs[i]->out[j] = pos2[out[vs[j]]];
                }
            }
            vcs[i]->x[vcs[i]->N - 2] = vcs[i]->y[vcs[i]->N - 2] = 0;
            vcs[i]->x[vcs[i]->N - 1] = vcs[i]->y[vcs[i]->N - 1] = 1;
        }
    }
    {
        for (unsigned int i = 0; i < packingB.size(); i++) {
            vector<int> &ps(packingB[i]);
            int maxID = -1;
            for (unsigned int j = 1; j < ps.size(); j++) {
                int v = ps[j];
                if (x2[v] < 0 || x2[v] == 2) {
                    maxID = std::max(maxID, pos1[v]);
                }
            }
            vcs[maxID]->packing.push_back(ps);
////            cout << __LINE__ << ": packing vector of size " << ps.size() << endl << flush;
////            cout << __LINE__ << ", " << this << ", " << depth << ": vcs.packing.size=" << vcs[maxID]->packing.size() << endl << flush;
        }
    }
    {
        for (int i = 0; i < modifiedN; i++) {
            shared_ptr<modified> mod = modifieds[i];
            int p = pos1[mod->removed[0]];
            if (p >= 0) {
                vcs[p]->modifieds[vcs[p]->modifiedN++] = mod;
////                cout << __LINE__ << ", " << this << ", " << depth << ": Setting modifieds[" << vcs[p]->modifiedN-1  << "]=" << modifieds[vcs[p]->modifiedN-1] << endl << flush;
            }
        }
    }
    vector<vector<int>> vss2(vss.size());
    {
        for (unsigned int i = 0; i < vss.size(); i++) vss2[i].resize(vcs[i]->N - 2, 0);
        for (int i = 0; i < N; i++) if (pos1[i] >= 0) vss2[pos1[i]][pos2[i]] = i;
    }
    int sum = crt;
    for (int i = 0; i < static_cast<int>(vss.size()) && opt > sum; i++) {
        branch_and_reduce_algorithm *vc = vcs[i];
        {
            vector<vector<int>> packing2;
////            list<vector<int>> packing2;
            for (vector<int> const &ps : vc->packing) {
                vector<int> &tmp = level;
                int p = 0;
                tmp[p++] = ps[0];
                for (unsigned int k = 1; k < ps.size(); k++) {
                    int v = ps[k];
                    if (pos1[v] == i) {
                        tmp[p++] = pos2[v];
                    } else {
                        assert(x2[v] == 0 || x2[v] == 1);
                        if (x2[v] == 0) tmp[0]--;
                    }
                }
                if (p - 1 < tmp[0]) return true;
                if (tmp[0] <= 0) continue;
                vector<int> copyOfTmp(tmp.begin(), tmp.begin() + p);
////                cout << __LINE__ << ": packing vector of size " << copyOfTmp.size() << endl << flush;
                packing2.emplace_back(std::move(copyOfTmp));
////            cout << __LINE__ << ", " << this << ", " << depth << ": vc.packing.size=" << packing.size() << endl << flush;
            }
            (vc->packing).swap(packing2);
        }
        {
            for (int j = 0; j < vc->modifiedN; j++) {
                shared_ptr<modified> mod = vc->modifieds[j];
                vector<int> removed(mod->removed.size());
                for (unsigned int k = 0; k < removed.size(); k++) {
                    int v = mod->removed[k];
                    assert(pos1[v] == i);
                    removed[k] = pos2[v];
                }
                if (dynamic_cast<fold*>(mod.get()) != nullptr) {
                    vector<int> vs(1, 0);
                    int v = mod->vs[0];
                    if (pos1[v] == i) {
                        vs[0] = pos2[v];
                    } else {
                        assert(x2[v] == 0 || x2[v] == 1);
                        vs[0] = vc->N - 2 + x2[v];
                    }
                    mod = make_shared<fold>(fold(removed, vs, this));
                } else {
                    vector<int> vs(mod->vs.size(), 0);
                    for (unsigned int k = 0; k < vs.size(); k++) {
                        int v = mod->vs[k];
                        if (pos1[v] == i) {
                            vs[k] = pos2[v];
                        } else {
                            assert(x2[v] == 0 || x2[v] == 1);
                            vs[k] = vc->N - 2 + x2[v];
                        }
                    }
                    mod = make_shared<alternative>(alternative(removed, vs, this, dynamic_cast<alternative*>(mod.get())->k));
                }
////                cout << __LINE__ << ", " << vc << ": Deleting " << vc->modifieds[j] << endl << flush;
////                delete vc->modifieds[j];
                vc->modifieds[j] = mod;
////                cout << __LINE__ << ", " << vc << ": Setting modifieds[" << j << "]=" << vc->modifieds[j] << endl << flush;
            }
        }
        vc->depth = depth + (vss.size() > 1 ? 1 : 0);
        if (debug >= 2 && depth <= maxDepth) {
            if (vss.size() == 1) fprintf(stderr, "%sshrink: %d -> %d (%d)\n", debugString().c_str(), n, vcs[i]->n, vcs[i]->N);
            else fprintf(stderr, "%sdecompose: %d (%d)\n", debugString().c_str(), vcs[i]->n, vcs[i]->N);
        }
        if (i + 1 == static_cast<int>(vss.size())) {
            vc->opt = std::min(vss[i].size(), static_cast<size_t>(opt - sum));
        }
        vc->reverse();
        for (int j = 0; j < vc->N; j++) assert(vc->y[j] == 0 || vc->y[j] == 1);
        vc->solve(t, time_limit);
        sum += vc->opt;
        for (int j = 0; j < vc->N - 2; j++) {
            x2[vss2[i][j]] = vc->y[j];
            assert(vc->y[j] == 0 || vc->y[j] == 1);
        }
    }
    if (opt > sum) {
        if (debug >= 2 && rootDepth <= maxDepth) fprintf(stderr, "%sopt: %d -> %d\n", debugString().c_str(), opt, sum);
        opt = sum;
        y = x;
        for (unsigned int i = 0; i < vss.size(); i++) {
            for (unsigned int j = 0; j < vss[i].size(); j++) y[vss[i][j]] = vcs[i]->y[j];
        }
        reverse();
    }

    for (branch_and_reduce_algorithm *pAlg : vcs) {
        delete pAlg;
        pAlg = nullptr;
    }

    return true;
}

bool branch_and_reduce_algorithm::reduce() {
    int oldn = rn;
    for (;;) {
        if (REDUCTION >= 0) deg1Reduction();
        // if (n > 100 && n * SHRINK >= rn && !outputLP && decompose()) return true;
        if (REDUCTION >= 0 && REDUCTION < 2 && dominateReduction()) continue;
        if (REDUCTION >= 2 && unconfinedReduction()) continue;
        if (REDUCTION >= 1 && lpReduction()) continue;
        if (REDUCTION >= 3) {
            int r = packingReduction();
            if (r < 0) return true;
            if (r > 0) continue;
        }
        if (REDUCTION >= 1 && fold2Reduction()) continue;
        if (REDUCTION >= 2 && twinReduction()) continue;
        if (REDUCTION >= 2 && funnelReduction()) continue;
        if (REDUCTION >= 2 && deskReduction()) continue;
        break;
    }
    if (debug >= 2 && depth <= maxDepth && oldn != rn) fprintf(stderr, "%sreduce: %d -> %d\n", debugString().c_str(), oldn, rn);
////    else if (debug >= 2 && depth <= maxDepth) fprintf(stderr, "%sreduce: %d -> %d\n", debugString().c_str(), oldn, rn);
    return false;
}

void branch_and_reduce_algorithm::rec(timer & t, double time_limit) {
    if (t.elapsed() >= time_limit) return;
    if (REDUCTION < 3) assert(packing.size() == 0);
    if (reduce()) return;
    if (lowerBound() >= opt) return;
    if (rn == 0) {
        if (debug >= 2 && rootDepth <= maxDepth) fprintf(stderr, "%sopt: %d -> %d\n", debugString().c_str(), opt, crt);
        opt = crt;
        y = x;
        reverse();
        return;
    }
    if (decompose(t, time_limit)) return;
    branching(t, time_limit);
}


int branch_and_reduce_algorithm::solve(timer & t, double time_limit) {
    if (t.elapsed() >= time_limit) return -1;
////    cout << "level.size=" << level.size() << endl << flush;
    // PrintState();
    if (LOWER_BOUND >= 2 && REDUCTION <= 0 && !outputLP) {
        cerr << "LP/cycle lower bounds require LP reduction." << endl << flush;
        assert(0);
    }
    rootDepth = depth;
    if (outputLP) {
        if (REDUCTION < 0) {
            lpReduction();
        } else {
            reduce();
        }
        printf("%.1f\n", crt + rn / 2.0);
        return opt;
    }
    rec(t, time_limit);
    if (debug >= 2 && depth <= maxDepth) fprintf(stderr, "%sopt: %d\n", debugString().c_str(), opt);
    if (t.elapsed() >= time_limit) return -1;
    else return opt;
}

std::string branch_and_reduce_algorithm::debugString() const {
    stringstream ins;
#ifdef PUT_TIME
    time_t rawtime;
    struct tm *timeinfo;

    time (&rawtime);
    timeinfo = localtime (&rawtime);
    ins << std::put_time(timeinfo, "%T") << "  ";
#else
    std::locale::global(std::locale("ja_JP.utf8"));
    std::time_t t = std::time(NULL);
    char mbstr[100];
    if (std::strftime(mbstr, sizeof(mbstr), "%T", std::localtime(&t))) {
        std::cout << mbstr << '\n';
    }
#endif
    for (int i = 0; i < depth && i < maxDepth; ++i) {
        ins << " ";
    }
////    printf ("Current local time and date: %s", asctime(timeinfo));
    return ins.str();

}

void branch_and_reduce_algorithm::PrintState() const
{
    cout << "State(" << this << "):" << endl << flush;
    cout << "adj=" << endl << flush;
    for (unsigned int j = 0; j < adj.size(); ++j) {
        cout << j << " : ";
        for (int const k : adj[j]) {
            cout << k << " ";
        }
        cout << endl;
    }
    cout << "N  =" << N << endl << flush;
    cout << "in =";
    for (int const i : in) {
        cout << i << " ";
    }
    cout << endl << flush;
    cout << "out=";
    for (int const i : out) {
        cout << i << " ";
    }
    cout << endl << flush;
}


void branch_and_reduce_algorithm::reduce_graph()
{
    int oldn = rn;
    for (;;) {
        if (REDUCTION >= 0) deg1Reduction();
////        if (n > 100 && n * SHRINK >= rn && !outputLP && decompose()) return true;
        if (REDUCTION >= 0 && REDUCTION < 2 && dominateReduction()) continue;
        if (REDUCTION >= 2 && unconfinedReduction()) continue;
        if (REDUCTION >= 1 && lpReduction()) continue;
        if (REDUCTION >= 3) {
            int r = packingReduction();
            if (r < 0) return;
            if (r > 0) continue;
        }
        if (REDUCTION >= 1 && fold2Reduction()) continue;
        if (REDUCTION >= 2 && twinReduction()) continue;
        if (REDUCTION >= 2 && funnelReduction()) continue;
        if (REDUCTION >= 2 && deskReduction()) continue;
        break;
    }
////    opt = crt;
    if (debug >= 2 && depth <= maxDepth && oldn != rn) fprintf(stderr, "%sreduce: %d -> %d\n", debugString().c_str(), oldn, rn);
    size_t low_degree_count(0);
    for (int v = 0; v < n; v++) if (x[v] < 0) {
        if (deg(v) <= 1) {
            low_degree_count++;
        }
    }
    cout << "There are " << low_degree_count << " degree 0 and 1 vertices left!" << endl << flush;
}

void branch_and_reduce_algorithm::initial_reduce_graph()
{
    reduce_graph();
    snapshotX = x;
    reductionSnapshotSize = modifieds.size();
}

size_t branch_and_reduce_algorithm::get_current_is_size() const {

    vector<int> x2(x);

    for (int i = modifiedN - 1; i >= 0; i--) {
////        cout << __LINE__ << ", " << this << ": modifieds[" << i << "] is " << modifieds[i] << endl << flush;
        modifieds[i]->reverse(x2);
    }

    size_t current_is_size(0);
    for (unsigned int i = 0; i < adj.size(); ++i) {
        if (x2[i] == 0) {
            current_is_size++;
        }
    }

    return current_is_size;
}

size_t branch_and_reduce_algorithm::get_current_is_size_with_folds() const {

    size_t folded_vertex_count(0);
    size_t current_is_size(0);
    for (int const i : x) {
        if (i == 0) current_is_size++;
        if (i == 2) folded_vertex_count++;
    }

    return current_is_size + folded_vertex_count/2;
}

bool branch_and_reduce_algorithm::folded_vertices_exist() const {

    vector<int> x2(x);

    for (int i = modifiedN - 1; i >= 0; i--) {
////        cout << __LINE__ << ", " << this << ": modifieds[" << i << "] is " << modifieds[i] << endl << flush;
        modifieds[i]->reverse(x2);
    }

    for (int const i : x2) {
        if (i == 2) return true;
    }

    return false;
}

vector<int> branch_and_reduce_algorithm::compute_maximal_is() {

    int vertexToForceInIndependentSet(0);
    while (vertexToForceInIndependentSet != -1) {
        reduce_graph();

        vertexToForceInIndependentSet = -1;
        for (unsigned int i = 0; i < x.size(); ++i) {
            if (x[i] == -1) { // status not determined
                vertexToForceInIndependentSet = i;
                break;
            }
        }

        // add vertex to independent set
        if (vertexToForceInIndependentSet != -1) {
            set(vertexToForceInIndependentSet, 0);
        }
    }

    vector<int> x2(x);

    for (int i = modifiedN - 1; i >= 0; i--) {
////        cout << __LINE__ << ", " << this << ": modifieds[" << i << "] is " << modifieds[i] << endl << flush;
        modifieds[i]->reverse(x2);
    }

    size_t current_is_size(0);
    for (unsigned int i = 0; i < adj.size(); ++i) {
        if (x2[i] == 0) {
            current_is_size++;
        }
    }

    return x2;
}

size_t branch_and_reduce_algorithm::compute_alternative_maximal_is_size() {

    int vertexToForceInIndependentSet(0);
    while (vertexToForceInIndependentSet != -1) {
        reduce_graph();

        vertexToForceInIndependentSet = -1;
        for (int i = 0; i < static_cast<int>(x.size()); ++i) {
            if (x[i] == -1) { // status not determined
                vertexToForceInIndependentSet = i;
                break;
            }
        }

        // add vertex to independent set
        if (vertexToForceInIndependentSet != -1) {
            set(vertexToForceInIndependentSet, 0);
        }
    }

    size_t numberOfFoldedVertices(0);
    size_t sizeOfIS(0);
    for (int const i : x) {
        if (i == 0) sizeOfIS++;
        if (i == 2) numberOfFoldedVertices++;
    }

    return sizeOfIS + numberOfFoldedVertices/2;
}

void branch_and_reduce_algorithm::restore_to_snapshot() {

    while (modifiedN > reductionSnapshotSize) {
        modifieds[--modifiedN]->restore();
        modifieds[modifiedN] = shared_ptr<modified>();
    }

    x = snapshotX;
}

size_t branch_and_reduce_algorithm::number_of_nodes_remaining() const {

    size_t node_count(0);
    for (int i : x) if (i == -1) node_count++;

    return node_count;
}

void branch_and_reduce_algorithm::force_into_independent_set(vector<NodeID> const &nodes) {

    for (NodeID const node : nodes) {
        assert(x[node] == -1); // should not have been assigned yet.
        if (x[node] != -1) {
            cout << "ERROR: invalid vertex selected for independent set!" << endl << flush;
        }
////        cout << "Switching node " << node << " from value " << x[node] << " to 0" << endl;
        set(node, 0);
    }
}

// input: vector mapping vertices in independent set to true

// WARNING: this is destructive, can only be applied once, and then the class is not
//          guaranteed to be a good state!
void branch_and_reduce_algorithm::extend_finer_is(vector<bool> & independent_set)
{

////    cout << "Extending finer independent set" << endl << flush;
    assert(independent_set.size() == adj.size());
    assert(independent_set.size() == x.size());

    // first fill in temporary vector with independent set vertices from finer set
    size_t new_independent_set_vertices(0);
    for (size_t index = 0; index < independent_set.size(); ++index) if (independent_set[index]) {
        assert(x[index] == -1); // it needs to be unassigned
        if (x[index] != -1) {
            cout << "ERROR: invalid vertex selected for independent set!" << endl << flush;
        }
        set(index, 0); // add to independent set
        new_independent_set_vertices++;
    }

    vector<int> x2(x);

    // undo reductions
    for (int i = modifiedN - 1; i >= 0; i--) {
////        cout << __LINE__ << ", " << this << ": modifieds[" << i << "] is " << modifieds[i] << endl << flush;
        modifieds[i]->reverse(x2);
    }

////    size_t independent_set_size(0);
    // update full independent set
    for (unsigned int i = 0; i < adj.size(); ++i) if (x2[i] == 0) {
        independent_set[i] = true;
////        independent_set_size++;
    }

////    if (independent_set_size != new_independent_set_vertices + get_current_is_size()) {
////        cout << "ERROR: incorrect original count for independent set size with reductions!" << endl << flush;
////    }
}

void branch_and_reduce_algorithm::get_solved_is(vector<bool> & independent_set) {
    for (unsigned int i = 0; i < y.size(); ++i) {
        if (y[i] == 0) independent_set[i] = true;
    }
}

void branch_and_reduce_algorithm::convert_adj_lists(graph_access & G, std::vector<NodeID> & reverse_mapping) const {
    // Number of nodes
    unsigned int const node_count = number_of_nodes_remaining();
    // Number of edges
    int m = 0;

    // Nodes -> Range
    std::vector<NodeID> mapping(adj.size(), UINT_MAX);

    // Get number of edges and reorder nodes
    unsigned int node_counter = 0;
    for (NodeID node = 0; node < adj.size(); ++node) if (x[node] < 0) {
        for (int const neighbor : adj[node]) if (x[neighbor] < 0) m++;
        mapping[node] = node_counter;
        reverse_mapping[node_counter] = node;
        node_counter++;
    }

    // Create the adjacency array
    std::vector<int> xadj(node_count + 1);
    std::vector<int> adjncy(m);
    unsigned int adjncy_counter = 0;
    for (unsigned int i = 0; i < node_count; ++i) {
        xadj[i] = adjncy_counter;
        for (int const neighbor : adj[reverse_mapping[i]]) {
            if (mapping[neighbor] == i) continue;
            if (mapping[neighbor] == UINT_MAX) continue;
            adjncy[adjncy_counter++] = mapping[neighbor];
        }
        std::sort(std::begin(adjncy) + xadj[i], std::begin(adjncy) + adjncy_counter);
    }
    xadj[node_count] = adjncy_counter;

    // Build the graph
    G.build_from_metis(node_count, &xadj[0], &adjncy[0]);

#if 0
    std::vector<std::set<int>> neighbors;
    neighbors.resize(G.number_of_nodes());
    // verify the graph
    forall_nodes(G, node) {
        forall_out_edges(G, edge, node) {
            NodeID neighbor = G.getEdgeTarget(edge);
            neighbors[node].insert(neighbor);    
        } endfor
    } endfor

    for (int vertex = 0; vertex < node_count; ++vertex) {
        size_t neighbors_in_subgraph(0);
        for (int const neighbor : adj[reverse_mapping[vertex]]) {
            bool const in_mapping(mapping[neighbor] != UINT_MAX);
            if (in_mapping) {
                bool const in_graph(neighbors[vertex].find(mapping[neighbor]) != neighbors[vertex].end());
                if (in_graph) neighbors_in_subgraph++;
            }
        }

        if (neighbors_in_subgraph != neighbors[vertex].size()) {
            cout << "ERROR: subgraph verification failed" << endl << flush;
        }
    }
#endif // DEBUG
}

