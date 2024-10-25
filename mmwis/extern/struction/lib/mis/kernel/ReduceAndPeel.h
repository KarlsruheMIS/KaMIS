//
// Created by alex on 16.11.19.
//

#ifndef LINEARTIMEKERNELIZATION_REDUCER_H
#define LINEARTIMEKERNELIZATION_REDUCER_H

#include <stack>
#include <algorithm>
#include <numeric>
#include <array>
#include "graph_access.h"
#include "timer.h"

namespace struction {

struct RestoreData {
    explicit RestoreData(std::vector<NodeID> &path) : path(path) {}
    std::vector<NodeID> path;
    std::vector<std::pair<NodeID, NodeWeight>> weights;
    std::array<NodeID, 4> unlink_nodes;
};
class ReduceAndPeel {
public:
    ReduceAndPeel(mmwis::graph_access &g) : g(g), degrees(g.number_of_nodes()), weight_diffs(g.number_of_nodes())
            , marked(g.number_of_nodes(), false), reduction_offset(0), active_nodes(g.number_of_nodes()), choices(g.number_of_nodes())
            , edges(g.number_of_edges()){
        std::iota(active_nodes.begin(), active_nodes.end(), 0);
        init(g);
    }

    void reduceAndPeel() {
        reduce();
        //calc_reduced_graph();
        while (peel()) {
            reduce();
        }
        calculate_IS();
    }

    void reduce() {
        while (true) {
            NodeID n;
            if (next_vertex(neighborhood_removal_vertices, n)) {
                reduce_neighborhood_vertex(n);
            } else if (next_vertex(deg_1_vertices, n)) {
                reduce_degree_one_vertex(n);
            } else if (next_vertex(deg_2_vertices, n)) {
                reduce_path(n);
            } else
                return;
        }
    }

    bool peel(double k=0.01) {
        update_active_nodes();
        if (active_nodes.empty()) return false;
        size_t lim = std::ceil(k * active_nodes.size());
        std::nth_element(active_nodes.begin(), active_nodes.begin() + lim, active_nodes.end(), [&](NodeID a, NodeID b) {
            //return degrees[a] > degrees[b] || (degrees[a] == degrees[b] && g.getNodeWeight(a) < g.getNodeWeight(b));
            //return (float) g.getNodeWeight(a) / degrees[a] < (float) g.getNodeWeight(b) / degrees[b];
            //return g.getNodeWeight(a) < g.getNodeWeight(b);
            //return (float) g.getNodeWeight(a) / (weight_diffs[a] - g.getNodeWeight(a)) < (float) g.getNodeWeight(b) / (weight_diffs[b] - g.getNodeWeight(b));
            //return weight_diffs[a] - g.getNodeWeight(a) < weight_diffs[b] - g.getNodeWeight(b);
            return weight_diffs[a] < weight_diffs[b];
        });
        for (size_t i = 0; i < lim; ++i)
            delete_vertex(active_nodes[i]);

        return true;
    }


    mmwis::graph_access &reducedGraph() {
        return reduced_graph;
    }

    size_t reductionOffset() {
        return reduction_offset;
    }

private:
    mmwis::graph_access &g;
    mmwis::graph_access reduced_graph;
    std::vector<size_t> degrees;
    std::vector<int> weight_diffs;
    std::stack<NodeID> deg_1_vertices;
    std::stack<NodeID> deg_2_vertices;
    std::stack<NodeID> neighborhood_removal_vertices;
    std::vector<bool> marked;
    std::vector<NodeID> active_nodes;
    size_t edges;

    std::vector<bool> choices;
    std::stack<RestoreData> restore_stack;

    std::vector<NodeID> path;

    size_t reduction_offset;

    void init(mmwis::graph_access &g) {
        forall_nodes(g, n)
        degrees[n] = g.getNodeDegree(n);
        weight_diffs[n] = get_weight_dif(n);
        try_append_vertex(n);
        endfor
    }

    int get_weight_dif(NodeID n) {
        int weight_dif = g.getNodeWeight(n);
        forall_out_edges(g, e, n)
        if (marked[g.getEdgeTarget(e)]) continue;
        weight_dif -= g.getNodeWeight(g.getEdgeTarget(e));
        endfor
        return weight_dif;
    }

    void set_node_weight(NodeID n, NodeWeight w) {
        int weight_dif = g.getNodeWeight(n) - w;

        g.setNodeWeight(n, w);
        weight_diffs[n] -= weight_dif;
        if (weight_diffs[n] >= 0)
            neighborhood_removal_vertices.push(n);
        forall_out_edges(g, e, n)
        weight_diffs[g.getEdgeTarget(e)] += weight_dif;
        if (weight_diffs[g.getEdgeTarget(e)] >= 0)
            neighborhood_removal_vertices.push(g.getEdgeTarget(e));
        endfor
    }

    void increment_reduction_offset(NodeWeight w) {
        reduction_offset += w;
        //std::cout << timer.elapsed() << "," << reduction_offset << std::endl;
    }

    inline void reduce_neighborhood_vertex(NodeID n) {
        if (weight_diffs[n] < 0) return;
        increment_reduction_offset(g.getNodeWeight(n));
        marked[n] = true;
        g.setPartitionIndex(n, true);
        delete_neighborhood(n);
    }

    inline void reduce_degree_one_vertex(NodeID n) {
        if (degrees[n] != 1) return;
        delete_vertex(n);
        increment_reduction_offset(g.getNodeWeight(n));
        forall_out_edges(g, e, n)
        NodeID u = g.getEdgeTarget(e);
        if (marked[u]) continue;
        if (g.getNodeWeight(n) >= g.getNodeWeight(u)) {
            delete_vertex(u);
            g.setPartitionIndex(n, true);
        } else {
            path.clear();
            path.push_back(u);
            path.push_back(n);
            path.push_back(u);
            append_restore_data({u});
            set_node_weight(u, g.getNodeWeight(u) - g.getNodeWeight(n));
        }

        break;
        endfor
    }

    void reduce_path(NodeID n) {
        if (degrees[n] != 2) return;
        find_max_path(n);

        NodeID v = path.front(), w = path.back();
        NodeWeight w_v = g.getNodeWeight(v), w_w = g.getNodeWeight(w);

        NodeWeight w_i_i,w_i_e,w_e_i,w_e_e;
        find_MIS_on_path(w_i_i, w_i_e, w_e_i, w_e_e);
        if (v == w) {
            NodeWeight w_e = w_e_e;
            NodeWeight w_i = w_i_i - g.getNodeWeight(v);
            increment_reduction_offset(degrees[v] == 2 ? std::max(w_i, w_e) : w_e);

            bool fold = degrees[v] != 2 && w_i > w_e;
            delete_path(path, fold, path.size() - 1);
            if (fold) {
                //fold
                append_restore_data({v});
                set_node_weight(v, w_i - w_e);
            } else {
                g.setPartitionIndex(v, w_i > w_e);
                apply(path);
            }
        } else {
            bool connected = are_connected(v, w);
            if (connected || w_i_i <= w_i_e || w_i_i <= w_e_i) {
                bool keep_v = w_i_e > w_e_e, keep_w = w_e_i > w_e_e;
                increment_reduction_offset(w_e_e);
                delete_path(path, keep_v, path.size() - keep_w);
                if (!keep_v && !keep_w) {
                    apply(path);
                } else {
                    if (!connected && keep_v && keep_w) {
                        relink(v, w, path[1], path[path.size() - 2]);
                        append_restore_data({v, w}, v, w, path[1], path[path.size() - 2]);
                    } else {
                        append_restore_data({v, w});
                    }
                    if (keep_v) set_node_weight(v, w_i_e - w_e_e);
                    if (keep_w) set_node_weight(w, w_e_i - w_e_e);
                }
            } else {
                int gap = static_cast<int>(w_e_e) - w_e_i - w_i_e + w_i_i;
                NodeID c = path[1];
                if (gap >= 0) {
                    if (path.size() <= 3) return;
                    append_restore_data({c, v, w}, c, w, path[2], path[path.size() - 2]);
                    increment_reduction_offset(w_e_e - gap);
                    delete_path(path, 2, path.size() - 1);
                    relink(c, w, path[2], path[path.size() - 2]);
                    set_node_weight(c, gap);
                    set_node_weight(v, w_i_i - w_e_i);
                    set_node_weight(w, w_i_i - w_i_e);
                } else {
                    if (path.size() <= 4) return;
                    NodeID d = path[path.size() - 2];
                    append_restore_data({c, d, v, w}, c, d, path[2], path[path.size() - 3]);
                    increment_reduction_offset(w_e_e + gap);
                    delete_path(path, 2, path.size() - 2);
                    relink(c, d, path[2], path[path.size() - 3]);
                    set_node_weight(c, -gap);
                    set_node_weight(d, -gap);
                    set_node_weight(v, w_i_e - w_e_e);
                    set_node_weight(w, w_e_i - w_e_e);
                }
            }
        }
    }

    bool next_vertex(std::stack<NodeID> &stack, NodeID &n) {
        for (;;) {
            if (stack.empty()) return false;
            n = stack.top();
            stack.pop();
            if (!marked[n]) return true;
        }
    }


    void find_MIS_on_path(NodeWeight &w_i_i, NodeWeight &w_i_e, NodeWeight &w_e_i, NodeWeight &w_e_e) {
        NodeWeight w_i = g.getNodeWeight(path[1]), w_e = 0;
        find_MIS_on_path(w_i, w_e, path);
        w_e_i = w_i;
        w_e_e = w_e;
        w_e = g.getNodeWeight(path[0]);
        w_i = 0;
        find_MIS_on_path(w_i, w_e, path);
        w_i_i = w_i;
        w_i_e = w_e;
    }

    template<bool track_choices=false>
    void find_MIS_on_path(NodeWeight &w_i, NodeWeight &w_e, std::vector<NodeID> &path) {
        for (size_t i = 2; i < path.size(); ++i) {
            if (track_choices)
                choices[i - 1] = w_i > w_e;
            NodeWeight next_e = std::max(w_e, w_i);
            w_i = w_e + g.getNodeWeight(path[i]);
            w_e = next_e;
        }
    }

    void find_max_path(NodeID n) {
        path.clear();
        path.push_back(n);

        NodeID last;
        NodeID current;
        NodeID next;
        forall_out_edges(g, e, n)
        last = n;
        current = g.getEdgeTarget(e);
        if (marked[current])
            continue;

        std::reverse(path.begin(), path.end());
        path.push_back(current);
        while (find_next_vertex_on_path(current, last, n, next)) {
            path.push_back(next);
            last = current;
            current = next;
        }
        if (path.front() == path.back())
            break;
        endfor
    }

    bool find_next_vertex_on_path(NodeID current, NodeID last, NodeID first, NodeID &next) {
        if (degrees[current] != 2 || current == first)
            return false;

        forall_out_edges(g, e, current)
        next = g.getEdgeTarget(e);
        if (!marked[next] && next != last)
            return true;
        endfor

        return false;
    }

    void append_restore_data(std::initializer_list<NodeID> weight_nodes, NodeID u=0, NodeID v=0, NodeID x=0, NodeID y=0) {
        restore_stack.emplace(path);
        restore_stack.top().unlink_nodes[0] = u;
        restore_stack.top().unlink_nodes[1] = v;
        restore_stack.top().unlink_nodes[2] = x;
        restore_stack.top().unlink_nodes[3] = y;
        for (NodeID n : weight_nodes) {
            restore_stack.top().weights.emplace_back(n, g.getNodeWeight(n));
        }
    }

    void try_append_vertex(NodeID n) {
        if (weight_diffs[n] >= 0)
            neighborhood_removal_vertices.push(n);
        else if (degrees[n] == 1)
            deg_1_vertices.push(n);
        else if (degrees[n] == 2) {
            deg_2_vertices.push(n);
        }
    }

    void delete_neighborhood(NodeID n) {
        forall_out_edges(g, e, n)
        delete_vertex(g.getEdgeTarget(e));
        endfor
    }

    void delete_path(std::vector<NodeID> &path, size_t from, size_t to) {
        for (size_t i = from + 1; i < to - 1; ++i) {
            marked[path[i]] = true;
        }
        edges -= 2 * std::max(0, static_cast<int>(to - from - 3));
        delete_vertex(path[from]);
        delete_vertex(path[to - 1]);
    }

    void delete_vertex(NodeID n) {
        if (marked[n])
            return;

        edges -= 2 * degrees[n];
        marked[n] = true;
        forall_out_edges(g, e, n)
        NodeID u = g.getEdgeTarget(e);
        if (marked[u]) continue;

        --degrees[u];
        weight_diffs[u] += g.getNodeWeight(n);
        try_append_vertex(u);
        endfor
    }

    //connect v and w; disconnect v from x and w from y
    void relink(NodeID v, NodeID w, NodeID x, NodeID y) {
        set_target(v, x, w);
        set_target(w, y, v);
        ++degrees[v];
        ++degrees[w];
        weight_diffs[v] -= g.getNodeWeight(w);
        weight_diffs[w] -= g.getNodeWeight(v);
        edges += 2;
    }

    //connect v and x, w and y, disconnect v from w
    void undo_relink(NodeID v, NodeID w, NodeID x, NodeID y) {
        set_target(v, w, x);
        set_target(w, v, y);
    }

    //disconnect n from orig_target, connect n to target
    void set_target(NodeID n, NodeID orig_target, NodeID target) {
        forall_out_edges(g, e, n)
        if (g.getEdgeTarget(e) == orig_target) {
            g.setEdgeTarget(e, target);
            return;
        }
        endfor
    }

    bool are_connected(NodeID v, NodeID w) const {
        forall_out_edges(g, e, v)
        if (g.getEdgeTarget(e) == w)
            return true;
        endfor
        return false;
    }

    void update_active_nodes() {
        for (size_t i = active_nodes.size(); i-- > 0;) {
            NodeID n = active_nodes[i];
            if (marked[n]) {
                active_nodes[i] = active_nodes.back();
                active_nodes.pop_back();
            }
        }
    }

    void calculate_IS() {
        while (restore_stack.size()) {
            auto &entry = restore_stack.top();
            for (auto &p : entry.weights)
                g.setNodeWeight(p.first, p.second);
            if (entry.unlink_nodes[0] != entry.unlink_nodes[1])
                undo_relink(entry.unlink_nodes[0], entry.unlink_nodes[1], entry.unlink_nodes[2], entry.unlink_nodes[3]);
            apply(entry.path);
            restore_stack.pop();
        }
    }

    void apply(std::vector<NodeID> &path) {
        NodeID v = path.front(), w = path.back();
        NodeWeight w_e = g.getPartitionIndex(v) ? g.getNodeWeight(v) : 0;
        NodeWeight w_i = g.getPartitionIndex(v) ? 0 : g.getNodeWeight(path[1]);
        find_MIS_on_path<true>(w_i, w_e, path);

        bool include = !g.getPartitionIndex(w) && choices[path.size() - 2];
        for (size_t i = path.size() - 2; i >= 1; --i) {
            g.setPartitionIndex(path[i], include);
            include = !include && choices[i - 1];
        }
    }

    void calc_reduced_graph() {
        std::vector<NodeID> map(g.number_of_nodes());
        NodeID id = 0;
        //size_t edges = 0;
        forall_nodes(g, n)
        map[n] = id;
        id += !marked[n];
        endfor
        reduced_graph.start_construction(id, edges);
        id = 0;
        forall_nodes(g, n)
        if (marked[n]) continue;

        reduced_graph.new_node();
        reduced_graph.setNodeWeight(id, g.getNodeWeight(n));
        forall_out_edges(g, e, n)
        NodeID v = g.getEdgeTarget(e);
        if (!marked[v]) {
            reduced_graph.new_edge(id, map[v]);
        }
        endfor
        ++id;
        endfor
        reduced_graph.finish_construction();
    }
};

}
#endif //LINEARTIMEKERNELIZATION_REDUCER_H

