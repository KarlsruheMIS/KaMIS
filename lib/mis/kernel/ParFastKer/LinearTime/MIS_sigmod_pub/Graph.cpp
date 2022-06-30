 /******************************************************************************
 * Copyright (C) 2019 Lijun Chang <ljchang.au@gmail.com>
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


// #include "Utility.h"
#include "mis/kernel/ParFastKer/LinearTime/MIS_sigmod_pub/Utility.h"
#include "mis/kernel/ParFastKer/LinearTime/MIS_sigmod_pub/Graph.h"
// #include "Graph.h"
#include <fstream>
#include <limits>
#include <algorithm>

using namespace std;

namespace LinearTime {
Graph::Graph() {

	n = m = 0;

	pstart = NULL;
	edges = NULL;
    reverse_mapping = NULL;
}

Graph::~Graph() {
	if(pstart != NULL) {
		delete[] pstart;
		pstart = NULL;
	}
	if(edges != NULL) {
		delete[] edges;
		edges = NULL;
	}
	if(reverse_mapping != NULL) {
		delete[] reverse_mapping;
		reverse_mapping = NULL;
	}
}

void Graph::write_kernel_to_adj_list(const char *is, const int *degree, const ui *pend, std::vector<std::vector<int>> &out_kernel) {
    reverse_mapping = new ui[n];
    ui *mapping = new ui[n];
    int kernal_size = 0;
    // int kernal_edges = 0;
    for(ui k = 0;k < n;k ++) {
        if(is[k]&&degree[k] > 0) {
            mapping[k] = kernal_size;
            ++ kernal_size;
            // for(ui j = pstart[k];j < pend[k];j ++) if(is[edges[j]]) ++ kernal_edges;
        }
    }

    out_kernel.resize(kernal_size);
    for(ui k = 0;k < n;k ++) {
        if(is[k]&&degree[k] > 0) {
            reverse_mapping[mapping[k]] = k;
            for(ui j = pstart[k];j < pend[k];j ++) {
                if(is[edges[j]]) {
                    out_kernel[mapping[k]].push_back(mapping[edges[j]]);
                }
            }
        }
    }
    delete[] mapping;

}

void Graph::write_graph_metis(const char *is, const int *degree, const ui *pend) {
    ui *mapping = new ui[n];
    int kernal_size = 0;
    int kernal_edges = 0;
    for(ui k = 0;k < n;k ++) {
        if(is[k]&&degree[k] > 0) {
            mapping[k] = kernal_size;
            ++ kernal_size;
            for(ui j = pstart[k];j < pend[k];j ++) if(is[edges[j]]) ++ kernal_edges;
        }
    }

    std::ofstream outputfile ("kernel.graph");

    outputfile << kernal_size << " " << kernal_edges / 2 << std::endl;

    for(ui k = 0;k < n;k ++) {
        if(is[k]&&degree[k] > 0) {
            for(ui j = pstart[k];j < pend[k];j ++) {
                if(is[edges[j]]) {
                    outputfile << mapping[edges[j]] + 1 << " ";
                }
            }
            outputfile << std::endl;
        }
    }
    outputfile.close();
    delete[] mapping;

}

void Graph::inputGraphAdjList(std::vector<std::vector<int>> &adj) {

    n = adj.size();
    m = 0;
    for(auto neighborarray : adj) {
        m += neighborarray.size();
    }
    if (pstart == NULL)
        pstart = new ui[n + 1];
    if (edges == NULL)
        edges = new ui[m];

    pstart[0] = 0;
    for(ui i = 0; i < adj.size(); ++i) {
        ui j = 0;
        for (j = 0; j < adj[i].size(); ++j) {
            edges[pstart[i] + j] = adj[i][j];
        }

        pstart[i + 1] = pstart[i] + j;
    }
}

void Graph::inputEdgeList(std::vector<ui> &_pstart, std::vector<ui> &_edges) {
    n = _pstart.size() - 1;
    m = _edges.size();

    if (pstart == NULL)
        pstart = new ui[n + 1];
    if (edges == NULL)
        edges = new ui[m];


    for(ui i = 0; i < n+1; ++i) {
        pstart[i] = _pstart[i];
    }
    for(ui j = 0; j < m; ++j) {
        edges[j] = _edges[j];
    }
}
void Graph::read_graph_metis() {
  std::string line;

  // open file for reading
  std::ifstream in(dir.c_str());
  if (!in) {
    std::cerr << "Error opening " << dir.c_str() << std::endl;
    exit(1);
  }

  std::getline(in, line);
  // skip comments
  while (line[0] == '%') {
    std::getline(in, line);
  }

  std::stringstream ss(line);
  ss >> n;
  ss >> m;

#ifdef UNSAFE_LONG
  if (2 * m > std::numeric_limits<long>::max() ||
      n > std::numeric_limits<long>::max()) {
    std::cerr << "The graph is too large. Currently only 64bit supported!"
              << std::endl;
    exit(0);
  }
#else
  if (2 * m > std::numeric_limits<int>::max() ||
      n > std::numeric_limits<int>::max()) {
    std::cerr << "The graph is too large. Currently only 32bit supported!"
              << std::endl;
    exit(0);
  }
#endif // UNSAFE_LONG

  m *= 2; // since we have forward and backward edges

  printf("\tn = %u; m = %u (undirected)\n", n, m);

  if (pstart == NULL)
    pstart = new ui[n + 1];
  if (edges == NULL)
    edges = new ui[m];

  ui i = 0;
  pstart[0] = 0;
  while (std::getline(in, line)) {
    if (line[0] == '%')
      continue;

    std::stringstream ss(line);

    ui target;
    ui j = 0;
    while (ss >> target) {
      edges[pstart[i] + j] = target - 1;
      j++;
    }

    pstart[i + 1] = pstart[i] + j;
    i++;

    if (in.eof())
      break;
  }
  in.close();
}

void Graph::read_graph() {
	FILE *f = open_file((dir + string("/b_degree.bin")).c_str(), "rb");

	int tt;
	fread(&tt, sizeof(int), 1, f);
	if(tt != (int)sizeof(int)) {
		printf("sizeof int is different: edge.bin(%d), machine(%d)\n", tt, (int)sizeof(int));
		return ;
	}
	fread(&n, sizeof(int), 1, f);
	fread(&m, sizeof(int), 1, f);

	printf("\tn = %u; m = %u (undirected)\n", n, m);

	ui *degree = new ui[n];
	fread(degree, sizeof(int), n, f);

#ifndef NDEBUG
	long long sum = 0;
	for(ui i = 0;i < n;i ++) sum += degree[i];
	if(sum != m) printf("WA input graph\n");
#endif

	fclose(f);

	f = open_file((dir + string("/b_adj.bin")).c_str(), "rb");

	if(pstart == NULL) pstart = new ui[n+1];
	if(edges == NULL) edges = new ui[m];

	ui *buf = new ui[n];

	pstart[0] = 0;
	for(ui i = 0;i < n;i ++) {
		if(degree[i] > 0) fread(buf, sizeof(int), degree[i], f);

		for(ui j = 0;j < degree[i];j ++) {
			edges[pstart[i] + j] = buf[j];

			if(buf[j] == i) printf("Self-loop!\n");
		}

		pstart[i+1] = pstart[i] + degree[i];
	}

	delete[] buf;
	buf = NULL;

	fclose(f);

	delete[] degree;
}

void Graph::check_is(const char *is, int count) {
	int cnt = 0;
	for(ui i = 0;i < n;i ++) if(is[i]) ++ cnt;
	if(count != cnt) printf("WA count in check_is! %d\n", cnt);

	int maximal = 1;
	for(ui i = 0;i < n;i ++) {
		if(is[i]) {
			for(ui j = pstart[i];j < pstart[i+1];j ++) if(is[edges[j]]) {
				printf("WA conflict is in check_is!\n");
			}
		}
		else if(maximal) {
			int find = 0;
			for(ui j = pstart[i];j < pstart[i+1];j ++) if(is[edges[j]]) {
				find = 1;
				break;
			}
			if(!find) {
				maximal = 0;
				//printf("%d:", i);
				//for(ui j = pstart[i];j < pstart[i+1];j ++) printf(" %d(%d),", edges[j], is[edges[j]]);
			}
		}
	}
	if(!maximal) printf("**** WA not maximal!\n");
}

void Graph::compute_upperbound(const char *is, char *fixed) {
	int delete_fixed = 0;
	if(fixed == NULL) {
		delete_fixed = 1;
		fixed = new char[n];
		memset(fixed, 0, sizeof(char)*n);
	}
	ui *degree = new ui[n];
	vector<pair<ui,ui> > vp;
	for(ui i = 0;i < n;i ++) if(is[i]&&!fixed[i]) {
		degree[i] = 0;
		for(ui j = pstart[i];j < pstart[i+1];j ++) if(!fixed[edges[j]]) ++ degree[i];
		vp.pb(mp(degree[i], i));
	}
	sort(vp.begin(), vp.end());

	ui *ids = new ui[n];
	ui *reverse_ids = new ui[n];
	ui *degree_starts = new ui[n];

	ui i = 0, j = 0;
	while(j < vp.size()) {
		degree_starts[i] = j;
		while(j < vp.size()&&vp[j].first == i) {
			ids[j] = vp[j].second;
			reverse_ids[vp[j].second] = j;
			++ j;
		}
		++ i;
	}
	degree_starts[i] = j;

	char *vis = new char[n];
	memset(vis, 0, sizeof(char)*n);

	int res = 0;
	for(ui i = 0;i < n;i ++) {
		if(is[i]&&fixed[i]) ++ res;
		if(fixed[i]) vis[i] = 1;
	}

	for(ui i = 0;i < vp.size();i ++) {
		ui u = ids[i];
		if(degree[u] > 0) res += degree[u];
		else ++ res;

		++ degree_starts[degree[u]];
		vis[u] = 1;

		for(ui j = pstart[u];j < pstart[u+1];j ++) if(!vis[edges[j]]){
			ui v = edges[j];
			vis[v] = 1;

			for(ui k = pstart[v];k < pstart[v+1];k ++) if(!vis[edges[k]]&&is[edges[k]]) {
				ui w = edges[k];
				ui deg = degree[w];

#ifndef NDEBUG
				if(ids[reverse_ids[w]] != w) printf("WWA!\n");
				if(reverse_ids[w] < degree_starts[deg]||reverse_ids[w] >= degree_starts[deg+1]) {
					printf("%d %d\n", i, deg);
					printf("WA in compute_upperbound! %d %d %d\n", degree_starts[deg], reverse_ids[w], degree_starts[deg+1]);
				}
				if(degree_starts[deg] <= i) printf("WWAA!\n");
#endif
				reverse_ids[ids[degree_starts[deg]]] = reverse_ids[w];
				swap(ids[degree_starts[deg]], ids[reverse_ids[w]]);
				reverse_ids[w] = degree_starts[deg];

				++ degree_starts[deg];
				-- degree[w];

				if(degree_starts[degree[w]] <= i) degree_starts[degree[w]] = i+1;
			}
		}
	}

	memset(vis, 0, sizeof(char)*n);
	int res2 = 0;

	for(ui i = 0;i < n;i ++) {
		if(is[i]&&fixed[i]) ++ res2;
		if(fixed[i]) vis[i] = 1;
	}

	for(ui i = 0;i < vp.size();i ++) degree[vp[i].second] = vp[i].first;

	for(ui i = 0;i < vp.size();i ++) {
		ui u = vp[i].second;
		res2 += degree[u];
		if(degree[u] == 0) ++ res2;

		for(ui j = pstart[u];j < pstart[u+1];j ++) if(!vis[edges[j]]) {
			ui v = edges[j];
			vis[v] = 1;
			for(ui k = pstart[v];k < pstart[v+1];k ++) if(is[edges[k]]) -- degree[edges[k]];
		}
	}

	printf("\tupper bound: %d(static) %d(dynamic)\n", res2, res);

	delete[] vis;
	delete[] degree;
	delete[] ids;
	delete[] reverse_ids;
	delete[] degree_starts;

	if(delete_fixed) delete[] fixed;
}

int Graph::general_swap(char *is, char *fixed) {
	int delete_fixed = 0;
	if(fixed == NULL) {
		fixed = new char[n];
		delete_fixed = 1;
		memset(fixed, 0, sizeof(char)*n);
	}

	vector<ui> nonfixed;
	for(ui i = 0;i < n;i ++) if(!fixed[i]) nonfixed.pb(i);

	ui *pend = new ui[n];
	for(ui i = 0;i < nonfixed.size();i ++) {
		ui u = nonfixed[i];
		pend[u] = pstart[u+1];
		ui j = pstart[u];
		while(true) {
			while(j < pend[u]&&!fixed[edges[j]]) ++ j;
			while(j < pend[u]&&fixed[edges[pend[u]-1]]) -- pend[u];

			if(j >= pend[u]) break;
			swap(edges[j], edges[pend[u]-1]);
		}
		sort(edges+pstart[u], edges+pend[u]);

#ifndef NDEBUG
		for(j = pend[u];j < pstart[u+1];j ++) if(!fixed[edges[j]]) printf("WA shrinking!\n");
#endif
	}

	int res = 0;
	for(ui i = 0;i < n;i ++) if(is[i]) ++ res;
	int old_res = res, increase = 0;

	int *isn = new int[n];
	for(ui i = 0;i < nonfixed.size();i ++) if(!is[nonfixed[i]]) {
		ui u = nonfixed[i];
		isn[u] = 0;
		for(ui j = pstart[u];j < pend[u];j ++) if(is[edges[j]]) ++ isn[u];
#ifndef NDEBUG
		if(isn[u] == 0) {
			for(ui j = pstart[u];j < pstart[u+1];j ++) {
				printf("(%d,is:%d,fixed:%d) ", edges[j], is[edges[j]], fixed[edges[j]]);
			}
			printf("WA in general_swap: not maximal! %d\n", u);
		}
#endif
	}

	int first_time = 1, two_three_swap = 1;
	while(two_three_swap) {
		vector<ui> candidates;
		for(ui i = 0;i < nonfixed.size();i ++) if(is[nonfixed[i]]) {
			ui u = nonfixed[i];
			int cnt = 0;
			for(ui j = pstart[u];j < pend[u];j ++) if(!is[edges[j]]&&isn[edges[j]] == 1) ++ cnt;
			if(cnt >= 2) candidates.pb(u);
		}

		while(!candidates.empty()) {
			ui u = candidates.back();
			candidates.pop_back();
			if(!is[u]) continue;

			vector<ui> neighbors;
			for(ui i = pstart[u];i < pend[u];i ++) if(!is[edges[i]]&&isn[edges[i]] == 1) neighbors.pb(edges[i]);

#ifndef NDEBUG
			for(ui i = 0;i < neighbors.size();i ++) {
				ui u = neighbors[i];
				if(fixed[u]) printf("WA1\n");
				if(pend[u] > pstart[u+1]) printf("WA2\n");
				for(ui j = pstart[u];j < pend[u];j ++) if(edges[j] >= n||fixed[edges[j]]) printf("WA3\n");
			}
#endif

			for(ui i = 0;i + 1 < neighbors.size();i ++) {
				ui v = neighbors[i];
				int find = 0;

				ui k = pstart[v];
				for(ui j = i+1;j < neighbors.size();j ++) {
					while(k < pend[v]&&edges[k] < neighbors[j]) ++ k;
					if(k >= pend[v]||edges[k] > neighbors[j]) {
						find = 1;
						break;
					}
				}

				if(find) {
					is[u] = 0; isn[u] = 0; -- res;
					for(ui j = pstart[u];j < pend[u];j ++) {
#ifndef NDEBUG
						if(is[edges[j]]) printf("conflict is!\n");
#endif
						-- isn[edges[j]];
						if(isn[edges[j]] == 1) {
							for(ui k = pstart[edges[j]];k < pend[edges[j]];k ++) if(is[edges[k]]) {
								candidates.pb(edges[k]);
							}
						}
					}

					for(ui j = i;j < neighbors.size();j ++) {
						v = neighbors[j];
						find = 0;
						for(ui k = pstart[v];k < pend[v];k ++) if(is[edges[k]]) {
							find = 1;
							break;
						}
						if(!find) {
							is[v] = 1; ++ res;
							candidates.pb(v);
							for(ui k = pstart[v];k < pend[v];k ++) ++ isn[edges[k]];
						}
					}
					break;
				}
			}
		}

		if(first_time) {
			increase = res-old_res;
			//printf(" %d", increase);
			first_time = 0;
		}

		two_three_swap = 0;
		for(ui i = 0;i < nonfixed.size();i ++) if(!is[nonfixed[i]]&&isn[nonfixed[i]] == 2) {
			ui u = nonfixed[i];
			ui x = n, y = n;
			for(ui j = pstart[u];j < pend[u];j ++) if(is[edges[j]]) {
				if(x == n) x = edges[j];
				else if(y == n) y = edges[j];
				else printf("WA isn[u] = 2!\n");
			}

			vector<ui> nx, ny;
			ui k = pstart[u];
			for(ui j = pstart[x];j < pend[x];j ++) if(isn[edges[j]] <= 2&&edges[j] != u) {
				while(k < pend[u]&&edges[k] < edges[j]) ++ k;
				if(k >= pend[u]||edges[k] > edges[j]) nx.pb(edges[j]);
			}
			k = pstart[u];
			for(ui j = pstart[y];j < pend[y];j ++) if(isn[edges[j]] <= 2&&edges[j] != u) {
				while(k < pend[u]&&edges[k] < edges[j]) ++ k;
				if(k >= pend[u]||edges[k] > edges[j]) ny.pb(edges[j]);
			}

			vector<ui> neighbors;
			ui j = 0;
			k = 0;
			while(j < nx.size()&&k < ny.size()) {
				if(nx[j] == ny[k]) {
					neighbors.pb(nx[j]);
					++ j; ++ k;
				}
				else if(nx[j] < ny[k]) {
					if(isn[nx[j]] == 1) neighbors.pb(nx[j]);
					++ j;
				}
				else {
					if(isn[ny[k]] == 1) neighbors.pb(ny[k]);
					++ k;
				}
			}
			while(j < nx.size()) {
				if(isn[nx[j]] == 1) neighbors.pb(nx[j]);
				++ j;
			}
			while(k < ny.size()) {
				if(isn[ny[k]] == 1) neighbors.pb(ny[k]);
				++ k;
			}

			for(ui ii = 0;ii + 1 < neighbors.size();ii ++) {
				ui v = neighbors[ii];
				int find = 0;

				ui k = pstart[v];
				for(ui j = ii+1;j < neighbors.size();j ++) {
					while(k < pend[v]&&edges[k] < neighbors[j]) ++ k;
					if(k >= pend[v]||edges[k] > neighbors[j]) {
						find = 1;
						break;
					}
				}
				if(find) {
					two_three_swap = 1;
					is[x] = 0; is[y] = 0;
					isn[x] = 0; isn[y] = 0;
					for(ui j = pstart[x];j < pend[x];j ++) -- isn[edges[j]];
					for(ui j = pstart[y];j < pend[y];j ++) -- isn[edges[j]];
					-- res;

					is[u] = 1;
					for(ui k = pstart[u];k < pend[u];k ++) ++ isn[edges[k]];

					for(ui j = ii;j < neighbors.size();j ++) {
						v = neighbors[j];
						find = 0;
						for(ui k = pstart[v];k < pend[v];k ++) if(is[edges[k]]) {
							find = 1;
							break;
						}
						if(!find) {
							is[v] = 1;
							++ res;
							for(ui k = pstart[v];k < pend[v];k ++) ++ isn[edges[k]];
						}
					}
					break;
				}
			}
		}
	}

	//printf(" %d\tMIS after general_swap %d\n", res-old_res-increase, res);
	printf("Improve1 (%d),", res-old_res);

	if(delete_fixed) delete[] fixed;
	delete[] isn;
	delete[] pend;

	return res;
}

void Graph::construct_degree_increase(ui *ids) {
	memset(ids, 0, sizeof(ui)*n);

	for(ui i = 0;i < n;i ++) {
		ui d = pstart[i+1] - pstart[i];
		++ ids[d];
	}
	for(ui i = 1;i < n;i ++) ids[i] += ids[i-1];

	ui *order = new ui[n];
	for(ui i = 0;i < n;i ++) {
		ui d = pstart[i+1] - pstart[i];
		order[i] = ids[d];
		-- ids[d];
	}

	for(ui i = 0;i < n;i ++) ids[order[i]-1] = i;

	delete[] order;
}

void Graph::greedy() {
#ifdef _LINUX_
	struct timeval start, end1, end;
	gettimeofday(&start, NULL);
#endif

	char *is = new char[n];
	for(ui i = 0;i < n;i ++) is[i] = 1;

	ui *ids = new ui[n];
	construct_degree_increase(ids);

	int res = 0;
	for(ui i = 0;i < n;i ++) if(is[ids[i]]) {
		ui u = ids[i];
		++ res;
		for(ui j = pstart[u];j < pstart[u+1];j ++) is[edges[j]] = 0;
	}

	printf("Greedy MIS: %d\n", res);

	delete[] ids;

#ifdef _LINUX_
	gettimeofday(&end1, NULL);

	long long mtime1, seconds1, useconds1;
	seconds1 = end1.tv_sec - start.tv_sec;
	useconds1 = end1.tv_usec - start.tv_usec;
	mtime1 = seconds1*1000000 + useconds1;
#endif

#ifndef NDEBUG
	check_is(is, res);
#endif

	//ils local_search(n, pstart, edges, is);
	//res = local_search.perform_ils(mtime1*100);

//#ifndef NDEBUG
//	check_is(is, res);
//#endif

	delete[] is;

#ifdef _LINUX_
	gettimeofday(&end, NULL);

	long long mtime, seconds, useconds;
	seconds = end.tv_sec - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = seconds*1000000 + useconds;

	printf("Process time: %lld, Swap time: %lld, Total time: %lld\n", mtime1, mtime-mtime1, mtime);
#endif
}

void Graph::greedy_dynamic() {
#ifdef _LINUX_
	struct timeval start, end1, end;
	gettimeofday(&start, NULL);
#endif

	char *is = new char[n];
	for(ui i = 0;i < n;i ++) is[i] = 1;

	ui *ids = new ui[n];
	construct_degree_increase(ids);

	ui *degree = new ui[n];
	for(ui i = 0;i < n;i ++) degree[i] = pstart[i+1]-pstart[i];

	ui *degree_starts = new ui[n+1];
	ui *reverse_ids = new ui[n];

	ui i = 0, j = 0;
	while(j < n) {
		degree_starts[i] = j;
		while(j < n&&degree[ids[j]] == i) {
			reverse_ids[ids[j]] = j;
			++ j;
		}
		++ i;
	}
	degree_starts[i] = j;

	int res = 0;
	for(ui ii = 0;ii < n;ii ++) {
		ui u = ids[ii];
		degree_starts[degree[u]] = ii+1;
		if(!is[u]) continue;

#ifndef NDEBUG
		if(reverse_ids[u] != ii) printf("WA1 in greedy_dynamic!\n");
#endif

		++ res;

		for(ui j = pstart[u];j < pstart[u+1];j ++) if(is[edges[j]]) {
			ui v = edges[j];
			is[v] = 0;
			for(ui k = pstart[v];k < pstart[v+1];k ++) if(is[edges[k]]&&edges[k] != u) {
				ui w = edges[k];

#ifndef NDEBUG
				ui deg = degree[w];
				if(ids[reverse_ids[w]] != w) printf("WWA!\n");
				if(reverse_ids[w] < degree_starts[deg]||reverse_ids[w] >= degree_starts[deg+1]) {
					printf("%d %d\n", ii, deg);
					printf("WA in greedy dynamic! %d %d %d\n", degree_starts[deg], reverse_ids[w], degree_starts[deg+1]);
				}
				if(degree_starts[deg] <= ii) printf("WWAA!\n");
#endif
				ui ds = degree_starts[degree[w]];
				reverse_ids[ids[ds]] = reverse_ids[w];
				swap(ids[ds], ids[reverse_ids[w]]);
				reverse_ids[w] = ds;

				++ degree_starts[degree[w]];
				-- degree[w];

				if(degree_starts[degree[w]] <= ii) degree_starts[degree[w]] = ii+1;
			}
		}
	}
	printf("greedy_dynamic MIS: %d\n", res);

	delete[] ids;
	delete[] degree_starts;
	delete[] reverse_ids;
	delete[] degree;

#ifdef _LINUX_
	gettimeofday(&end1, NULL);

	long long mtime1, seconds1, useconds1;
	seconds1 = end1.tv_sec - start.tv_sec;
	useconds1 = end1.tv_usec - start.tv_usec;
	mtime1 = seconds1*1000000 + useconds1;
#endif

#ifndef NDEBUG
	check_is(is, res);
#endif

	//ils local_search(n, pstart, edges, is);
	//res = local_search.perform_ils(mtime1*factor);

#ifndef NDEBUG
	check_is(is, res);
#endif

	delete[] is;

#ifdef _LINUX_
	gettimeofday(&end, NULL);

	long long mtime, seconds, useconds;
	seconds = end.tv_sec - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = seconds*1000000 + useconds;

	printf("Process time: %lld, Swap time: %lld, Total time: %lld\n", mtime1, mtime-mtime1, mtime);
#endif
}

void Graph::degree_one_kernal_and_remove_max_degree() {
#ifdef _LINUX_
	struct timeval start, end1, end;
	gettimeofday(&start, NULL);
#endif

	char *is = new char[n];
	for(ui i = 0;i < n;i ++) is[i] = 1;

	int *bin_head = new int[n];
	int *bin_next = new int[n];
	int *degree = new int[n];
	memset(bin_head, -1, sizeof(int)*n);

	vector<ui> degree_ones, S;

	int max_d = 0, res = 0;
	for(ui i = 0;i < n;i ++) {
		degree[i] = pstart[i+1] - pstart[i];
		bin_next[i] = bin_head[degree[i]];
		bin_head[degree[i]] = i;

		if(degree[i] == 0) ++ res;
		if(degree[i] == 1) degree_ones.pb(i);
		if(degree[i] > max_d) max_d = degree[i];
	}

	char *fixed = new char[n];
	memset(fixed, 0, sizeof(char)*n);

	int kernal_size = 0, first_time = 1;
	while(!degree_ones.empty()||max_d >= 2) {
		while(!degree_ones.empty()) {
			ui u = degree_ones.back();
			degree_ones.pop_back();
			if(!is[u]||degree[u] != 1) continue;

			int cnt = 0;
			for(int j = pstart[u];j < pstart[u+1];j ++) if(is[edges[j]]) {
				++ cnt;
				res += delete_vertex(edges[j], is, degree, degree_ones);
			}
#ifndef NDEBUG
			if(cnt > 1) printf("WA degree one! %d\n", degree[u]);
#endif
		}

		if(first_time) {
			first_time = 0;
			for(ui k = 0;k < n;k ++) {
				if(is[k]&&degree[k] > 0) ++ kernal_size;
				else fixed[k] = 1;
			}
		}

		while(degree_ones.empty()) {
			while(max_d >= 2&&bin_head[max_d] == -1) -- max_d;
			if(max_d < 2) break;

			int v = -1;
			for(v = bin_head[max_d];v != -1;) {
				int tmp = bin_next[v];
				if(is[v]&&degree[v] > 0) {
					if(degree[v] < max_d) {
						bin_next[v] = bin_head[degree[v]];
						bin_head[degree[v]] = v;
					}
					else {
						S.pb(v);
						res += delete_vertex(v, is, degree, degree_ones);

						bin_head[max_d] = tmp;
						break;
					}
				}
				v = tmp;
			}
			if(v == -1) bin_head[max_d] = -1;
		}
	}

	for(int i = S.size()-1;i >= 0;i --) {
		ui u = S[i];

		int ok = 1;
		for(ui i = pstart[u];i < pstart[u+1];i ++) if(is[edges[i]]) {
			ok = 0;
			break;
		}
		if(ok) {
			is[u] = 1;
			++ res;
		}
	}

	//printf("Degree_one MIS: %d (kernal |V|: %d, inexact reduction: %lu)\n", res, kernal_size, S.size());

	delete[] bin_head;
	delete[] bin_next;
	delete[] degree;

#ifdef _LINUX_
	gettimeofday(&end1, NULL);

	long long mtime1, seconds1, useconds1;
	seconds1 = end1.tv_sec - start.tv_sec;
	useconds1 = end1.tv_usec - start.tv_usec;
	mtime1 = seconds1*1000000 + useconds1;
#endif

#ifndef NDEBUG
	check_is(is, res);
#endif

	//if(!S.empty()) general_swap(is, fixed);

	delete[] is;
	delete[] fixed;

#ifdef _LINUX_
	gettimeofday(&end, NULL);

	long long mtime, seconds, useconds;
	seconds = end.tv_sec - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = seconds*1000000 + useconds;

	printf("Process time: %lld, Swap time: %lld, Total time: %lld\n", mtime1, mtime-mtime1, mtime);
#endif
}

int Graph::delete_vertex(ui v, char *is, int *degree, vector<ui> &degree_ones, vector<ui> &degree_twos) {
	is[v] = 0;
	int res = 0;
	for(int k = pstart[v];k < pstart[v+1];k ++) if(is[edges[k]]) {
		int w = edges[k];
		-- degree[w];
		if(degree[w] == 0) ++ res;
		else if(degree[w] == 1) degree_ones.pb(w);
		else if(degree[w] == 2) degree_twos.pb(w);
	}
	return res;
}

/*
void Graph::delete_vertex(ui v, char *is, int *degree, vector<ui> &degree_ones, vector<ui> &degree_twos, vector<ui> &degree_threes) {
	// printf("delete %d\n", v);
	is[v] = 0;
	for(int k = pstart[v];k < pstart[v+1];k ++) if(edges[k] != n&&is[edges[k]]) {
		int w = edges[k];
		-- degree[w];
		if(degree[w] == 1) degree_ones.pb(w);
		else if(degree[w] == 2) degree_twos.pb(w);
		else if(degree[w] == 3) degree_threes.pb(w);
	}
}
*/

int Graph::delete_vertex(ui v, const ui *pend, char *is, int *degree, vector<ui> &degree_ones, vector<ui> &degree_twos) {
	is[v] = 0;
	int res = 0;
	for(int k = pstart[v];k < pend[v];k ++) if(is[edges[k]]) {
		int w = edges[k];
		-- degree[w];
		if(degree[w] == 0) ++ res;
		else if(degree[w] == 1) degree_ones.pb(w);
		else if(degree[w] == 2) degree_twos.pb(w);
	}
	return res;
}

int Graph::delete_vertex(ui v, char *is, int *degree, vector<ui> &degree_ones) {
	is[v] = 0;
	int res = 0;
	for(int k = pstart[v];k < pstart[v+1];k ++) if(is[edges[k]]) {
		int w = edges[k];
		-- degree[w];
		if(degree[w] == 0) ++ res;
		else if(degree[w] == 1) degree_ones.pb(w);
	}
	return res;
}

int Graph::exist_edge(ui u, ui v, const ui *pend) {
	if(pend[u]-pstart[u] < pend[v]-pstart[v]) {
		for(ui i = pstart[u];i < pend[u];i ++) {
			if(edges[i] == v) return 1;
		}
		return 0;
	}
	for(ui i = pstart[v];i < pend[v];i ++) {
		if(edges[i] == u) return 1;
	}
	return 0;
}

int Graph::find_other_endpoint(ui u, ui v, char *is) {
	for(ui j = pstart[u];j < pstart[u+1];j ++) {
		if(edges[j] != v&&is[edges[j]]) return edges[j];
	}
#ifndef NDEBUG
	int cnt = 0;
	for(ui j = pstart[u];j < pstart[u+1];j ++) {
		if(edges[j] != n&&is[edges[j]]) {
			printf(" %d,", edges[j]);
			++ cnt;
		}
	}
	printf("WA in find_other_endpoint! degree(%d), exclude(%d)\n", cnt, v);
#endif
	return -1;
}

ui Graph::edge_rewire(ui u, const ui *pend, ui v, ui w) {
#ifndef NDEBUG
	for(ui i = pstart[u];i < pend[u];i ++) if(edges[i] == w) printf("WA preexist edge!\n");
#endif

	for(ui i = pstart[u];i < pend[u];i ++) if(edges[i] == v) {
		edges[i] = w;
		return i;
	}
	printf("WA in edge_rewire!\n");
	return 0;
}

void Graph::shrink(ui u, ui &end, const char *is) {
	ui i = pstart[u];
	while(true) {
		while(i < end&&is[edges[i]]) ++ i;
		while(i < end&&!is[edges[end-1]]) -- end;

		if(i >= end) break;
		swap(edges[i], edges[end-1]);
	}
}

void Graph::shrink(ui u, ui &end, const char *is, ui *tri) {
	ui i = pstart[u];
	while(true) {
		while(i < end&&is[edges[i]]) ++ i;
		while(i < end&&!is[edges[end-1]]) -- end;

		if(i >= end) break;
		swap(edges[i], edges[end-1]);
		swap(tri[i], tri[end-1]);
	}
}

int Graph::degree_two_kernal_and_remove_max_degree_without_contraction(std::vector<std::vector<int>> &out_kernel) {
#ifndef NDEBUG
	ui *tmp_edges = new ui[m];
	for(ui i = 0;i < m;i ++) tmp_edges[i] = edges[i];
	ui *tmp_start = new ui[n+1];
	for(ui i = 0;i <= n;i ++) tmp_start[i] = pstart[i];
#endif

#ifdef _LINUX_
	struct timeval start, end1, end;
	gettimeofday(&start, NULL);
#endif

	S = vector<pair<ui,ui> >(0) ;
	is = new char[n];
	for(ui i = 0;i < n;i ++) is[i] = 1;

	int *bin_head = new int[n];
	int *bin_next = new int[n];
	int *degree = new int[n];
	memset(bin_head, -1, sizeof(int)*n);

	vector<ui> degree_ones, degree_twos;
	vector<pair<pair<ui,ui>, ui> > modified_edges;

	int max_d = 0, res = 0;
	for(ui i = 0;i < n;i ++) {
		degree[i] = pstart[i+1] - pstart[i];
		bin_next[i] = bin_head[degree[i]];
		bin_head[degree[i]] = i;

		if(degree[i] == 0) ++ res;
		else if(degree[i] == 1) degree_ones.pb(i);
		else if(degree[i] == 2) degree_twos.pb(i);

		if(degree[i] > max_d) max_d = degree[i];
	}

	fixed = new char[n];
	memset(fixed, 0, sizeof(char)*n);

	ui *pend = new ui[n];
	for(ui i = 0;i < n;i ++) pend[i] = pstart[i+1];

	int kernal_size = 0, inexact = 0, first_time = 1, S_size = (int)S.size();
	int kernal_edges = 0;
	// while(!degree_ones.empty()||!degree_twos.empty()) {
		while(!degree_ones.empty()||!degree_twos.empty()) {
			while(!degree_ones.empty()) {
				ui u = degree_ones.back();
				degree_ones.pop_back();
				if(!is[u]||degree[u] != 1) continue;

				int cnt = 0;
				for(int j = pstart[u];j < pend[u];j ++) if(is[edges[j]]) {
					++ cnt;
					res += delete_vertex(edges[j], pend, is, degree, degree_ones, degree_twos);
				}
#ifndef NDEBUG
				if(cnt > 1) printf("WA degree two kernel! %d\n", degree[u]);
#endif
			}

			while(!degree_twos.empty()&&degree_ones.empty()) {
				ui u = degree_twos.back();
				degree_twos.pop_back();
				if(!is[u]||degree[u] != 2) continue;

				shrink(u, pend[u], is);
#ifndef NDEBUG
				if(pend[u] != pstart[u] + 2) printf("WA degree two! %d\n", u);
#endif
				ui u1 = edges[pstart[u]], u2 = edges[pstart[u]+1];

				ui pre = u, cnt = 1;
				while(u1 != u&&degree[u1] == 2) {
					++ cnt;
					shrink(u1, pend[u1], is);
#ifndef NDEBUG
					if(pend[u1] != pstart[u1] + 2) printf("WA degree two! %d\n", u1);
#endif
					int tmp = u1;
					if(edges[pstart[u1]] != pre) u1 = edges[pstart[u1]];
					else u1 = edges[pstart[u1]+1];
					pre = tmp;
				}
#ifndef NDEBUG
				if(u1 != u&&degree[u1] <= 2) printf("WAxx!\n");
#endif
				if(u1 == u) {
					res += delete_vertex(u, pend, is, degree, degree_ones, degree_twos);
#ifndef NDEBUG
					if(degree_ones.empty()) printf("WWAAS\n");
#endif
					continue;
				}

				pre = u;
				while(degree[u2] == 2) {
					++ cnt;
					shrink(u2, pend[u2], is);
#ifndef NDEBUG
					if(pend[u2] != pstart[u2] + 2) printf("WA degree two! %d\n", u2);
#endif
					int tmp = u2;
					if(edges[pstart[u2]] != pre) u2 = edges[pstart[u2]];
					else u2 = edges[pstart[u2]+1];
					pre = tmp;
				}
				if(u1 == u2) {
					res += delete_vertex(u1, pend, is, degree, degree_ones, degree_twos);
					continue;
				}

				shrink(u1, pend[u1], is);
				shrink(u2, pend[u2], is);

				if(cnt%2 == 1) {
					if(exist_edge(u1, u2, pend)) {
						res += delete_vertex(u1, pend, is, degree, degree_ones, degree_twos);
						res += delete_vertex(u2, pend, is, degree, degree_ones, degree_twos);
#ifndef NDEBUG
						if(degree_ones.empty()) printf("WWAAS2\n");
#endif
					}
					else if(cnt > 1) {
						ui idx = pstart[pre];
						if(edges[idx] == u2) ++ idx;
						u = edges[idx];
						edges[idx] = u1;
						if(!first_time) modified_edges.pb(mp(mp(pre,u), u1));

						u2 = pre;
						while(u != u1) {
							is[u] = 0;
							ui tmp = u;
							if(edges[pstart[u]] == pre) u = edges[pstart[u]+1];
							else u = edges[pstart[u]];
							S.pb(mp(tmp, u));
							pre = tmp;
						}

						edge_rewire(u1, pend, pre, u2);
						if(!first_time) modified_edges.pb(mp(mp(u1,pre), u2));
					}
				}
				else {
					ui v2 = pre, v1 = pre;
					pre = u2;
					while(v1 != u1) {
						is[v1] = 0;
						int tmp = v1;
						if(edges[pstart[v1]] == pre) v1 = edges[pstart[v1]+1];
						else v1 = edges[pstart[v1]];
						S.pb(mp(tmp, v1));
						pre = tmp;
					}
					v1 = pre;
					if(exist_edge(u1, u2, pend)) {
						-- degree[u1];
						-- degree[u2];
#ifndef NDEBUG
						if(degree[u1] <= 1||degree[u2] <= 1) printf("WA xx\n");
#endif
						if(degree[u1] == 2) degree_twos.pb(u1);
						if(degree[u2] == 2) degree_twos.pb(u2);
					}
					else {
						edge_rewire(u1, pend, v1, u2);
						edge_rewire(u2, pend, v2, u1);
						if(!first_time) {
							modified_edges.pb(mp(mp(u1,v1), u2));
							modified_edges.pb(mp(mp(u2,v2), u1));
						}
					}
				}
			}
		}

		if(first_time) {
			S_size = (int)S.size();
			first_time = 0;
			for(ui k = 0;k < n;k ++) {
				if(is[k]&&degree[k] > 0) {
					++ kernal_size;
					for(ui j = pstart[k];j < pend[k];j ++) if(is[edges[k]]) ++ kernal_edges;
				}
				else fixed[k] = 1;
			}
		}

        // INEXACT REDUCTIONS

		// while(degree_ones.empty()&&degree_twos.empty()) {
		// 	while(max_d >= 3&&bin_head[max_d] == -1) -- max_d;
		// 	if(max_d < 3) break;

		// 	int v = -1;
		// 	for(v = bin_head[max_d];v != -1;) {
		// 		int tmp = bin_next[v];
		// 		if(is[v]&&degree[v] > 0) {
		// 			if(degree[v] < max_d) {
		// 				bin_next[v] = bin_head[degree[v]];
		// 				bin_head[degree[v]] = v;
		// 			}
		// 			else {
		// 				S.pb(mp(v,n)); ++ inexact;

		// 				res += delete_vertex(v, pend, is, degree, degree_ones, degree_twos);

		// 				bin_head[max_d] = tmp;
		// 				break;
		// 			}
		// 		}
		// 		v = tmp;
		// 	}
		// 	if(v == -1) bin_head[max_d] = -1;
		// }
	// }

    //UNDOING STUFF

	// for(int i = S.size()-1;i >= 0;i --) {
	// 	ui u1 = S[i].first, u2 = S[i].second;
	// 	assert(!is[u1]);

	// 	if(u2 != n) {
	// 		if(!is[u2]) {
	// 			is[u1] = 1;
	// 			++ res;
	// 		}
	// 		continue;
	// 	}

	// 	int ok = 1;
	// 	for(ui i = pstart[u1];i < pstart[u1+1];i ++) if(is[edges[i]]) {
	// 		ok = 0;
	// 		break;
	// 	}
	// 	if(ok) {
	// 		is[u1] = 1;
	// 		++ res;
	// 	}
	// }

#ifdef _LINUX_
	gettimeofday(&end1, NULL);

	long long mtime1, seconds1, useconds1;
	seconds1 = end1.tv_sec - start.tv_sec;
	useconds1 = end1.tv_usec - start.tv_usec;
	mtime1 = seconds1*1000000 + useconds1;
#endif

	//printf("Degree_two_path MIS: %d (kernal (|V|,|E|): (%d,%d), inexact reduction: %d)\n", res, kernal_size, kernal_edges, inexact);


    long offset = 0;
    for (ui i = 0; i < n; i++) {
        if (fixed[i]) {
            if (is[i]) offset++;
            continue;
        }
    }
    //if(S.size() % 2 != 0) {
        //printf("ERROR! S.size() = %d", S.size());
    //}
    offset += S.size() / 2;
    //printf("Offset = %d\n", offset);
    // write_graph_metis(is, degree, pend);
    write_kernel_to_adj_list(is, degree, pend, out_kernel);

	delete[] bin_head;
	delete[] bin_next;
	delete[] degree;
	delete[] pend;


	for(int i = (int)modified_edges.size()-1;i >= 0;i --) {
		ui u = modified_edges[i].first.first, u1 = modified_edges[i].first.second, u2 = modified_edges[i].second;
		for(ui j = pstart[u];j < pstart[u+1];j ++) if(edges[j] == u2) {
			edges[j] = u1;
			break;
		}
	}

	/*if(inexact) {
		res = ARW(is, fixed, mtime1, time_limit);
		for(int i = S_size-1;i >= 0;i --) {
			ui u1 = S[i].first, u2 = S[i].second;

			if(!is[u2]) is[u1] = 1;
			else is[u1] = 0;
		}
	}*/

#ifndef NDEBUG
	//compute_upperbound(is, fixed);
	swap(edges, tmp_edges);
	delete[] tmp_edges;
	swap(pstart, tmp_start);
	delete[] tmp_start;
	check_is(is, res);
#endif

	// delete[] is;
	// delete[] fixed;

#ifdef _LINUX_
	gettimeofday(&end, NULL);

	long long mtime, seconds, useconds;
	seconds = end.tv_sec - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = seconds*1000000 + useconds;

	//printf("Process time: %lld, Swap time: %lld, Total time: %lld\n", mtime1, mtime-mtime1, mtime);
#endif
    return offset;
}

void Graph::UndoReductions(std::vector<bool> &in_out_is) {

    for(unsigned int i = 0; i < in_out_is.size(); ++i) {
        is[reverse_mapping[i]] = in_out_is[i] ? 1 : 0;
    }
    
    for(int i = S.size()-1;i >= 0;i --) {
        ui u1 = S[i].first, u2 = S[i].second;

        if(!is[u2]) is[u1] = 1;
        else is[u1] = 0;
    }
    in_out_is.resize(n);
    for(unsigned int i = 0; i < n; ++i) {
        in_out_is[i] = is[i];
    }
	delete[] is;
	delete[] fixed;
}

int Graph::compute_triangle_counts(ui *tri, ui *pend, char *adj, char *is, int *degree, char *dominate, vector<ui> &dominated) {
	vector<ui> vs;
	for(ui i = 0;i < n;i ++) if(is[i]&&degree[i] > 0) vs.pb(i);

	ui *ids = new ui[vs.size()];
	memset(ids, 0, sizeof(ui)*vs.size());
	for(ui i = 0;i < vs.size();i ++) ++ ids[degree[vs[i]]];
	for(ui i = 1;i < vs.size();i ++) ids[i] += ids[i-1];

	int *order = new int[n];
	memset(order, -1, sizeof(int)*n);
	for(ui i = 0;i < vs.size();i ++) order[vs[i]] = (-- ids[degree[vs[i]]]);

	for(ui i = 0;i < vs.size();i ++) ids[order[vs[i]]] = vs[i];

	int res = 0;
	for(int i = (int)vs.size()-1;i >= 0;i --) {
		ui u = ids[i];
		shrink(u, pend[u], is, tri);

		for(ui j = pstart[u];j < pend[u];j ++) adj[edges[j]] = 1;

		if(!dominate[u]) {
			adj[u] = 1;
			for(ui j = pstart[u];j < pend[u];j ++) if(degree[edges[j]] <= degree[u]) {
				int tri_cnt = 0;
				ui v = edges[j];
				for(ui k = pstart[v];k < pend[v];k ++) if(is[edges[k]]) {
					if(!adj[edges[k]]) break;
					++ tri_cnt;
				}
				if(tri_cnt == degree[v]) {
					dominate[u] = 1;
					break;
				}
			}
			adj[u] = 0;
		}

		if(dominate[u]) {
			is[u] = 0;
			for(ui j = pstart[u];j < pend[u];j ++) if((-- degree[edges[j]]) == 0) ++ res;
			for(ui j = pstart[u];j < pend[u];j ++) if(order[edges[j]] > i) {
				ui v = edges[j];
				for(ui k = pstart[v];k < pend[v];k ++) if(is[edges[k]]) {
					ui w = edges[k];
					if(adj[w]) {
						assert(tri[k] > 0);
						-- tri[k];
					}
					if(!dominate[v]&&tri[k]+1 == degree[w]) {
						dominate[v] = 1;
						dominated.pb(v);
					}
					if(!dominate[w]&&tri[k]+1 == degree[v]) {
						dominate[w] = 1;
						if(order[w] > i) dominated.pb(w);
					}
#ifndef NDEBUG
					if(order[w] == i) printf("%d %d %d %d\n", w, order[w], u, i);
#endif
					assert(degree[v] > 1||dominate[w]);
					assert(degree[w] > 1||dominate[v]);
				}
			}
		}
		else {
#ifndef NDEBUG
			for(ui j = pstart[u];j < pend[u];j ++) if(degree[edges[j]] == 1) printf("WA dominate!\n");
#endif
			for(ui j = pstart[u];j < pend[u];j ++) {
				ui v = edges[j];
				tri[j] = 0;
				for(ui k = pstart[v];k < pend[v];k ++) if(adj[edges[k]]) ++ tri[j];
				assert(tri[j]+1 != degree[v]);
				if(!dominate[v]&&tri[j]+1 == degree[u]) {
					dominate[v] = 1;
					if(order[v] > i) dominated.pb(v);
				}
			}
		}

		for(ui j = pstart[u];j < pend[u];j ++) adj[edges[j]] = 0;
	}

	delete[] order;
	delete[] ids;

	return res;
}

int Graph::dominated_check(ui u, ui *pend, char *is, ui *tri, int *degree) {
	for(ui i = pstart[u];i < pend[u];i ++) if(is[edges[i]]&&tri[i]+1 == degree[edges[i]]) return 1;
	return 0;
}

int Graph::delete_vertex(ui u, ui *pend, char *is, vector<ui> &degree_twos, ui *tri, char *adj, int *degree, char *dominate, vector<ui> &dominated) {
	is[u] = 0;
	int res = 0;
	shrink(u, pend[u], is, tri);
	for(ui i = pstart[u];i < pend[u];i ++) {
		ui v = edges[i];
		adj[v] = 1;
		-- degree[v];

		if(degree[v] == 0) ++ res;
		else if(degree[v] == 2) degree_twos.pb(v);
	}

	for(ui i = pstart[u];i < pend[u];i ++) {
		ui v = edges[i];

		for(ui j = pstart[v];j < pend[v];j ++) if(is[edges[j]]) {
			ui w = edges[j];

#ifndef NDEBUG
			if(adj[w]&&tri[j] == 0) printf("WA delete_vertex!\n");
#endif
			if(adj[w]) -- tri[j];

			if(tri[j]+1 == degree[v]&&!dominate[w]) {
				dominate[w] = 1;
				dominated.pb(w);
			}
		}
	}

	for(ui i = pstart[u];i < pend[u];i ++) adj[edges[i]] = 0;

	return res;
}

void Graph::update_triangle(ui u1, ui u2, ui *pend, char *is, char *adj, ui *tri, int *degree, char *dominate, vector<ui> &dominated) {
	int cnt = 0;
	shrink(u1, pend[u1], is, tri);
	shrink(u2, pend[u2], is, tri);

	for(ui i = pstart[u1];i < pend[u1];i ++) adj[edges[i]] = 1;

	for(ui i = pstart[u2];i < pend[u2];i ++) if(adj[edges[i]]) {
		ui v = edges[i];
		tri[i] ++;
		++ cnt;
		if(tri[i]+1 == degree[u2]&&!dominate[v]) {
			dominate[v] = 1;
			dominated.pb(v);
		}
		if(tri[i]+1 == degree[v]&&!dominate[u2]) {
			dominate[u2] = 1;
			dominated.pb(u2);
		}

		for(ui j = pstart[v];j < pend[v];j ++) if(edges[j] == u2) {
			tri[j] ++;
			break;
		}
	}

	for(ui i = pstart[u1];i < pend[u1];i ++) {
		adj[edges[i]] = 0;
		if(edges[i] == u2) tri[i] = cnt;
	}

	if(cnt+1 == degree[u1]&&!dominate[u2]) {
		dominate[u2] = 1;
		dominated.pb(u2);
	}
	if(cnt+1 == degree[u2]&&!dominate[u1]) {
		dominate[u1] = 1;
		dominated.pb(u1);
	}

	for(ui i = pstart[u2];i < pend[u2];i ++) adj[edges[i]] = 1;

	for(ui i = pstart[u1];i < pend[u1];i ++) if(adj[edges[i]]) {
		ui v = edges[i];
		tri[i] ++;
		if(tri[i]+1 == degree[u1]&&!dominate[v]) {
			dominate[v] = 1;
			dominated.pb(v);
		}
		if(tri[i]+1 == degree[v]&&!dominate[u1]) {
			dominate[u1] = 1;
			dominated.pb(u1);
		}

		for(ui j = pstart[v];j < pend[v];j ++) if(edges[j] == u1) {
			tri[j] ++;
			break;
		}
	}

	for(ui i = pstart[u2];i < pend[u2];i ++) {
		adj[edges[i]] = 0;
		if(edges[i] == u1) tri[i] = cnt;
	}
}

void Graph::get_two_neighbors(ui u, char *is, ui &u1, ui &u2) {
	for(ui i = pstart[u];i < pstart[u+1];i ++) if(is[edges[i]]) {
		if(u1 == n) {
			u1 = edges[i];
			swap(edges[pstart[u]], edges[i]);
		}
		else {
			assert(u2 == n);
			u2 = edges[i];
			swap(edges[pstart[u]+1], edges[i]);
			break;
		}
	}
}

ui Graph::get_other_neighbor(ui u, char *is, ui u1) {
	ui idx = 0;
	for(ui i = pstart[u];i < pstart[u+1];i ++) if(is[edges[i]]) {
		swap(edges[pstart[u]+idx], edges[i]);
		if((++ idx) == 2) break;
	}
	assert(edges[pstart[u]] == u1||edges[pstart[u]+1] == u1);
	if(edges[pstart[u]] == u1) return edges[pstart[u]+1];
	return edges[pstart[u]];
}

int Graph::exist_edge(ui u1, ui u2) {
	if(pstart[u1+1]-pstart[u1] < pstart[u2+1]-pstart[u2]) {
		for(ui i = pstart[u1];i < pstart[u1+1];i ++) if(edges[i] == u2) return 1;
		return 0;
	}
	for(ui i = pstart[u2];i < pstart[u2+1];i ++) if(edges[i] == u1) return 1;
	return 0;
}

void Graph::edge_rewire(ui u, ui u1, ui u2) {
	for(ui i = pstart[u];i < pstart[u+1];i ++) if(edges[i] == u1) {
		edges[i] = u2;
		break;
	}
}

int Graph::remove_degree_one_two(vector<ui> &degree_ones, vector<ui> &degree_twos, char *is, int *degree, vector<pair<ui,ui> > &S) {
	int res = 0;

	while(!degree_ones.empty()||!degree_twos.empty()) {
		while(!degree_ones.empty()) {
			ui u = degree_ones.back();
			degree_ones.pop_back();
			if(!is[u]||degree[u] != 1) continue;

			for(int j = pstart[u];j < pstart[u+1];j ++) if(is[edges[j]]) {
				res += delete_vertex(edges[j], is, degree, degree_ones, degree_twos);
			}
		}

		while(!degree_twos.empty()&&degree_ones.empty()) {
			ui u = degree_twos.back();
			degree_twos.pop_back();
			if(!is[u]||degree[u] != 2) continue;

			ui u1 = n, u2 = n;
			get_two_neighbors(u, is, u1, u2);
			assert(u1 != n&&u2 != n);

			ui pre = u, cnt = 1;
			while(u1 != u&&degree[u1] == 2) {
				++ cnt;
				ui tmp = get_other_neighbor(u1, is, pre);
				pre = u1;
				u1 = tmp;
			}
			if(u1 == u) {
				res += delete_vertex(u, is, degree, degree_ones, degree_twos);
				assert(!degree_ones.empty());
				continue;
			}

			pre = u;
			while(degree[u2] == 2) {
				++ cnt;
				ui tmp = get_other_neighbor(u2, is, pre);
				pre = u2;
				u2 = tmp;
			}
			if(u1 == u2) {
				res += delete_vertex(u1, is, degree, degree_ones, degree_twos);
				assert(!degree_ones.empty());
				continue;
			}

			if(cnt%2 == 1) {
				if(exist_edge(u1, u2)) {
					res += delete_vertex(u1, is, degree, degree_ones, degree_twos);
					res += delete_vertex(u2, is, degree, degree_ones, degree_twos);
					assert(!degree_ones.empty());
				}
				else if(cnt > 1) {
					u = get_other_neighbor(pre, is, u2);
					edge_rewire(pre, u, u1);

					u2 = pre;
					while(u != u1) {
						is[u] = 0;
						ui tmp = get_other_neighbor(u, is, pre);
						S.pb(mp(u, tmp));
						pre = u;
						u = tmp;
					}

					edge_rewire(u1, pre, u2);
				}
			}
			else {
				ui v2 = pre, v1 = pre;
				pre = u2;
				while(v1 != u1) {
					is[v1] = 0;
					ui tmp = get_other_neighbor(v1, is, pre);
					S.pb(mp(v1,tmp));
					pre = v1;
					v1 = tmp;
				}
				v1 = pre;
				if(exist_edge(u1, u2)) {
					if((-- degree[u1]) == 2) degree_twos.pb(u1);
					if((-- degree[u2]) == 2) degree_twos.pb(u2);

					assert(degree[u1] > 1&&degree[u2] > 1);
				}
				else {
					edge_rewire(u1, v1, u2);
					edge_rewire(u2, v2, u1);
				}
			}
		}
	}

	return res;
}

int Graph::initial_dominance_and_degree_two_remove(vector<ui> &degree_ones, vector<ui> &degree_twos, char *is, int *degree, char *adj, vector<pair<ui,ui> > &S) {
	int res = 0;
	res += remove_degree_one_two(degree_ones, degree_twos, is, degree, S);

	// sort vertices in degree increasing order
	ui *ids = new ui[n];
	memset(ids, 0, sizeof(ui)*n);
	for(ui i = 0;i < n;i ++) ++ ids[degree[i]];
	for(ui i = 1;i < n;i ++) ids[i] += ids[i-1];

	ui *order = new ui[n];
	for(ui i = 0;i < n;i ++) {
		order[i] = ids[degree[i]];
		-- ids[degree[i]];
	}
	for(ui i = 0;i < n;i ++) ids[order[i]-1] = i;

	// compute dominance for vertices in degree decreasing order
	for(int i = n-1;i >= 0;i --) {
		int u = ids[i];
		if(!is[u]||degree[u] <= 0) continue;

		for(ui j = pstart[u];j < pstart[u+1];j ++) if(is[edges[j]]) adj[edges[j]] = 1;
		adj[u] = 1;

		int dominate = 0;
		for(ui j = pstart[u];j < pstart[u+1];j ++) if(is[edges[j]]&&degree[edges[j]] <= degree[u]) {
			int tri_cnt = 0;
			ui v = edges[j];
			for(ui k = pstart[v];k < pstart[v+1];k ++) if(is[edges[k]]) {
				if(!adj[edges[k]]) break;
				++ tri_cnt;
			}
			if(tri_cnt == degree[v]) {
				dominate = 1;
				break;
			}
		}

		for(ui j = pstart[u];j < pstart[u+1];j ++) adj[edges[j]] = 0;
		adj[u] = 0;

		if(dominate) {
			res += delete_vertex(u, is, degree, degree_ones, degree_twos);
			res += remove_degree_one_two(degree_ones, degree_twos, is, degree, S);
		}
	}

	delete[] ids;
	delete[] order;

	return res;
}

int Graph::lp_reduction(ui *ids, ui ids_n, char *is, int *degree) {
	// printf("begin lp_reduction\n");
	int *new_id = new int[n];
	memset(new_id, -1, sizeof(int)*n);
	for(ui i = 0;i < ids_n;i ++) new_id[ids[i]] = i;
	for(ui i = 0;i < ids_n;i ++) {
		ui u = ids[i];
		for(ui j = pstart[u];j < pstart[u+1];j ++) edges[j] = new_id[edges[j]];
	}
	delete[] new_id;
	new_id = NULL;

	ui *new_pstart = new ui[ids_n+1];
	for(ui i = 0;i < ids_n;i ++) new_pstart[i] = pstart[ids[i]];
	new_pstart[ids_n] = pstart[ids[ids_n-1]+1];

	int res = 0;
	int *x = new int[ids_n];
	int *y = new int[ids_n];
	memset(x, -1, sizeof(int)*ids_n);
	memset(y, -1, sizeof(int)*ids_n);

	char *used = new char[2*ids_n];
	int *level = new int[2*ids_n];
	ui *que = new ui[2*ids_n];
	ui *iter = new ui[ids_n];

	while(true) {
		memset(used, 0, sizeof(char)*2*ids_n);
		ui que_n = 0;
		for(ui u = 0;u < ids_n;u ++) if(x[u] == -1){
			level[u] = 0;
			used[u] = 1;
			que[que_n++] = u;
		}
		int find = 0;
		for(ui i = 0;i < que_n;i ++) {
			ui u = que[i];
			iter[u] = new_pstart[u+1];
			for(ui j = new_pstart[u];j < new_pstart[u+1];j ++) if(!used[ids_n+edges[j]]) {
				used[ids_n+edges[j]] = 1;
				int v = y[edges[j]];
				if(v == -1) find = 1;
				else {
					assert(!used[v]);
					used[v] = 1;
					level[v] = level[u]+1;
					que[que_n++] = v;
				}
			}
		}
		if(!find) break;

		for(ui i = 0;i < ids_n;i ++) if(x[i] == -1) {
			que_n = 0;
			que[que_n++] = i;
			while(que_n > 0) {
				ui u = que[que_n-1];
				if(iter[u] <= new_pstart[u]) {
					-- que_n;
					if(que_n > 0) -- iter[que[que_n-1]];
				}
				else {
					ui v = edges[iter[u]-1];
					if(y[v] == -1) {
						for(ui j = 0;j < que_n;j ++) {
							u = que[j];
							v = edges[iter[u]-1];
							x[u] = v; y[v] = u;
						}
						que_n = 0;
					}
					else if(level[y[v]] > level[u]) que[que_n++] = y[v];
					else -- iter[u];
				}
			}
		}
	}

	for(ui u = 0;u < ids_n;u ++) if(!used[u]&&used[ids_n+u]) {
		is[ids[u]] = 0;
		for(ui j = new_pstart[u];j < new_pstart[u+1];j ++) if(is[ids[edges[j]]]) {
			if( (-- degree[ids[edges[j]]]) == 0) ++ res;
		}
	}

#ifndef NDEBUG
	for(ui i = 0;i < ids_n;i ++) {
		ui u = ids[i];
		if(used[i]&&!used[ids_n+i]&&(!is[u]||degree[u] > 0)) printf("WA\n");
	}

	for(ui i = 0;i < ids_n;i ++) if((used[i]&&used[i+ids_n])||(!used[i]&&!used[i+ids_n])) {
		if(x[i] == -1||y[i] == -1) printf("WA half-integral solution! x[i]: %d, y[i]: %d\n", x[i], y[i]);
	}
#endif

	memset(used, 0, sizeof(char)*2*ids_n);
	for(ui i = 0;i < ids_n;i ++) iter[i] = new_pstart[i+1];

	ui level_n = 0;
	for(ui u = 0;u < ids_n;u ++) if(!used[u]&&is[ids[u]]&&degree[ids[u]] > 0) {
		ui que_n = 0;
		que[que_n ++] = u;
		used[u] = 1;
		while(que_n > 0) {
			u = que[que_n-1];
			if(u < ids_n) {
				assert(x[u] != -1&&is[ids[x[u]]]&&degree[ids[x[u]]] > 0);
				//assert(x[u] != -1); assert(is[ids[x[u]]]); assert(degree[ids[x[u]]] > 0);
				if(!used[x[u]+ids_n]) {
					used[x[u]+ids_n] = 1;
					que[que_n ++] = x[u]+ids_n;
				}
				else {
					level[level_n++] = u;
					-- que_n;
				}
			}
			else {
				u -= ids_n;
				ui v = -1;
				while(iter[u] > new_pstart[u]) {
					v = edges[iter[u]-1];
					if(!used[v]&&is[ids[v]]) {
						assert(degree[ids[v]] > 0);
						break;
					}
					-- iter[u];
				}
				if(iter[u] <= new_pstart[u]) -- que_n;
				else {
					used[v] = 1;
					que[que_n++] = v;
				}
			}
		}
	}

	memset(used, 0, sizeof(char)*2*ids_n);
	char *in = new char[ids_n];
	memset(in, 0, sizeof(char)*ids_n);

	for(int i = level_n - 1;i >= 0;i --) {
		ui u = level[i];
		if(used[u]||!is[ids[u]]||degree[ids[u]] == 0) continue;

		ui que_n = 0;
		que[que_n ++] = u;
		used[u] = 1;
		int ok = 1;
		for(ui j = 0;j < que_n;j ++) {
			u = que[j];
			if(ok&&used[u < ids_n? (u+ids_n) : (u-ids_n)]) ok = 0;
			if(u < ids_n) {
				for(ui k = new_pstart[u];k < new_pstart[u+1];k ++) {
					ui v = edges[k];
					if(!is[ids[v]]||degree[ids[v]] <= 0) continue;

					v += ids_n;
					if(!used[v]) {
						in[v-ids_n] = 1;
						used[v] = 1;
						que[que_n++] = v;
					}
					else if(!in[v-ids_n]) ok = 0;
				}
			}
			else {
				u -= ids_n;
				assert(y[u] != -1&&is[ids[y[u]]]&&degree[ids[y[u]]] > 0);
				if(!used[y[u]]) {
					used[y[u]] = 1;
					que[que_n++] = y[u];
				}
			}
		}

		for(ui j = 0;j < que_n;j ++) if(que[j] >= ids_n) in[que[j]-ids_n] = 0;

		if(ok) {
			for(ui j = 0;j < que_n;j ++) {
				u = que[j];
				if(u >= ids_n) {
					u -= ids_n;
					is[ids[u]] = 0;
					for(ui k = new_pstart[u];k < new_pstart[u+1];k ++) if(is[ids[edges[k]]]) {
						if( (-- degree[ids[edges[k]]]) == 0) ++ res;
					}
				}
			}
		}
	}

	delete[] new_pstart;
	delete[] iter;
	delete[] level;
	delete[] que;
	delete[] used;
	delete[] x;
	delete[] y;
	delete[] in;

	for(ui i = 0;i < ids_n;i ++) {
		ui u = ids[i];
		for(ui j = pstart[u];j < pstart[u+1];j ++) edges[j] = ids[edges[j]];
	}
	//printf("finish lp reduction!\n");

	return res;
}

void Graph::degree_two_kernal_dominate_lp_and_remove_max_degree_without_contraction() {
#ifndef NDEBUG
	ui *tmp_edges = new ui[m];
	for(ui i = 0;i < m;i ++) tmp_edges[i] = edges[i];
	ui *tmp_pstart = new ui[n+1];
	for(ui i = 0;i <= n;i ++) tmp_pstart[i] = pstart[i];
#endif

#ifdef _LINUX_
	struct timeval start, end1, end;
	gettimeofday(&start, NULL);
#endif

	char *is = new char[n];
	for(ui i = 0;i < n;i ++) is[i] = 1;
	char *adj = new char[n];
	memset(adj, 0, sizeof(char)*n);

	int res = 0;
	vector<ui> degree_ones, degree_twos;
	int *degree = new int[n];
	for(ui i = 0;i < n;i ++) {
		degree[i] = pstart[i+1]-pstart[i];
		if(degree[i] == 0) ++ res;
		else if(degree[i] == 1) degree_ones.pb(i);
		else if(degree[i] == 2) degree_twos.pb(i);
	}

	vector<pair<ui,ui> > S;
	vector<pair<pair<ui,ui>, ui> > modified_edges;

	res += initial_dominance_and_degree_two_remove(degree_ones, degree_twos, is, degree, adj, S);

#ifdef _LINUX_
	struct timeval end_t1;
	gettimeofday(&end_t1, NULL);
#endif

#ifndef NDEBUG
	int new_n = 0;
	for(ui i = 0;i < n;i ++) {
		//if(is[i]&&degree[i] == 1) printf("WA degree one exist!\n");
		if(is[i]&&degree[i] > 0) ++ new_n;
	}
	printf("initial dominance and degree two remove: %d -> %d\n", n, new_n);
#endif

	ui *ids = new ui[n];
	ui ids_n = 0;

	ui new_m = 0;
	for(ui i = 0;i < n;i ++) {
		ui tmp = pstart[i];
		pstart[i] = new_m;
		if(!is[i]||degree[i] <= 0) continue;

		ids[ids_n ++] = i;
		for(ui j = tmp;j < pstart[i+1];j ++) if(is[edges[j]]) {
			assert(degree[edges[j]] > 0);
			edges[new_m ++] = edges[j];
		}
	}
	pstart[n] = new_m;

	// if(ids_n > 0) res += lp_reduction(ids, ids_n, is, degree);

#ifdef _LINUX_
	struct timeval end_t2;
	gettimeofday(&end_t2, NULL);
#endif

#ifndef NDEBUG
	int new_n2 = 0;
	for(ui i = 0;i < n;i ++) {
		if(is[i]&&degree[i] > 0) ++ new_n2;
	}
	printf("lp reduction: %d -> %d\n", new_n, new_n2);
#endif

	assert(degree_ones.empty()&&degree_twos.empty());
	for(ui i = 0;i < ids_n;i ++) if(is[ids[i]]&&degree[ids[i]] > 0) {
		if(degree[ids[i]] == 1) degree_ones.pb(ids[i]);
		else if(degree[ids[i]] == 2) degree_twos.pb(ids[i]);
	}

	res += remove_degree_one_two(degree_ones, degree_twos, is, degree, S);

	delete[] ids; ids = NULL;

	new_m = 0;
	for(ui i = 0;i < n;i ++) {
		ui tmp = pstart[i];
		pstart[i] = new_m;
		if(!is[i]||degree[i] <= 0) continue;

		for(ui j = tmp;j < pstart[i+1];j ++) if(is[edges[j]]) {
			assert(degree[edges[j]] > 0);
			edges[new_m ++] = edges[j];
		}
	}
	pstart[n] = new_m;

	ui *pend = new ui[n];
	int max_dd = 0;
	for(ui i = 0;i < n;i ++) {
		pend[i] = pstart[i+1];
		if(pend[i] - pstart[i] > max_dd) max_dd = pend[i] - pstart[i];
	}
	//printf("max_d: %d, edges: %u\n", max_dd, new_m/2);

	int delete_tri = 0;
	ui *tri = NULL;
	if(new_m <= m/2) tri = edges+new_m;
	else {
		tri = new ui[pstart[n]];
		delete_tri = 1;
	}
	char *dominate = new char[n];
	memset(dominate, 0, sizeof(char)*n);
	vector<ui> dominated;

	res += compute_triangle_counts(tri, pend, adj, is, degree, dominate, dominated);

#ifdef _LINUX_
	struct timeval end_t3;
	gettimeofday(&end_t3, NULL);

	long long mtime11, seconds11, useconds11;
	seconds11 = end_t1.tv_sec - start.tv_sec;
	useconds11 = end_t1.tv_usec - start.tv_usec;
	mtime11 = seconds11*1000000 + useconds11;
	//printf("time1: %lld,", mtime11);
	seconds11 = end_t2.tv_sec - end_t1.tv_sec;
	useconds11 = end_t2.tv_usec - end_t1.tv_usec;
	mtime11 = seconds11*1000000 + useconds11;
	//printf(" time2: %lld,", mtime11);
	seconds11 = end_t3.tv_sec - end_t2.tv_sec;
	useconds11 = end_t3.tv_usec - end_t2.tv_usec;
	mtime11 = seconds11*1000000 + useconds11;
	//printf(" time3: %lld\n", mtime11);
#endif

#ifndef NDEBUG
	for(ui i = 0;i < n;i ++) if(is[i]&&!dominate[i]) {
		for(ui j = pstart[i];j < pstart[i+1];j ++) if(is[edges[j]]&&degree[edges[j]] == 1) {
			printf("WA degree one not dominate!\n");
		}
	}
#endif

	int *bin_head = new int[ids_n];
	int *bin_next = new int[n];
	memset(bin_head, -1, sizeof(int)*ids_n);

	assert(degree_twos.empty());
	int max_d = 0;
	for(ui i = 0;i < n;i ++) if(is[i]&&degree[i] > 0) {
		bin_next[i] = bin_head[degree[i]];
		bin_head[degree[i]] = i;

		if(degree[i] == 2) degree_twos.pb(i);
		if(degree[i] > max_d) max_d = degree[i];
	}

	char *fixed = new char[n];
	memset(fixed, 0, sizeof(char)*n);

	int kernal_size = 0, inexact = 0, first_time = 1, S_size = (int)S.size();
	int kernal_edges = 0;
	while(!dominated.empty()||!degree_twos.empty()) {
		while(!dominated.empty()||!degree_twos.empty()) {
			while(!dominated.empty()) {
				ui u = dominated.back();
				dominated.pop_back();
				if(!is[u]||degree[u] == 0) continue;

				if(!dominated_check(u, pend, is, tri, degree)) dominate[u] = 0;
				else res += delete_vertex(u, pend, is, degree_twos, tri, adj, degree, dominate, dominated);
			}

			while(!degree_twos.empty()&&dominated.empty()) {
				ui u = degree_twos.back();
				degree_twos.pop_back();
				if(!is[u]||degree[u] != 2) continue;

				shrink(u, pend[u], is, tri);
				assert(pend[u] == pstart[u] + 2);
				ui u1 = edges[pstart[u]], u2 = edges[pstart[u]+1];

				ui pre = u, cnt = 1;
				while(u1 != u&&degree[u1] == 2) {
					++ cnt;
					shrink(u1, pend[u1], is, tri);
					assert(pend[u1] == pstart[u1] + 2);
					int tmp = u1;
					if(edges[pstart[u1]] != pre) u1 = edges[pstart[u1]];
					else u1 = edges[pstart[u1]+1];
					pre = tmp;
				}

#ifndef NDEBUG
				if(u1 != u&&degree[u1] <= 2) {
					printf("%d u1(%d) degree_u1(%d) WAxx!\n", pre, u1, degree[u1]);
					printf("%d:", u1);
					for(ui k = pstart[u1];k < pend[u1];k ++) if(edges[k] != n&&is[edges[k]]) {
						printf(" %d(tri:%d,dominate:%d)", edges[k], tri[k], dominate[edges[k]]);
					}
					printf("\n");
				}
#endif
				if(u1 == u) {
					res += delete_vertex(u, pend, is, degree_twos, tri, adj, degree, dominate, dominated);
					assert(!dominated.empty());
					continue;
				}

				pre = u;
				while(degree[u2] == 2) {
					++ cnt;
					shrink(u2, pend[u2], is, tri);
					assert(pend[u2] == pstart[u2] + 2);
					int tmp = u2;
					if(edges[pstart[u2]] != pre) u2 = edges[pstart[u2]];
					else u2 = edges[pstart[u2]+1];
					pre = tmp;
				}
				if(u1 == u2) {
					res += delete_vertex(u1, pend, is, degree_twos, tri, adj, degree, dominate, dominated);
					assert(!dominated.empty());
					continue;
				}

				if(cnt%2 == 1) {
					if(exist_edge(u1,u2, pend)) {
						res += delete_vertex(u1, pend, is, degree_twos, tri, adj, degree, dominate, dominated);
						res += delete_vertex(u2, pend, is, degree_twos, tri, adj, degree, dominate, dominated);
						assert(!dominated.empty());
					}
					else if(cnt > 1) {
						ui idx = pstart[pre];
						if(edges[idx] == u2) ++ idx;
						assert(degree[pre] == 2&&tri[idx] == 0);
						u = edges[idx];
						edges[idx] = u1;
						if(!first_time) modified_edges.pb(mp(mp(pre,u), u1));

						u2 = pre;
						while(u != u1) {
							is[u] = 0;
							int tmp = u;
							if(edges[pstart[u]] == pre) u = edges[pstart[u]+1];
							else u = edges[pstart[u]];
							pre = tmp;
							S.pb(mp(pre,u));
						}
						idx = edge_rewire(u1, pend, pre, u2);
						assert(tri[idx] == 0);
						if(!first_time) modified_edges.pb(mp(mp(u1,pre),u2));
					}
				}
				else {
					ui v2 = pre, v1 = pre;
					pre = u2;
					while(v1 != u1) {
						is[v1] = 0;
						int tmp = v1;
						if(edges[pstart[v1]] == pre) v1 = edges[pstart[v1]+1];
						else v1 = edges[pstart[v1]];
						pre = tmp;
						S.pb(mp(pre,v1));
					}
					v1 = pre;
					if(exist_edge(u1, u2, pend)) {
						if((-- degree[u1]) == 2) degree_twos.pb(u1);
						if((-- degree[u2]) == 2) degree_twos.pb(u2);

						assert(degree[u1] > 1&&degree[u2] > 1);

						for(ui k = pstart[u1];k < pend[u1];k ++) {
							if(is[edges[k]]&&!dominate[edges[k]]&&tri[k]+1 == degree[u1]) {
								dominate[edges[k]] = 1;
								dominated.pb(edges[k]);
							}
						}
						for(ui k = pstart[u2];k < pend[u2];k ++) {
							if(is[edges[k]]&&!dominate[edges[k]]&&tri[k]+1 == degree[u2]) {
								dominate[edges[k]] = 1;
								dominated.pb(edges[k]);
							}
						}
					}
					else {
						ui idx = edge_rewire(u1, pend, v1, u2);
						assert(tri[idx] == 0);
						idx = edge_rewire(u2, pend, v2, u1);
						assert(tri[idx] == 0);

						if(!first_time) {
							modified_edges.pb(mp(mp(u1,v1),u2));
							modified_edges.pb(mp(mp(u2,v2),u1));
						}

						update_triangle(u1, u2, pend, is, adj, tri, degree, dominate, dominated);
					}
				}
			}
		}

		if(first_time) {
			first_time = 0;
			S_size = (int)S.size();
			for(ui k = 0;k < n;k ++) {
				if(is[k]&&degree[k] > 0) {
					++ kernal_size;
					for(ui j = pstart[k];j < pend[k];j ++) if(is[edges[j]]) ++ kernal_edges;
				}
				else fixed[k] = 1;
			}
		}

        // INEXACT REDUCTIONS

		// while(dominated.empty()&&degree_twos.empty()) {
		// 	while(max_d >= 3&&bin_head[max_d] == -1) -- max_d;
		// 	if(max_d < 3) break;

		// 	int v = -1;
		// 	for(v = bin_head[max_d];v != -1;) {
		// 		int tmp = bin_next[v];
		// 		if(is[v]&&degree[v] > 0) {
		// 			if(degree[v] < max_d) {
		// 				bin_next[v] = bin_head[degree[v]];
		// 				bin_head[degree[v]] = v;
		// 			}
		// 			else {
		// 				S.pb(mp(v,n)); ++ inexact;
		// 				res += delete_vertex(v, pend, is, degree_twos, tri, adj, degree, dominate, dominated);

		// 				bin_head[max_d] = tmp;
		// 				break;
		// 			}
		// 		}
		// 		v = tmp;
		// 	}
		// 	if(v == -1) bin_head[max_d] = -1;
		// }
	}

    // UNDO DEG 2 PATH AND INEXACT REDUCTIONS

	ui UB = 0;

// 	for(int i = S.size()-1;i >= 0;i --) {
// 		ui u1 = S[i].first, u2 = S[i].second;
// 		assert(is[u1] == 0);

// 		if(u2 != n) {
// 			if(!is[u2]) {
// 				is[u1] = 1;
// #ifndef NDEBUG
// 				for(ui j = tmp_pstart[u1];j < tmp_pstart[u1+1];j ++) if(is[tmp_edges[j]]) {
// 					printf("WA conflict1!\n");
// 				}
// #endif
// 				++ res;
// 			}
// 			continue;
// 		}

// 		int ok = 1;
// 		for(ui j = pstart[u1];j < pstart[u1+1];j ++) if(is[edges[j]]) {
// 			ok = 0;
// 			break;
// 		}
// 		if(ok) {
// 			is[u1] = 1;
// 			++ res;

// #ifndef NDEBUG
// 			for(ui j = tmp_pstart[u1];j < tmp_pstart[u1+1];j ++) if(is[tmp_edges[j]]) {
// 				printf("WA conflict2!\n");
// 			}
// #endif
// 		}
// 		else ++ UB;
// 	}

	printf("Degree_two_dominate_lp MIS: %d (kernal (|V|,|E|): (%d,%d), inexact reduction: %d, UB: %d)\n", res, kernal_size, kernal_edges, inexact, res+UB);

	delete[] bin_head;
	delete[] bin_next;
	delete[] degree;

	delete[] pend;
	if(delete_tri) delete[] tri;
	delete[] dominate;
	delete[] adj;

#ifdef _LINUX_
	gettimeofday(&end1, NULL);

	long long mtime1, seconds1, useconds1;
	seconds1 = end1.tv_sec - start.tv_sec;
	useconds1 = end1.tv_usec - start.tv_usec;
	mtime1 = seconds1*1000000 + useconds1;
#endif


	for(int i = (int)modified_edges.size()-1;i >= 0;i --) {
		ui u = modified_edges[i].first.first, u1 = modified_edges[i].first.second, u2 = modified_edges[i].second;
		for(ui j = pstart[u];j < pstart[u+1];j ++) if(edges[j] == u2) {
			edges[j] = u1;
			break;
		}
	}

	/*
	if(inexact) {
		res = ARW(is, fixed, mtime1, time_limit);
		for(int i = S_size-1;i >= 0;i --) {
			ui u1 = S[i].first, u2 = S[i].second;

			if(!is[u2]) is[u1] = 1;
			else is[u1] = 0;
		}
	}
	*/

#ifndef NDEBUG
	compute_upperbound(is, fixed);
	swap(edges, tmp_edges);
	delete[] tmp_edges;
	swap(pstart, tmp_pstart);
	delete[] tmp_pstart;

	check_is(is, res);
#endif

	delete[] is;
	delete[] fixed;

#ifdef _LINUX_
	gettimeofday(&end, NULL);

	long long mtime, seconds, useconds;
	seconds = end.tv_sec - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = seconds*1000000 + useconds;

	printf("Process time: %lld, Swap time: %lld, Total time: %lld\n", mtime1, mtime-mtime1, mtime);
#endif
}

int Graph::delete_vertex(ui v, char *is, int *degree, int *head, Edge *es, int *bin_head, int *bin_next, int *bin_pre, vector<ui> &degree_ones, vector<ui> &degree_twos) {
	is[v] = 0;

	int res = 0;
	for(int k = head[v];k != -1;k = es[k].next) if(is[es[k].id]) {
		int w = es[k].id;

		if(bin_pre[w] != -1) bin_next[bin_pre[w]] = bin_next[w];
		else {
#ifndef NDEBUG
			if(bin_head[degree[w]] != w) printf("WWAA1 bin_head[degree[w]]: %d, w = %d, degree[w]: %d\n", bin_head[degree[w]], w, degree[w]);
#endif
			bin_head[degree[w]] = bin_next[w];
		}
		if(bin_next[w] != -1) bin_pre[bin_next[w]] = bin_pre[w];

		-- degree[w];

		bin_next[w] = bin_head[degree[w]];
		bin_pre[w] = -1;
		if(bin_head[degree[w]] != -1) bin_pre[bin_head[degree[w]]] = w;
		bin_head[degree[w]] = w;

		if(degree[w] == 0) ++ res;
		else if(degree[w] == 1) degree_ones.pb(w);
		else if(degree[w] == 2) degree_twos.pb(w);

#ifndef NDEBUG
		int cntx = 0;
		for(int ii = head[w];ii != -1;ii = es[ii].next) if(is[es[ii].id]) ++ cntx;
		if(cntx != degree[w]) {
			//for(int ii = pstart[w];ii < pstart[w+1];ii ++) printf("%d %d\n", edges[ii], vis[edges[ii]]);

			printf("WA degree1 %d: %d %d!\n", w, cntx, degree[w]);
		}
#endif
	}
	return res;
}

void Graph::degree_two_kernal_and_remove_max_degree_with_contraction() {
	int *head = new int[n];
	Edge *es = new Edge[m];
	int cnt = 0;
	memset(head, -1, sizeof(int)*n);
	for(int i = 0;i < n;i ++) {
		for(int j = pstart[i];j < pstart[i+1];j ++) if(edges[j] > i) {
			int a = i, b = edges[j];

			es[cnt].id = b;
			es[cnt].duplicate = cnt+1;
			es[cnt].next = head[a];
			head[a] = cnt ++;

			es[cnt].id = a;
			es[cnt].duplicate = cnt-1;
			es[cnt].next = head[b];
			head[b] = cnt ++;
		}
	}

#ifdef _LINUX_
	struct timeval start, end1, end;
	gettimeofday(&start, NULL);
#endif

	char *is = new char[n];
	for(ui i = 0;i < n;i ++) is[i] = 1;

	int *bin_head = new int[n];
	int *bin_next = new int[n];
	int *bin_pre = new int[n];
	int *degree = new int[n];
	memset(bin_head, -1, sizeof(int)*n);

	int *includes = new int[n];
	int *excludes = new int[n];
	memset(includes, -1, sizeof(int)*n);
	memset(excludes, -1, sizeof(int)*n);

	vector<ui> degree_ones, degree_twos;

	int max_d = 0, res = 0;
	for(ui i = 0;i < n;i ++) {
		degree[i] = pstart[i+1] - pstart[i];
		bin_next[i] = bin_head[degree[i]];
		bin_pre[i] = -1;
		if(bin_head[degree[i]] != -1) bin_pre[bin_head[degree[i]]] = i;
		bin_head[degree[i]] = i;

		if(degree[i] == 0) ++ res;
		else if(degree[i] == 1) degree_ones.pb(i);
		else if(degree[i] == 2) degree_twos.pb(i);
		if(degree[i] > max_d) max_d = degree[i];
	}

	char *in = new char[n];
	memset(in, 0, sizeof(char)*n);
	char *fixed = new char[n];
	memset(fixed, 0, sizeof(char)*n);
	char *contracted = new char[n];
	memset(contracted, 0, sizeof(char)*n);

	int kernal_size = 0, inexact = 0, first_time = 1;
	while(!degree_ones.empty()||!degree_twos.empty()||max_d >= 3) {
		while(!degree_ones.empty()||!degree_twos.empty()) {
			while(!degree_ones.empty()) {
				int u = degree_ones.back();
				degree_ones.pop_back();
				if(!is[u]||degree[u] != 1) continue;

#ifndef NDEBUG
				if(degree[u] > 1) printf("WW degree one1!\n");
#endif
				cnt = 0;
				for(int j = head[u];j != -1;j = es[j].next) if(is[es[j].id]) {
					int v = es[j].id;
#ifndef NDEBUG
					if((++ cnt) > 1) printf("WA degree one!\n");
#endif
					res += delete_vertex(v, is, degree, head, es, bin_head, bin_next, bin_pre, degree_ones, degree_twos);
				}
			}
			if(first_time) {
				first_time = 0;
				for(ui j = 0;j < n;j ++) {
					if(is[j]&&degree[j] > 0) ++ kernal_size;
					else fixed[j] = 1;
				}
			}
			while(degree_ones.empty()&&!degree_twos.empty()) {
				int u = degree_twos.back();
				degree_twos.pop_back();
				if(!is[u]||degree[u] != 2) continue;

				int v1 = -1, v2 = -1;
				for(int j = head[u];j != -1;j = es[j].next) if(is[es[j].id]) {
					if(v1 == -1) v1 = es[j].id;
					else if(v2 == -1) v2 = es[j].id;
					else printf("WA degree two!\n");
				}

#ifndef NDEBUG
				if(v2 == -1) printf("WWWAAA\n");
#endif

				char find = 0;
				for(int j = head[v1];j != -1;j = es[j].next) if(es[j].id == v2) find = 1;

				if(find) {
					res += delete_vertex(v1, is, degree, head, es, bin_head, bin_next, bin_pre, degree_ones, degree_twos);
					res += delete_vertex(v2, is, degree, head, es, bin_head, bin_next, bin_pre, degree_ones, degree_twos);
				}
				else {
					contracted[u] = contracted[v2] = 1;
					is[u] = is[v2] = 0;
					++ res;
					int tmp = v1;
					while(includes[tmp] != -1) tmp = includes[tmp];
					includes[tmp] = v2;
					if(excludes[v1] == -1) excludes[v1] = u;
					else {
						tmp = excludes[v1];
						while(includes[tmp] != -1) tmp = includes[tmp];
						includes[tmp] = u;
					}

					if(bin_pre[v1] != -1) bin_next[bin_pre[v1]] = bin_next[v1];
					else {
#ifndef NDEBUG
						if(bin_head[degree[v1]] != v1) printf("WWAAx\n");
#endif
						bin_head[degree[v1]] = bin_next[v1];
					}
					if(bin_next[v1] != -1) bin_pre[bin_next[v1]] = bin_pre[v1];

					int pre_degree = degree[v1];
					-- degree[v1];

					int last = -1;
					for(int j = head[v1];j != -1;j = es[j].next) {
						in[es[j].id] = 1;
						last = j;
					}
#ifndef NDEBUG
					if(last == -1) printf("WA last\n");
#endif
					for(int j = head[v2];j != -1;j = es[j].next) if(is[es[j].id]) {
						if(!in[es[j].id]) {
							es[last].next = j;
							last = j;
							++ degree[v1];
							es[es[j].duplicate].id = v1;
						}
						else {
							int w = es[j].id;

							if(bin_pre[w] != -1) bin_next[bin_pre[w]] = bin_next[w];
							else {
#ifndef NDEBUG
								if(bin_head[degree[w]] != w) {
									printf("WWAAxx\n");
									printf("%d %d %d\n", w, degree[w], bin_head[degree[w]]);
								}
#endif
								bin_head[degree[w]] = bin_next[w];
							}
							if(bin_next[w] != -1) bin_pre[bin_next[w]] = bin_pre[w];

							-- degree[w];

							bin_next[w] = bin_head[degree[w]];
							bin_pre[w] = -1;
							if(bin_head[degree[w]] != -1) bin_pre[bin_head[degree[w]]] = w;
							bin_head[degree[w]] = w;
							if(degree[w] == 1) degree_ones.pb(w);
							else if(degree[w] == 2) degree_twos.pb(w);
						}
					}
					es[last].next = -1;

					for(int j = head[v1];j != -1;j = es[j].next) in[es[j].id] = 0;

					if(degree[v1] > max_d) max_d = degree[v1];
					bin_next[v1] = bin_head[degree[v1]];
					bin_pre[v1] = -1;
					if(bin_head[degree[v1]] != -1) bin_pre[bin_head[degree[v1]]] = v1;
					bin_head[degree[v1]] = v1;

					if(degree[v1] == 1&&pre_degree != 1) degree_ones.pb(v1);
					else if(degree[v1] == 2&&pre_degree != 2) degree_twos.pb(v1);
#ifndef NDEBUG
					int cntx = 0;
					for(int ii = head[v1];ii != -1;ii = es[ii].next) if(is[es[ii].id]) ++ cntx;
					if(cntx != degree[v1]) {
						//printf("**Adj %d:", v1);
						//for(int ii = head[v1];ii != -1;ii = es[ii].next) printf("\t%d,%d",es[ii].id, vis[es[ii].id]);
						//printf("\n");

						printf("WA degree4 %d %d %d!\n", cntx, degree[v1], pstart[v1+1]-pstart[v1]);
					}
#endif
				}
			}
		}

		while(degree_ones.empty()&&degree_twos.empty()&&max_d >= 3) {
			if(bin_head[max_d] == -1) -- max_d;
			else {
				int u = bin_head[max_d];
				bin_head[max_d] = bin_next[u];
				if(bin_next[u] != -1) bin_pre[bin_next[u]] = -1;
#ifndef NDEBUG
				if(bin_next[u] == u) printf("WAWA\n");
#endif
				if(!is[u]) continue;

				++ inexact;
				delete_vertex(u, is, degree, head, es, bin_head, bin_next, bin_pre, degree_ones, degree_twos);
			}
		}
	}

	for(ui i = 0;i < n;i ++) if(!contracted[i]) {
		vector<ui> backtrack;
		backtrack.pb(i);
		while(!backtrack.empty()) {
			ui u = backtrack.back();
			backtrack.pop_back();

			if(is[u]) {
				if(includes[u] != -1) {
#ifndef NDEBUG
					if(is[includes[u]]) printf("WA ss1\n");
#endif
					is[includes[u]] = 1;
					backtrack.pb(includes[u]);
				}
				if(excludes[u] != -1) {
#ifndef NDEBUG
					if(is[excludes[u]]) printf("WA ss2\n");
#endif
					is[excludes[u]] = 0;
					backtrack.pb(excludes[u]);
				}
			}
			else {
				if(includes[u] != -1) {
#ifndef NDEBUG
					if(is[includes[u]]) printf("WA tt1\n");
#endif
					is[includes[u]] = 0;
					backtrack.pb(includes[u]);
				}
				if(excludes[u] != -1) {
#ifndef NDEBUG
					if(is[excludes[u]]) printf("WA tt2\n");
#endif
					is[excludes[u]] = 1;
					backtrack.pb(excludes[u]);
				}
			}
		}
	}

	for(ui i = 0;i < n;i ++) if(!is[i]) {
		int find = 0;
		for(ui j = pstart[i];j < pstart[i+1];j ++) if(is[edges[j]]) {
			find = 1;
			break;
		}
		if(!find) {
			++ res;
			is[i] = 1;
		}
	}

	printf("Degree_two_c MIS: %d (kernal |V|: %d, inexact reduction: %d)\n", res, kernal_size, inexact);

	delete[] includes;
	delete[] excludes;
	delete[] head;
	delete[] es;
	delete[] bin_head;
	delete[] bin_pre;
	delete[] bin_next;
	delete[] degree;
	delete[] in;

#ifdef _LINUX_
	gettimeofday(&end1, NULL);

	long long mtime1, seconds1, useconds1;
	seconds1 = end1.tv_sec - start.tv_sec;
	useconds1 = end1.tv_usec - start.tv_usec;
	mtime1 = seconds1*1000000 + useconds1;
	//printf("Processing Time: %lld\n", mtime1);
#endif

#ifndef NDEBUG
	check_is(is, res);
	//compute_upperbound(is, fixed);
#endif

	//if(inexact) general_swap(is, fixed);

	delete[] is;
	delete[] fixed;

#ifdef _LINUX_
	gettimeofday(&end, NULL);

	long long mtime, seconds, useconds;
	seconds = end.tv_sec - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = seconds*1000000 + useconds;

	printf("Process time: %lld, Swap time: %lld, Total time: %lld\n", mtime1, mtime-mtime1, mtime);
#endif
}
}
