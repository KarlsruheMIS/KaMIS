# KaMIS v3.0 — Karlsruhe Maximum Independent Sets

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![C++11/14](https://img.shields.io/badge/C++-11%2F14-blue.svg?style=flat)
[![CMake](https://img.shields.io/badge/build-CMake-blue)](https://cmake.org/)
[![Linux](https://img.shields.io/badge/platform-Linux-blue)](https://github.com/KarlsruheMIS/KaMIS)
[![macOS](https://img.shields.io/badge/platform-macOS-blue)](https://github.com/KarlsruheMIS/KaMIS)
[![Homebrew](https://img.shields.io/badge/homebrew-available-orange)](https://github.com/KarlsruheMIS/homebrew-kamis)
[![GitHub Stars](https://img.shields.io/github/stars/KarlsruheMIS/KaMIS)](https://github.com/KarlsruheMIS/KaMIS/stargazers)
[![GitHub Issues](https://img.shields.io/github/issues/KarlsruheMIS/KaMIS)](https://github.com/KarlsruheMIS/KaMIS/issues)
[![GitHub Last Commit](https://img.shields.io/github/last-commit/KarlsruheMIS/KaMIS)](https://github.com/KarlsruheMIS/KaMIS/commits)
[![arXiv](https://img.shields.io/badge/arXiv-1509.00764-b31b1b)](https://arxiv.org/abs/1509.00764)
[![arXiv](https://img.shields.io/badge/arXiv-1708.06151-b31b1b)](https://arxiv.org/abs/1708.06151)
[![arXiv](https://img.shields.io/badge/arXiv-1810.10834-b31b1b)](https://arxiv.org/abs/1810.10834)
[![arXiv](https://img.shields.io/badge/arXiv-2008.05180-b31b1b)](https://arxiv.org/abs/2008.05180)
[![arXiv](https://img.shields.io/badge/arXiv-2208.13645-b31b1b)](https://arxiv.org/abs/2208.13645)
[![Agent-Ready](https://img.shields.io/badge/agent--ready-yes-brightgreen)](https://github.com/KarlsruheMIS/KaMIS)
[![Heidelberg University](https://img.shields.io/badge/Heidelberg-University-red)](https://www.uni-heidelberg.de)

## Description ##
This is the open source project KaMIS - Karlsruhe Maximum Independent Sets. Given a graph G=(V,E), the goal of the maximum independent set problem is to compute a maximum cardinality set of vertices I, such that no vertices in the set are adjacent to one another. Such a set is called a maximum independent set. The problem is NP-hard and particularly difficult to solve in large sparse graphs. 

<p align="center">
<img src="./img/isolateclique.png"
  alt="framework overview"
  width="700" height="245">
</p>


Main project site:
<http://KarlsruheMIS.github.io>

## Install via Homebrew

```console
brew install KarlsruheMIS/kamis/kamis
```

Then run directly:
```console
redumis network.graph --output independent_set.txt --time_limit 60 --console_log
```

## Installation (from source) ##
As a first step, please run *git submodule update  --init --recursive*. Then compile the source by running *compile_withcmake.sh*. The binaries can then be found in the folder *deploy*.  To compile the programs you need g++, OpenMP and cmake installed. 

To convert a graph from DIMACS to METIS format or sort its edges you can use the python scripts in the *misc* folder.

This version of our framework currently contains the following algorithms:
* struction -- a branch and reduce algorithm with increasing transformations
* mmwis -- run an iterative reduce and evolution algorithm for the weighted problem
* redumis -- run an evolutionary algorithm on a reduced graph 
* onlinemis -- local search pruned with reductions
* weighted_branch_reduce -- a branch and reduce algorithm for weighted maximum independent sets
* weighted_local_search -- a local search algorithm for weighted maximum independent sets
* If you want to use the solver that won the vertex cover track of the PACE Challenge, go here <https://github.com/KarlsruheMIS/pace-2019>

Furthermore, the framework contains tools to make life a little bit easier:
* sort_adjacencies -- takes a graph file and sorts the neighborhoods of vertices (this is required by our algorithms) 
* graphchecker -- check if the graph file you gave to algorithm is in the correct format

## NEW in v3.0: 
*mmwis (Memetic Maxmimum Weight Independent Set):* Our iterative reduce and evolution algorithm to solve the maximum weight independent set problem.

*struction:* Our new branch and reduce algorithm using increasing transformations. 


## Usage ReduMIS ##
`redumis FILE [options]`.    

### Options ###
This is a brief overview of the most important options.
For a full description, please take a look at the user guide.

`FILE`
Path to graph file that you want the maximum independent set for.

`--help`
Print help.

`--console_log`
Write the log to the console.

`--output=<string>`
Path to store the resulting independent set.

`--seed=<int>`
Seed to use for the random number generator.

`--config=<string>`
Config to use for the evolutionary algorithm [standard|social].

`--time_limit=<double>`
Time limit until the algorithm terminates.

## Usage OnlineMIS ##
`online_mis FILE [options]`.    

### Options ###
This is a brief overview of the most important options.
For a full description, please take a look at the user guide.

`FILE`
Path to graph file that you want the maximum independent set for.

`--help`
Print help.

`--console_log`
Write the log to the console.

`--output=<string>`
Path to store the resulting independent set.

`--seed=<int>`
Seed to use for the random number generator.

`--time_limit=<double>`
Time limit until the algorithm terminates.

`--adaptive_greedy`
Use adaptive greedy solution

## Usage Weighted Branch and Reduce ##
`weighted_branch_reduce FILE [options]`.    
`weighted_local_search FILE [options]`.    

### Options ###
This is a brief overview of the most important options.
For a full description, please take a look at the user guide.

`FILE`
Path to graph file that you want the maximum independent set for.

`--help`
Print help.

`--console_log`
Write the log to the console.

`--output=<string>`
Path to store the resulting independent set.

`--seed=<int>`
Seed to use for the random number generator.

`--time_limit=<double>`
Time limit until the algorithm terminates.

`--weight_source=<string>`
Choose how the weights are assigned. Can be either: file (default), hybrid, uniform, geometric.

`--reduction_style=<string>`
Choose the type of reductions appropriate for the input graph. Can be either: normal/sparse (default), dense/osm.

## Usage MMWIS ##
`mmwis FILE [options]`.    


### Options ###
This is a brief overview of the most important options.
For a full description, please take a look at the user guide.

`FILE`
Path to graph file that you want the maximum independent set for.

`--help`
Print help.

`--console_log`
Write the log to the console.

`--output=<string>`
Path to store the resulting independent set.

`--seed=<int>`
Seed to use for the random number generator.

`--config=<string>`
Config to use for the evolutionary algorithm [mmwis|mmwiss].

`--time_limit=<double>`
Time limit until the algorithm terminates.

## Usage struction ##
`struction FILE [options]`.    


### Options ###
This is a brief overview of the most important options.
For a full description, please take a look at the user guide.

`FILE`
Path to graph file that you want the maximum independent set for.

`--help`
Print help.

`--console_log`
Write the log to the console.

`--output=<string>`
Path to store the resulting independent set.

`--seed=<int>`
Seed to use for the random number generator.

`--time_limit=<double>`
Time limit until the algorithm terminates.

`--cyclicStrong`
Use cyclicStrong instead of cyclicFast (default).

## Usage Sort Adjacencies ##
`sort_adjacencies FILE`.    

The program reads a Metis file, sorts the neighborhood of each node and prints the graph to the console.

## Usage Graph Checker ##
`graphchecker FILE`.    

The program reads a Metis file and checks the file for correctness.


## License
The project is released under MIT. However, some files used for kernelization are released under the BSD 3-clause license. See the respective files for their license.
If you publish results using our algorithms, please acknowledge our work by quoting one or more of the following papers:
```
@article{DBLP:journals/heuristics/LammSSSW17,
  author    = {Sebastian Lamm and
               Peter Sanders and
               Christian Schulz and
               Darren Strash and
               Renato F. Werneck},
  title     = {Finding near-optimal independent sets at scale},
  journal   = {J. Heuristics},
  volume    = {23},
  number    = {4},
  pages     = {207--229},
  year      = {2017},
  url       = {https://doi.org/10.1007/s10732-017-9337-x},
  doi       = {10.1007/s10732-017-9337-x},
  timestamp = {Fri, 27 Dec 2019 21:13:52 +0100},
  biburl    = {https://dblp.org/rec/journals/heuristics/LammSSSW17.bib},
  bibsource = {dblp computer science bibliography, https://dblp.org}
}
```

If you use fast kernelization routines (note that this is the default), the please also cite the following:
```
@article{DBLP:journals/jea/Hespe0S19,
  author    = {Demian Hespe and
               Christian Schulz and
               Darren Strash},
  title     = {Scalable Kernelization for Maximum Independent Sets},
  journal   = {{ACM} Journal of Experimental Algorithmics},
  volume    = {24},
  number    = {1},
  pages     = {1.16:1--1.16:22},
  year      = {2019},
  url       = {https://doi.org/10.1145/3355502},
  doi       = {10.1145/3355502},
  timestamp = {Fri, 27 Mar 2020 08:38:35 +0100},
  biburl    = {https://dblp.org/rec/journals/jea/Hespe0S19.bib},
  bibsource = {dblp computer science bibliography, https://dblp.org}
}
```

If you use OnlineMIS, then please also cite the following:
```
@inproceedings{DBLP:conf/wea/DahlumLS0SW16,
  author    = {Jakob Dahlum and
               Sebastian Lamm and
               Peter Sanders and
               Christian Schulz and
               Darren Strash and
               Renato F. Werneck},
  title     = {Accelerating Local Search for the Maximum Independent Set Problem},
  booktitle = {15th International Symposium on Experimental Algorithms {SEA}},
  pages     = {118--133},
  year      = {2016},
  series    = {Lecture Notes in Computer Science},
  volume    = {9685},
  publisher = {Springer},
  url       = {https://doi.org/10.1007/978-3-319-38851-9\_9}
}
```

If you use the weighted independents set algorithms, please also cite the following: 
```
@inproceedings{DBLP:conf/alenex/Lamm0SWZ19,
        author    = {Sebastian Lamm and
                     Christian Schulz and
                     Darren Strash and
                     Robert Williger and
                     Huashuo Zhang},
        title     = {Exactly Solving the Maximum Weight Independent Set Problem on Large Real-World Graphs},
        booktitle = {Proceedings of the Twenty-First Workshop on Algorithm Engineering and Experiments, {ALENEX} 2019},
        pages     = {144--158},
        year      = {2019},
        url       = {https://doi.org/10.1137/1.9781611975499.12},
        doi       = {10.1137/1.9781611975499.12},
        publisher = {{SIAM}},
        year      = {2019}
}
```

If you use struction, please also cite the following:
```
@inproceedings{DBLP:conf/alenex/GellnerLSSZ21,
  author       = {Alexander Gellner and
                  Sebastian Lamm and
                  Christian Schulz and
                  Darren Strash and
                  Bogd{\'{a}}n Zav{\'{a}}lnij},
  title        = {Boosting Data Reduction for the Maximum Weight Independent Set Problem
                  Using Increasing Transformations},
  booktitle    = {Proceedings of the 23rd Symposium on Algorithm Engineering and Experiments,
                  {ALENEX} 2021, Virtual Conference, January 10-11, 2021},
  pages        = {128--142},
  publisher    = {{SIAM}},
  year         = {2021},
  url          = {https://doi.org/10.1137/1.9781611976472.10},
  doi          = {10.1137/1.9781611976472.10},
  biburl       = {https://dblp.org/rec/conf/alenex/GellnerLSSZ21.bib},
  bibsource    = {dblp computer science bibliography, https://dblp.org}
}
```

If you use mmwis, please also cite the following:
```
@article{DBLP:journals/jgaa/GrossmannLSS24,
  author       = {Ernestine Gro{\ss}mann and
                  Sebastian Lamm and
                  Christian Schulz and
                  Darren Strash},
  title        = {Finding Near-Optimal Weight Independent Sets at Scale},
  journal      = {J. Graph Algorithms Appl.},
  volume       = {28},
  number       = {1},
  pages        = {439--473},
  year         = {2024},
  url          = {https://doi.org/10.7155/jgaa.v28i1.2997},
  doi          = {10.7155/JGAA.V28I1.2997},
  biburl       = {https://dblp.org/rec/journals/jgaa/GrossmannLSS24.bib},
  bibsource    = {dblp computer science bibliography, https://dblp.org}
}
```