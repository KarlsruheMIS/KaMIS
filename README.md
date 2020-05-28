# KaMIS #
## Description ##
Framework for finding high-quality independent sets of large sparse graphs.

Main project site:
http://algo2.iti.kit.edu/kamis/

## Installation ##
Compile the source by running *compile_withcmake.sh*. The binaries can then be found in the folder *deploy*.
To compile the programs you need to have Argtable, g++, OpenMP and cmake installed.

To convert a graph from DIMACS to METIS format or sort its edges you can use the python scripts in the *misc* folder.

## Usage ##
`redumis FILE [options]`.    

### Options ###
This is a brief overview of the most important options.
For a full description, please take a look at the user guide.

`FILE`
Path to graph file that you want the maximum independent set for.

`-help`
Print help.

`-console_log`
Write the log to the console.

`-output=<string>`
Path to store the resulting independent set.

`-seed=<int>`
Seed to use for the random number generator.

`-config=<string>`
Config to use for the evolutionary algorithm [standard|social].

`-time_limit=<double>`
Time limit until the algorithm terminates.

## License
The project is released under MIT. However, some files used for kernelization are released under the BSD 3-clause license. See the respective files for their license.
If you publish results using our algorithms, please acknowledge our work by quoting one or more of the following papers:
```
@inproceedings{kamis2015,
             AUTHOR = {Lamm, S. and Sanders, P. and Schulz, C. and Strash, D. and Werneck, R. F.},
             TITLE = {{Finding Near-Optimal Independent Sets at Scale}},
             BOOKTITLE = {18th Meeting on Algorithm Engineering and Exerpimentation (ALENEX'16)},
             YEAR = {2016}
}
```

If you use parallel kernelization routines (note that this is the default), the please also cite the following:
```
@article{kamis2019JV,
             AUTHOR = {Hespe, D. and Schulz, C. and Strash, D.},
             TITLE = {{Scalable Kernelization for Maximum Independent Sets}},
             JOURNAL = {Journal of Experimental Algorithms},
             PUBLISHER = {ACM},
             PAGES = {1.16:1--1.16:22},
             VOLUME = {24},
             NUMBER = {1},
             YEAR = {2019}
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
