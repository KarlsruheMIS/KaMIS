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
