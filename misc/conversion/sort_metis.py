# Converts a DIMACS-graph into the METIS format
import sys, os, re


def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    return [atoi(c) for c in re.split("(\d+)", text)]


def sort_metis(filename):
    if not os.path.isfile(filename):
        print("File not found.")
        sys.exit(0)

    number_nodes = 0
    number_edges = 0
    edges_counted = 0
    adjacency = []

    print("Reading the file.")

    with open(filename) as f:
        node = 0
        for line in f:
            args = line.strip().split()
            if node == 0:
                number_nodes = int(args[0])
                number_edges = int(args[1])
                print(
                    "Graph has "
                    + str(number_nodes)
                    + " nodes and "
                    + str(number_edges)
                    + " edges"
                )
                adjacency = [[] for _ in range(0, number_nodes + 1)]
                adjacency[0] = [str(number_nodes), str(number_edges)]
            else:
                adjacency[node] = args
                edges_counted += len(args)
            node += 1

    if edges_counted < number_edges:
        print("Found less edges than specified")
        sys.exit(0)

    filepath = os.path.splitext(filename)
    new_file = filepath[0] + "-sorted.graph"

    print("Writing new file to " + new_file)

    with open(new_file, "w") as f:
        node = 0
        for neighbors in adjacency:
            if node != 0:
                neighbors.sort(key=natural_keys)
            if not neighbors:
                f.write(" ")
            else:
                f.write(" ".join(neighbors))
            f.write("\n")
            node += 1

    print("Finished converting.")


if __name__ == "__main__":
    sort_metis(sys.argv[1])
