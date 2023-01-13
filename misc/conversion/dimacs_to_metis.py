# Converts a DIMACS-graph into the METIS format
import sys, os, re


def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    return [atoi(c) for c in re.split("(\d+)", text)]


def convert_to_metis(filename):
    if not os.path.isfile(filename):
        print("File not found.")
        sys.exit(0)

    number_nodes = 0
    number_edges = 0
    edges_counted = 0
    adjacency = []

    with open(filename) as f:
        for line in f:
            args = line.strip().split()
            if "cf" in args or "co" in args or "c" in args:
                continue
            if len(args) == 4:
                type, source, target, omit = args
            else:
                type, source, target = args
            if type == "p":
                number_nodes = source
                number_edges = target
                print(
                    "Graph has "
                    + number_nodes
                    + " nodes and "
                    + number_edges
                    + " edges"
                )
                adjacency = [[] for _ in range(0, int(number_nodes) + 1)]
            elif type == "e" or type == "a":
                edge_added = False
                if not target in adjacency[int(source)]:
                    adjacency[int(source)].append(target)
                    edge_added = True
                if not source in adjacency[int(target)]:
                    adjacency[int(target)].append(source)
                    edge_added = True
                if edge_added:
                    edges_counted += 1
            else:
                print("Could not read line.")

    adjacency[0].append(number_nodes)
    # adjacency[0].append(number_edges)
    adjacency[0].append(str(edges_counted))

    print("Writing new file")

    filepath = os.path.splitext(filename)
    new_file = filepath[0] + "-sorted.graph"

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
    convert_to_metis(sys.argv[1])
