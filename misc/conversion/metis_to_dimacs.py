# Converts a METIS-graph into the DIMACS format
import sys, os


def convert_to_dimacs(filename):
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
                adjacency[0] = args
            else:
                adjacency[node] = args
                edges_counted += len(args)
            node += 1

    print("Writing new file")

    filepath = os.path.splitext(filename)
    new_file = filepath[0] + ".dimacs"

    with open(new_file, "w") as f:
        node = 0
        for neighbors in adjacency:
            if node == 0:
                spec = " ".join(neighbors)
                f.write("p " + spec)
                f.write("\n")
            else:
                for neighbor_node in neighbors:
                    if neighbor_node != "":
                        f.write("e " + str(node) + " " + neighbor_node)
                        f.write("\n")
                        if str(node) in adjacency[int(neighbor_node)]:
                            adjacency[int(neighbor_node)].remove(str(node))
            node += 1

    print("Finished converting.")


if __name__ == "__main__":
    convert_to_dimacs(sys.argv[1])
