#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

int main(int argc, char* argv[]) {
	if (argc != 4) {
		std::cout << "Usage: graph_filename weights_filename output_filename" << std::endl;
		std::cout << "The program merges to files one graph and one weight file to produce a valid weighted Metis file." << std::endl;
		return 0;
	}

	std::ifstream graph_file(argv[1]);
	std::ifstream weights_file(argv[2]);
	std::ofstream output_file(argv[3]);

	if (!graph_file || !weights_file || !output_file) {
		std::cerr << "Error on opening one of the files!" << std::endl;
		return 1;
	}

	std::string graph_line, weights_line, tmp_str;

	// ignore header
	std::getline(graph_file, graph_line);

	output_file << graph_line << " 10\n";

	while (std::getline(graph_file, graph_line) && std::getline(weights_file, weights_line)) {
		std::stringstream weights_sstr;
		weights_sstr << weights_line;
		weights_sstr >> tmp_str; // ignore node ID
		weights_sstr >> tmp_str; // read weight as XX.0

		output_file << tmp_str.substr(0, tmp_str.find('.'));

		if (!graph_line.empty()) {
			output_file << " " << graph_line;
		}

		output_file << "\n";
	}

	graph_file.close();
	weights_file.close();
	output_file.close();
	
	return 0;
}
