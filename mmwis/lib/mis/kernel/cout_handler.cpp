#include "cout_handler.h"

using namespace mmwis;

// Static member initialization
int cout_handler::disable_count = 0;
std::stringstream cout_handler::buffered_output;
std::streambuf* cout_handler::cout_rdbuf_backup = std::cout.rdbuf();


cout_handler::cout_handler() {};

cout_handler::~cout_handler() {};

void cout_handler::disable_cout() {
		disable_count++;

		if (disable_count == 1) {
			buffered_output.str(std::string());
			buffered_output.clear();
			std::cout.rdbuf(buffered_output.rdbuf());
		}
}

void cout_handler::enable_cout() {
		if (disable_count == 0)
			return;

		disable_count--;

		if (disable_count == 0)
			std::cout.rdbuf(cout_rdbuf_backup);
}
