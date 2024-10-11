/******************************************************************************
* cout_handler.h
*
*****************************************************************************/
#pragma once

// system includes
#include <cstddef>
#include <string>
#include <iostream>
#include <sstream>

namespace mmwis {
class cout_handler {

public:
	cout_handler();
    ~cout_handler();

	void disable_cout();
	void enable_cout();

private:
    static int disable_count;
    static std::stringstream buffered_output;
    static std::streambuf* cout_rdbuf_backup;
};
}

