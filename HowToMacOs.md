# MacOS Build Instructions

Here are extended instructions for how to build KaMIS from source on MacOS.

This procedure has been tested using:
- MacOS (v10.15.6)
- gcc (v10.1.0)
- KaMIS (v2.0)

1. Install Homebrew (https://brew.sh/)
   ```bash
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
   ```
2. Install gcc (which includes g++)
   Even though MacOS claims it has g++ installed the default OS call of `g++` will actually result in a call to `clang`  and _not_ g++.
   ```bash
   brew install gcc
   ```
   As of this post will install g++ v10 which install the executable as `/usr/local/bin/g++-10`

3. Install OpenMP
   ```bash
   brew install libomp
   ```
4. Ensure that `/usr/local/bin` is part of the user PATH environment variable.
   Brew installs both `g++` and the OpenMP library to `/usr/local/bin` so it needs to be discoverable by CMake during compilation. 
   This can be done by adding the following to your your `~/.bash_profile` (or whatever your shell profile is)
   ```bash
   export PATH=/usr/local/bin:$PATH
   ```
   Then executing:
   ```bash
   source ~/.bash_profile 
   ```
5. Compile KaMIS using the g++-10 compile, instead of clang.
   The default settings will have KaMIS try to compile with Clang, so instead we need to pass an extra variable to CMake to tell it to compile with g++ installed above. To do this we need to edit the base build script `~/.compile_withcmake.sh` to change the line:
   ```bash
   cmake ../
   ```
   to
   ```bash
   cmake ../ -DCMAKE_C_COMPILER=/usr/local/bin/gcc-10 -DCMAKE_CXX_COMPILER=/usr/local/bin/g++-10
   ```
6. Now you can run the base compile script as normal:
   ```bash
   ./compile_withcmake.sh
   ```
7. That’s it! You have now built KaMIS on MacOS.
