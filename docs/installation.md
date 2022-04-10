# Installation

## Including the library

You will need:
- a C++17 compliant compiler (CI verifies MSVC, Clang and GCC)
- the dependencies:
    - [Eigen](https://eigen.tuxfamily.org)
    - [XSimd](https://github.com/xtensor-stack/xsimd)

The library is header-only, therefore it's enough to put the headers in the include path.

The dependencies, Eigen and XSimd, also need to be in the include path and linked against. DSPBB does not incorporate or automatically download or find these libraries, you will have to do that yourself.

### Without dependencies

Using DSPBB without the aformentioned dependcies is on the todo list. This would mean losing some features, but making installation simpler.

### Package managers

I will put DSPBB on package management libraries as soon as it reaches that development stage.

## Compiling tests and examples

You will need:
- a C++17 compliant compiler
- CMake
- conan

Once you have `conan` in the path, just configure the main `CMakeLists.txt` and you are good to go. Be sure to set the options to compile the examples and benchmarks.