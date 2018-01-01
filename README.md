# Get the latest release

You can get the source code and binaries for the latest release here: [![GitHub release](https://img.shields.io/github/release/skrakau/PureCLIP.svg)](https://github.com/skrakau/PureCLIP/releases/latest)

Alternatively, you can find instructions how to build the latest development version below.

# Build instructions

Clone the repository

    $ git clone https://github.com/skrakau/PureCLIP.git
    $ cd PureCLIP

Create a build directory, configure the build and compile

    $ mkdir build
    $ cd build
    $ cmake ../src
    $ make

Requirements

 - C++14 compliant compiler
 - GSL
 - cmake 3.0 or newer


# Documentation

Please have a look at PureCLIPs [documentation](http://pureclip.readthedocs.io/en/latest/).

# Citation

Krakau S, Richard H, Marsico A: PureCLIP: Capturing target-specific protein-RNA interaction footprints from single-nucleotide CLIP-seq data. Genome Biology 2017; 18:240; doi: https://doi.org/10.1186/s13059-017-1364-2
