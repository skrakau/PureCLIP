[![Build Status](https://travis-ci.org/skrakau/PureCLIP.svg?branch=master)](https://travis-ci.org/skrakau/PureCLIP)

# Get the latest release

You can get the source code and binaries for the latest release here: [![GitHub release](https://img.shields.io/github/release/skrakau/PureCLIP.svg)](https://github.com/skrakau/PureCLIP/releases/latest)

Alternatively, you can build the latest development version or use the Galaxy platform.

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

# Galaxy

PureCLIP has been integrated into the European Galaxy server https://usegalaxy.eu/, an open, web-based platform for accessible, reproducible, and transparent computational biological research and is available [here](https://usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/iuc/pureclip/pureclip/1.0.4).

Thanks to the Freiburg Galaxy Team!


# Documentation

Please have a look at PureCLIPs [documentation](http://pureclip.readthedocs.io/en/latest/).

# Citation

Krakau S, Richard H, Marsico A: PureCLIP: Capturing target-specific protein-RNA interaction footprints from single-nucleotide CLIP-seq data. Genome Biology 2017; 18:240; https://doi.org/10.1186/s13059-017-1364-2
