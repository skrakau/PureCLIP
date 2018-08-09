[![Build Status](https://travis-ci.org/skrakau/PureCLIP.svg?branch=master)](https://travis-ci.org/skrakau/PureCLIP) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/pureclip/README.html) [![GitHub release](https://img.shields.io/github/release/skrakau/PureCLIP.svg)](https://github.com/skrakau/PureCLIP/releases/latest)

PureCLIP is a tool to detect protein-RNA interaction footprints from single-nucleotide CLIP-seq data, such as iCLIP and eCLIP.

# Installation

You can install PureCLIP from Bioconda, from the release tarballs and from source.

## Bioconda 

Using Conda with an activated [Bioconda](http://bioconda.github.io) channel is the easiest way to install PureCLIP:

    $ conda install pureclip
    
## Release Tarballs

Alternatively, you can get the source code and binaries for macOS and Linux [here](https://github.com/skrakau/PureCLIP/releases/latest).

# Galaxy: use PureCLIP online

PureCLIP has also been integrated into the European Galaxy server https://usegalaxy.eu/, an open, web-based platform for accessible, reproducible, and transparent computational biological research and is available [here](https://usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/iuc/pureclip/pureclip/1.0.4).

Thanks to the Freiburg Galaxy Team!

# Build from source

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

Krakau S, Richard H, Marsico A: PureCLIP: Capturing target-specific protein-RNA interaction footprints from single-nucleotide CLIP-seq data. Genome Biology 2017; 18:240; https://doi.org/10.1186/s13059-017-1364-2
