Quantum KITE is an open-source software suite for accurate real space evaluations of electronic structure and response functions of large-scale tight-binding (TB) models of complex molecules and crystals with up to multi billions of atoms (orbitals).

KITE is written in C++ with advanced code optimization at the design level (including automated multithreaded lattice decomposition) and allow users to change defaults to enable the best possible use of resources. KITE's user-friendly interface and processing tools are based on <a href="http://docs.pybinding.site/en/stable/">Pybinding</a>, a scientific Python package for TB calculations.

KITE is also a quantum transport simulator. Multi-orbital bond disorder can be defined at the interface level and added to the system according to pre-defined probability distributions, allowing to simulate  the behavior of realistic disordered systems.  

*This is our pre-release BETA version*. The fully-debugged official release is scheduled for 01/11/2018. 

For currently available functionalities and to-do-list please refer to the README file. Please share your feedback and bug reports. We would love to hear your ideas for future improvements (contact email: support at quantum-kite.com).

#Installation KITE

KITE runs on UNIX®-based systems (including Mac OS X) and requires <a href="http://docs.pybinding.site/en/stable/">Pybinding</a>, <a href="eigen.tuxfamily.org">Eigen</a> C++ template library and <a href="https://www.hdfgroup.org/">HDF5</a> support for multi-dimensional datasets. These packages are available from public domain servers (see below).

- <a href="#linux">Linux </a>
	- <a href="#ubuntu">Ubuntu Installation</a>
	- <a href="#source">Compiling Libraries From Source Code</a>
		
- <a href="#macosx"> Mac OS X </a>	

##Linux:
<a name="linux"></a>

###Ubuntu installation
<a name="ubuntu"></a>

The instructions below were tested with Ubuntu LTS realese 16.04; for other Linux distributions please refer to *Compiling libraries from Source Code*.

In order to compile our source code, the compiler must be up to date (e.g. GCC 4.8.1 or newer). The required dependencies are:

* <a href="http://eigen.tuxfamily.org">Eigen3</a>

* <a href="https://www.hdfgroup.org/">HDF5</a> (version 1.8.13 or newer) 

* <a href="http://docs.pybinding.site/en/stable/">Pybinding</a>

Eigen (Eigen3) features various linear algebra tools. To install it, run:

~~~bash
sudo apt-get install libeigen3-dev
~~~

Hierarchical Data Format (HDF5) is used to store the inputs/outputs of the program. To install HDF5, run:

~~~bash
sudo apt-get install h5utils
sudo apt-get install libhdf5-dev
~~~

Calculations on KITE are configured using a python script which interfaces with Pybinding. To install Pybinding, you will need **CMake**  and **pip**:

~~~bash
sudo apt-get install cmake
sudo apt-get install pip3
~~~

Pybinding also requires the SciPy packages but  pip will resolve all the SciPy dependencies automatically:

~~~bash
pip3 install pybinding
~~~

Alternativelly, you might prefer to follow the instructions on <a href="http://docs.pybinding.site/en/stable/">Pybinding</a> webpages.


**After successfully installing these libraries, you will be ready to compile KITE.**


1. Fetch the source code from the Git Hub repository

~~~bash
git clone (...)
git checkout develop
~~~

2. Execute the Makefile inside the KITE folder to compile the code

~~~bash
make
~~~

3. That’s it! To compile the post-processing program run

~~~bash
cd tools/KITE
make
~~~

To generate the input file, try one of our examples. In the KITE folder, run

~~~bash
python example1.py
~~~

It creates a file names example.h5 that is used as an input for KITE:

~~~bash
./KITEx example1.h5
~~~

This first example calculates the density of states of pure graphene. To obtain the data file, you need to postprocess the output:

~~~bash
./tools/KITE/KITEx-tools example1.h5
~~~

For more details refer to the KITE [Documentation](https://quantum-kite.com/Documentation/).

###Compiling Libraries From Source Code
<a name="source"></a>

In order to compile our source code, the compiler must be up to date (e.g. GCC 4.8.1 or newer). The required dependencies are: 

* <a href="http://eigen.tuxfamily.org">Eigen3</a>

* <a href="https://www.hdfgroup.org/">HDF5</a> (version 1.8.13 or newer) 

* <a href="http://docs.pybinding.site/en/stable/">Pybinding</a>


The detailed installation steps are described on the respective websites. The GCC (g++) compiler must support C++11 features and OpenMP parallelization (version 4.8.1 or newer). To check the version installed on your computer you can use:

~~~bash
g++ --version
~~~

You can retrieve the latest stable release of Eigen3 from <a href="http://eigen.tuxfamily.org/">eigen.tuxfamily.org</a>. Unzip the file and copy the Eigen folder to /usr/include/. We suggest the installation of Pybinding using Miniconda: <a href="http://docs.pybinding.site/en/stable/install/quick.html">http://docs.pybinding.site/en/stable/install/quick.html</a>. 

##Mac OS X:
<a name="macosx"></a>

In order to install KITE on a MAC OS X system, you will need Xcode command-line tools from Apple Developer. You can download Xcode from the Apple Store. Alternatively, you can install the command tools directly from the terminal:

~~~bash
xcode-select --install
~~~

You will also need an open-source software package management system. Here, we provide detailed instructions using [Homebrew](https://brew.sh/). To install Homebrew run the following command on the terminal:

~~~bash
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
~~~

Follow the instructions provided. Next, install the C++ compiler:

~~~bash
brew install gcc@6
~~~

KITE has the following dependencies:

* <a href="http://eigen.tuxfamily.org">Eigen3</a>

* <a href="https://www.hdfgroup.org/">HDF5</a> (version 1.8.13 or newer) 

* <a href="http://docs.pybinding.site/en/stable/">Pybinding</a>

Eigen (Eigen3) features various linear algebra tools. To install it, run:

~~~bash
brew install eigen
~~~

Hierarchical Data Format (HDF5) is used to store the inputs/outputs of the program. To install HDF5, run:

~~~bash
brew install hdf5 --cc=gcc-6
~~~

Calculations on KITE are configured using a python script which interfaces with Pybinding. To install the Pybinding package, we suggest using Homebrew. (Alternatively, you can proceed with the [installation suggested by Pybinding](http://docs.pybinding.site/en/stable/install/quick.html), with the use of Miniconda.)

~~~bash
brew install cmake
~~~

~~~bash
brew install python
~~~

Last, install Pybinding with pip:

~~~bash
pip3 install pybinding
~~~
**After successfully installing these libraries, you are now ready to compile KITE.**

**IMPORTANT: the last version of matlibplot is having issues with pybinding. Until this is resolved, use:

~~~bash
pip3 install matplotlib==2.1.1 
~~~

1. Fetch the source code from the Git Hub repository

~~~bash
git clone (...)
git checkout develop
~~~

2. Execute the Makefile inside the KITE folder to compile the code

~~~bash
make
~~~

3. That’s it! To compile the post-processing program run

~~~bash
cd tools/KITE/
make
~~~

To generate the input file, try one of our examples. In the KITE folder, run

~~~bash
python example1.py
~~~

It creates a file names example.h5 that is used as an input for KITE:

~~~bash
./KITEx example1.h5
~~~

This first example calculates the density of states of pure graphene. To obtain the data file, you need to postprocess the output:

~~~bash
./tools/KITE/KITEx-tools example1.h5
~~~

For more details refer to the KITE [Documentation](https://quantum-kite.com/Documentation/).