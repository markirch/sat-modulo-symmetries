#!/bin/bash

function msg() {
	echo -e "\n--SMS-- build-and-install.sh: $1\n"
}

usage () {
cat <<EOF
usage: build-and-install.sh [-h|--help][-g|--debug][-e|--expert]

-h | --help    print this command line option summary
-g | --debug   compile in debug mode, with assertion checking and symbols
-s | --glasgow compile with the Glasgow subgraph solver
-l | --local   install locally for the current user (into \$HOME/.local/)
EOF
}

debug=0
glasgow=0
loc_inst=0

while [ $# -gt 0 ]
do
  case $1 in
    -h|--help) usage; exit 0;;
    -g|--debug) debug=1;;
    -s|--glasgow) glasgow=1;;
    -l|--local) loc_inst=1;;
    *) die "invalid option '$1' (try '-h')";;
  esac
  shift
done

CMAKE_BUILD_DIR="build"
CMAKE_FLAGS="-B$CMAKE_BUILD_DIR -S."

if [ $debug = 1 ]; then
	CMAKE_FLAGS="$CMAKE_FLAGS -DCMAKE_BUILD_TYPE=Debug"
fi

if [ $loc_inst = 1 ]; then
	CMAKE_FLAGS="$CMAKE_FLAGS -DCMAKE_INSTALL_PREFIX=$HOME/.local"
fi

if [ $glasgow = 1 ]; then
	# get the Glasgow subgraph solver
	if [ ! -d "glasgow-subgraph-solver/" ]; then
		msg "Attempting to clone the Glasgow Subgraph Solver from GitHub"
		git clone https://github.com/ciaranm/glasgow-subgraph-solver
		cd glasgow-subgraph-solver
	else
		cd glasgow-subgraph-solver
		msg "Attempting to pull updates to the Glasgow Subgraph Solver from GitHub"
		git pull
	fi

	#CXX_VERSION_MAJOR=$(c++ --version | head -1 | awk '{print $NF}' | cut -d. -f1)
	#if ((CXX_VERSION_MAJOR < 10)); then
	#	# prior versions of compilers specify c++20 as c++2a
	#	sed -i "s/-std=c++20/-std=c++2a/" main.mk
	#fi

	msg "Building the Glasgow Subgraph Solver"
	cmake -Bbuild -S.
	if cmake --build build -j2; then
		msg "Glasgow solver built successfully"
	else
		msg "Unable to build the Glasgow solver, exiting"
		exit 1
	fi

	# revert the change
	#sed -i "s/-std=c++2a/-std=c++20/" main.mk
	# climb back out
	cd -

	# set CMake flags
	CMAKE_FLAGS="$CMAKE_FLAGS -DGLASGOW=1"
else
	msg "Building without the Glasgow Subgraph Solver"
fi


msg "Configuring the build directory..."
cmake $CMAKE_FLAGS

msg "building SMS"
if cmake --build "$CMAKE_BUILD_DIR" -j2; then
	msg "SMS built successfully"
else
	msg "Unable to build SMS, exiting"
	exit 1;
fi

msg "installing SMS binaries"
if [ $loc_inst = 1 ]; then
	cmake --install build
else
	sudo cmake --install build
fi

msg "installing the Python wrapper"
pip install .

