#!/bin/bash


if [ -t 1 ]; then
	reset="\033[0m"
	 bold="\033[1m"
	  red="\033[91m"
	green="\033[32m"
fi

msg() {
	msgcol=$green
	if [ "${2-""}" = "err" ] ; then
		msgcol=$red
	fi
	echo -e "\n${bold}--SMS-- build-and-install.sh:${reset}${msgcol} $1${reset}\n"
}


usage () {
cat <<EOF
usage: build-and-install.sh [-h|--help][-g|--debug][-s|--glasgow][-l|--local][-c|--cmake][-p|--pip]

-h | --help     print this command line option summary
-g | --debug    compile in debug mode, with assertion checking and symbols
-s | --glasgow  compile with the Glasgow subgraph solver
-l | --local    install locally for the current user (into \$HOME/.local/)
-c | --cmake   	use this cmake command (default: cmake)
-p | --pip   	use this pip command (default: pip)
EOF
}

debug=0
glasgow=0
loc_inst=0

CMAKE_CMD="cmake"
CMAKE_BUILD_DIR="build"
CMAKE_FLAGS=(-B "$CMAKE_BUILD_DIR" "-S" ".")
CMAKE_BUILD_PARALLEL_LEVEL=$(nproc --all)
export CMAKE_BUILD_PARALLEL_LEVEL
CONF_FLAGS=("-fPIC")

CADICAL_DIR="cadical_sms/"

PIP_CMD="pip"

while [ $# -gt 0 ]
do
  case $1 in
    -h|--help)    usage; exit 0;;
    -g|--debug)   debug=1;;
    -s|--glasgow) glasgow=1;;
    -l|--local)   loc_inst=1;;
    -c|--cmake)   if [ $# -eq 1 ]; then msg "expecting cmake command after '$1'" "err" && exit 1; else shift; CMAKE_CMD="$1"; fi;;
    -p|--pip)     if [ $# -eq 1 ]; then msg "expecting pip command after '$1'"   "err" && exit 1; else shift; PIP_CMD="$1"; fi;;
    *) msg "invalid option '$1' (try '-h')" "err" && exit 1;;
  esac
  shift
done

if [ $debug = 1 ]; then
	msg "Build type set to DEBUG"
	CMAKE_FLAGS+=("-DCMAKE_BUILD_TYPE=Debug")
	CONF_FLAGS+=("-g")
else
	msg "Build type set to RELEASE"
fi

if [ $loc_inst = 1 ]; then
	CMAKE_FLAGS+=("-DCMAKE_INSTALL_PREFIX=$HOME/.local")
fi

if [ ! -d "$CADICAL_DIR" ]; then
	msg "CaDiCaL expected at $CADICAL_DIR, but the directory does not exist" "err"
	exit 1
else
	msg "Building CaDiCaL"
	if [ -z "$(ls "$CADICAL_DIR")" ]; then
		echo "empty"
		git submodule init
		git submodule update
	fi
	cd "$CADICAL_DIR" || { msg "cannot cd to '$CADICAL_DIR', exiting" "err" ; exit ; }
	if [ ! -d "build" ] || [ ! -f "build/makefile" ] || [ ! -f "makefile" ]; then
		./configure "${CONF_FLAGS[@]}" || { msg "CaDiCaL 'configure' failed, exiting" "err" ; exit ; }
	fi
	make -j"$CMAKE_BUILD_PARALLEL_LEVEL" || { msg "CaDiCaL 'make' failed, exiting" "err" ; exit ; }
	cd .. || exit
fi

if [ $glasgow = 1 ]; then
	# get the Glasgow subgraph solver
	if [ ! -d "glasgow-subgraph-solver/" ]; then
		msg "Attempting to clone the Glasgow Subgraph Solver from GitHub"
		git clone https://github.com/ciaranm/glasgow-subgraph-solver
		cd glasgow-subgraph-solver || exit
	else
		cd glasgow-subgraph-solver || exit
	#	msg "Attempting to pull updates to the Glasgow Subgraph Solver from GitHub"
	#	git pull
	fi

	#CXX_VERSION_MAJOR=$(c++ --version | head -1 | awk '{print $NF}' | cut -d. -f1)
	#if ((CXX_VERSION_MAJOR < 10)); then
	#	# prior versions of compilers specify c++20 as c++2a
	#	sed -i "s/-std=c++20/-std=c++2a/" main.mk
	#fi

	msg "Building the Glasgow Subgraph Solver"
	"$CMAKE_CMD" "${CMAKE_FLAGS[@]}"
	if "$CMAKE_CMD" --build "$CMAKE_BUILD_DIR"; then
		msg "Glasgow solver built successfully"
	else
		msg "Unable to build the Glasgow solver, exiting" "err"
		exit 1
	fi

	# revert the change
	#sed -i "s/-std=c++2a/-std=c++20/" main.mk
	# climb back out
	cd - || exit

	# set CMake flags
	CMAKE_FLAGS+=("-DGLASGOW=1")
else
	msg "Building without the Glasgow Subgraph Solver"
fi


msg "Configuring SMS build directory..."
"$CMAKE_CMD" "${CMAKE_FLAGS[@]}"

msg "building SMS"
if "$CMAKE_CMD" --build "$CMAKE_BUILD_DIR"; then
	msg "SMS built successfully"
else
	msg "Unable to build SMS, exiting" "err"
	exit 1;
fi

msg "installing SMS binaries"
if [ $loc_inst -eq 1 ]; then
	"$CMAKE_CMD" --install "$CMAKE_BUILD_DIR"
else
	sudo "$CMAKE_CMD" --install "$CMAKE_BUILD_DIR"
fi

msg "installing PySMS"
"$PIP_CMD" install .
