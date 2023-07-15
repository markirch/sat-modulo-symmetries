#!/bin/bash
set -e

function msg() {
	echo -e "\n-- build-and-install-local.sh: $1\n"
}

msg "Configuring the build directory..."
cmake -Bbuild -H. -DCMAKE_INSTALL_PREFIX=$HOME/.local

msg "building SMS"
cmake --build build -j 2

msg "installing SMS binaries"
cmake --install build

msg "installing Python wrapper"
pip install .

# g++ ./scripts/other/createDimacsHeader.cpp -o ./scripts/other/createDimacs
