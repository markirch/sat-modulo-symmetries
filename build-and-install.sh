#!/bin/bash

function msg() {
	echo -e "\n-- build-and-install.sh: $1\n"
}

msg "Configuring the build directory..."
cmake -Bbuild -H. # -DGLASGOW=1

msg "building SMS"
cmake --build build -j2

msg "installing SMS binaries"
sudo cmake --install build

msg "installing Python wrapper"
pip install .

# g++ ./scripts/other/createDimacsHeader.cpp -o ./scripts/other/createDimacs
