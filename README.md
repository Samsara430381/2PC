# OPAF: Optimized Secure Two-Party Computation Protocols for Nonlinear Activation Functions in Recurrent Neural Network
## Introduction
This directory contains the code for the 2PC protocols from "OPAF: Optimized Secure Two-Party Computation Protocols for Nonlinear Activation Functions in Recurrent Neural Network".

## Setup
For setting up the code:
1. Set-up [EzPC/SCI](https://github.com/mpc-msri/EzPC/tree/master/SCI)
2. Replace the "CMakeList.txt" file in the SCI directory
3. Move the "mytests" directory into the SCI directory
4. Move the "MyProtocol" directory into the SCI/src directory

## Compilation
```
cd EzPC/SCI
mkdir mybuild && cd mybuild
cmake ..
make -j 16
```

## Running mytests
On successful compilation, the binaries will be created in `mybuild/bin/`.

Run mytests as follows to make sure everything works as intended:

`./<mytest> r=1 [port=port] & ./<mytest> r=2 [port=port]`

# Acknowledgements
Our implementation includes code from the following external repository:
 - [EzPC/SCI]([https://github.com/emp-toolkit/emp-tool/tree/c44566f40690d2f499aba4660f80223dc238eb03/emp-tool](https://github.com/mpc-msri/EzPC/tree/master/SCI)) for 
fundamental building blocks.
