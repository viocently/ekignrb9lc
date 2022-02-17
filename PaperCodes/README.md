# Reference codes for the paper "Stretching Cube Attacks: Improved Methods to Recover Massive Superpolies"

## 1. Main components of the codes
The codes for Trivium, Grain-128AEAD, Kreyvium and Acorn are provided in the folder named by the algorithm. For example, the folder "Trivium" contains several important files:

trivium.cpp : this file defines functions for the monomial propagation of Trivium.

utrivium.cpp : this file defines functions for the NBDP propagation of Trivium.

newtrivium.cpp : this file defines functions for the core monomial propagation of Trivium, together with the MILP Model "CMP-MPModelTrivium" mentioned in the paper.

triviumtrack.cpp : this file is the implementation of our flag technique.

triviumcallback.cpp : this file shows how to construct the MILP Model "NBDP-MPModelTrivium" for the term expander. It uses the callback functions in Gurobi.

## 2. How to use codes
We only describe how to execute our code on the linux platform. Still, we take Trivium as an example.

1. First we need to install Gurobi and configure the required environment variables such as "GUROBI_HOME" and "LD_LIBRARY_PATH".

2. Switch to the folder "Trivium" and edit "main.cpp" to set the *cubeIndex* and *ROUND* you are interested in.

3. Type `sh exec.sh` in the console to generate the executable program "trivium". At the same time, three folders will be created in the current directory, namely 
"STATE", "LOG" and "TERM". "LOG" contains log files; "STATE" saves some important data so that the program can continue to run after a crash; "TERM" stores the results 
produced by the execution of the program.

4. Type `./trivium` in the console to run the program. While the program is running, you can check the status of the program by viewing the files in the LOG.

## 3. How to generate the final superpoly
We introduce how to generate the final superpoly after the program terminates. This part may be different for different algorithms.

For Kreyvium and Acorn, when the propram terminates, a file "res.txt" will be gerenated in the folder "TERM". This file contains the final superpoly.
To analyze the final superpoly, you can take the following steps:
1. Switch to the folder "Kreyvium" or "Acorn".
2. Type `g++ evaluateRes.cpp -o -std=c++11 -O3 -o evaluateRes` in the console to compile the source file "evaluateRes.cpp".
3. Type `./evaluateRes ./TERM/res.txt` in the console to start the analysis of the superpoly. When it finishes, information about the superpoly will be output.

For Trivium or Grain, to obtain the final superpoly, you can take the following steps:
1. Switch to the folder "Trivium" ("Grain").
2. Type `g++ evaluateResTrivium.cpp (evaluateResGrain.cpp) -o -std=c++17 -O3 -lpthread -o evaluateRes` in the console to compile the source file "evaluateResTrivium.cpp" ("evaluateResGrain.cpp").
3. Type `./evaluateRes ./TERM` in the console to start the analysis of the superpoly. When it finishes, a file "res.txt" containing the final superpoly will be generated in the folder "TERM". Also, information about the superpoly will be output.

Note that the header file "dynamic_bitset.hpp" included in "evaluateResTrivium.cpp" is from Boost C++ Libraires, which can be downloaded from (https://www.boost.org/).












