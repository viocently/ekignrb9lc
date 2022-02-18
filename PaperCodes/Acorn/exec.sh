mkdir LOG
mkdir TERM
mkdir STATE
g++  main.cpp Acornframework.cpp newAcorn.cpp Acorn.cpp Acorncallback.cpp Acorntrack.cpp uAcorn.cpp  -o Acorn -std=c++11 -O3 -lm -lpthread -I/$GUROBI_HOME/include/ -L/$GUROBI_HOME/lib -lgurobi_c++ -lgurobi91 -lm