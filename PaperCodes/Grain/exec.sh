mkdir LOG
mkdir STATE
mkdir TERM
g++  main.cpp grain.cpp  newgrain.cpp graincallback.cpp graintrack.cpp grainframework.cpp ugrain.cpp -o grain -std=c++11 -O3 -lm -lpthread -I/$GUROBI_HOME/include/ -L/$GUROBI_HOME/lib -lgurobi_c++ -lgurobi91 -lm