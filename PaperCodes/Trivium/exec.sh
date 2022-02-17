mkdir LOG
mkdir TERM
mkdir STATE
g++  main.cpp deg.cpp trivium.cpp newtrivium.cpp triviumcallback.cpp triviumtrack.cpp triviumframework.cpp utrivium.cpp -o trivium -std=c++11 -O3 -lm -lpthread -I/$GUROBI_HOME/include/ -L/$GUROBI_HOME/lib -lgurobi_c++ -lgurobi91 -lm




