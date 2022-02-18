mkdir LOG
mkdir STATE
mkdir TERM
g++  main.cpp kreyviumframework.cpp kreyvium.cpp deg.cpp newkreyvium.cpp kreyviumcallback.cpp kreyviumtrack.cpp ukreyvium.cpp  -o kreyvium -std=c++11 -O3 -lm -lpthread -I/$GUROBI_HOME/include/ -L/$GUROBI_HOME/lib -lgurobi_c++ -lgurobi91 -lm