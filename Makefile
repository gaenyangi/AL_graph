all: compile

compile: graph.cpp
	g++ graph.cpp -o graph

run: graph
	./graph input_sp.txt input_mst.txt output_sp.txt output_mst.txt
clean: graph
	rm graph
