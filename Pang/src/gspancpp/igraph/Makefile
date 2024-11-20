igraph_test: igraph_test.o
	g++ -o igraph_test igraph_test.o -g -L/usr/local/lib -ligraph 

igraph_test.o: igraph_test.cpp
	g++ -c -g -I/usr/local/include/igraph igraph_test.cpp -o igraph_test.o

clean:
	rm igraph_test.o igraph_test


