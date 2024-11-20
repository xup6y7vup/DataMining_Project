#include <igraph.h>
#include <stdio.h>
#include <vector>
#include <iostream>

using namespace std;

int main(void)
{
    FILE *graphfile;
    graphfile = fopen("graphfile.dot", "w");
    vector<igraph_t> graphs;
    igraph_t graph1, graph2;
    igraph_vector_long_t g1_vcolor, g2_vcolor, g1_ecolor, g2_ecolor;
    igraph_bool_t result;
    long int color = 1;
    // create graph 1
    igraph_empty(&graph1, 0, false);
    igraph_vector_long_init(&g1_vcolor, 0);
    igraph_vector_long_init(&g1_ecolor, 0);
    igraph_add_vertices(&graph1, 1, 0); // add vertices
    igraph_vector_long_push_back(&g1_vcolor, color); // assign color
    igraph_add_vertices(&graph1, 1, 0); 
    igraph_vector_long_push_back(&g1_vcolor, color);
    igraph_add_vertices(&graph1, 1, 0); 
    igraph_vector_long_push_back(&g1_vcolor, color);
    igraph_add_vertices(&graph1, 1, 0); 
    igraph_vector_long_push_back(&g1_vcolor, color);
    igraph_add_vertices(&graph1, 1, 0); 
    igraph_vector_long_push_back(&g1_vcolor, color);
    igraph_add_edge(&graph1, 0, 1); // add edge
    igraph_vector_long_push_back(&g1_ecolor, color - 1); // assign color
    igraph_add_edge(&graph1, 1, 2);
    igraph_vector_long_push_back(&g1_ecolor, 1);
    igraph_add_edge(&graph1, 2, 3);
    igraph_vector_long_push_back(&g1_ecolor, color + 1);
    igraph_add_edge(&graph1, 3, 4);
    igraph_vector_long_push_back(&g1_ecolor, color + 2);
    igraph_add_edge(&graph1, 4, 0);
    igraph_vector_long_push_back(&g1_ecolor, color + 3);
    
    // create graph 2
    igraph_empty(&graph2, 0, false);
    igraph_vector_long_init(&g2_vcolor, 0);
    igraph_vector_long_init(&g2_ecolor, 0);
    igraph_add_vertices(&graph2, 1, 0);
    igraph_vector_long_push_back(&g2_vcolor, color);
    igraph_add_vertices(&graph2, 1, 0);
    igraph_vector_long_push_back(&g2_vcolor, color);
    igraph_add_vertices(&graph2, 1, 0);
    igraph_vector_long_push_back(&g2_vcolor, color);
    igraph_add_vertices(&graph2, 1, 0);
    igraph_vector_long_push_back(&g2_vcolor, color);
    igraph_add_vertices(&graph2, 1, 0);
    igraph_vector_long_push_back(&g2_vcolor, color);
    igraph_add_vertices(&graph2, 1, 0);
    igraph_vector_long_push_back(&g2_vcolor, color);
    igraph_add_vertices(&graph2, 1, 0);
    igraph_vector_long_push_back(&g2_vcolor, color);
    igraph_add_vertices(&graph2, 1, 0);
    igraph_vector_long_push_back(&g2_vcolor, color);
    igraph_add_vertices(&graph2, 1, 0);
    igraph_vector_long_push_back(&g2_vcolor, color);
    igraph_add_edge(&graph2, 0, 1);
    igraph_vector_long_push_back(&g2_ecolor, color - 1);
    igraph_add_edge(&graph2, 0, 2);
    igraph_vector_long_push_back(&g2_ecolor, color - 1);
    igraph_add_edge(&graph2, 0, 4);
    igraph_vector_long_push_back(&g2_ecolor, color + 3);
    igraph_add_edge(&graph2, 0, 5);
    igraph_vector_long_push_back(&g2_ecolor, color - 1);
    igraph_add_edge(&graph2, 0, 6);
    igraph_vector_long_push_back(&g2_ecolor, color - 1);
    igraph_add_edge(&graph2, 1, 2);
    igraph_vector_long_push_back(&g2_ecolor, color);
    igraph_add_edge(&graph2, 1, 7);
    igraph_vector_long_push_back(&g2_ecolor, color - 1);
    igraph_add_edge(&graph2, 2, 3);
    igraph_vector_long_push_back(&g2_ecolor, color + 1);
    igraph_add_edge(&graph2, 2, 8);
    igraph_vector_long_push_back(&g2_ecolor, color - 1);
    igraph_add_edge(&graph2, 3, 4);
    igraph_vector_long_push_back(&g2_ecolor, color + 2);
    // igraph_write_graph_dot(&graph2, graphfile);
    // fclose(graphfile);

    graphs.push_back(graph1);
    graphs.push_back(graph2);
    igraph_subisomorphic_vf2(&graphs[1], &graphs[0], &g2_vcolor, &g1_vcolor, NULL, NULL, &result, NULL, NULL);
    cout << (result ? "true" : "false") << endl;
    return 0;
}
