gSpanSubs -- Frequent Subgraph Pattern Mining with Substructural Constraints
add by Lei Zhao

gSpanCORK  --  Subgraph mining routine for discriminative subgraphs

gSpan implementation by Taku Kudo (from http://www.kyb.mpg.de/bs/people/nowozin/gboost/)
Conversion to feature selection procedure by Marisa Thoma

This package comes with three pre-compiled binaries in the bin/ directory
  Win 32bit: gSpanCORK.exe,
  Unix 32bit: gSpanCORK32bit, 
  Unix 64bit: gSpanCORK64bit.
In order to build the program yourself, call make in the src/ directory.
data/ contains a small example dataset from the predictive toxicology challenge.

INPUT:

  - A list of GRAPHS in DIMACS format, i.e.
    t # <graph_id>
    v <vertex_id> <vertex_label>
    ...
    e <edge_from> <edge_to> <edge_label>
    ...
    <next_graph_or_end_of_file>
    
    The vertices' indices must start at 0; all labels must be positive integers.
    
  - A class label file, row-wise containing the class labels (the first number expression per line is taken as class label) for each of the GRAPHS in the order of occurrence. Only two different classes are allowed.
  
  
OUTPUT:

  Selected Subgraphs: The selected subgraphs are printed to stdout in DIMACS format; if '-e' is set, the subgraphs are encoded as DFS codes (see also X. Yan and J. Han. "gSpan: Graph-based substructure pattern mining". In ICDM, pages 721–724, 2002):
    <forward_edge>  ::= <edge_label> (<starting_vertex_id>f<target_vertex_label>) 
    <backward_edge> ::= <edge_label> (b<target_vertex_id>) 
    <dfs_code>      ::= <starting_vertex_label> {<forward_edge> | <backward_edge>}

  Traceback: The input graphs which match the selected subgraphs can be printed out as well (option -w <where>). For <where>=1, only the matching graph ids are returned, for <where>=2, the number of occurrences of the current subgraph in the matching graphs are appended after a colon. <graph_id>:<frequency>
    

USAGE:

Call
    bin/gscork.exe -h


EXAMPLES:
  
  As an example see the ptc dataset with the graphs in data/ptc.txt and the class labels in data/ptc.lab
  
    # Get the top 10 subgraphs for a 10% frequency bound in DIMACS format 
    # without any additional information.
    bin/gSpanCORK -f CORK10 -l data/ptc.lab -m 34 < data/ptc.txt
    
    # Get the top 10 subgraphs without a frequency bound as DFS codes 
    # without any additional information.
    bin/gSpanCORK -f CORK10 -l data/ptc.lab -e < data/ptc.txt
    
    # Get all interesting subgraphs in DIMACS format with a traceback of the matching
    # input graphs for each subgraph together with its frequency in the input graph.
    bin/gSpanCORK -f CORK -l data/ptc.lab -w 2 < data/ptc.txt

    # Get all interesting subgraphs in DIMACS format for a 10% frequency bound, running
    # gSpan_CORK adapted to the input graphs which are yet to be resolved.
    # Note that the returned frequencies are <= the actual frequencies as the dataset
    # shrinks dynamically with the number of resolved subgraphs.
    bin/gSpanCORK -f FCORK -l data/ptc.lab -m 34 < data/ptc.txt
    
    # Get the top 10 subgraphs as DFS codes with additional correspondence information
    # printed to stderr.
    bin/gSpanCORK -f CORK10 -l data/ptc.lab -e -v < data/ptc.txt


This code was retrieved from http://www.dbs.ifi.lmu.de/~thoma/pub/sam2010/sam2010.zip

For citations and more information, please refer to
  Marisa Thoma, Hong Cheng, Arthur Gretton, Jiawei Han, Hans-Peter Kriegel, Alex Smola, Le Song, Philip Yu, Xifeng Yan, Karsten Borgwardt. "Discriminative frequent subgraph mining with optimality guarantees", Statistical Analysis and Data Mining (2010)

The datasets used in this paper are available under http://www.dbs.ifi.lmu.de/~thoma/pub/sam2010/data.zip  (23.4MB).

For further questions, please contact Marisa Thoma: thoma@dbs.ifi.lmu.de

  
   This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
