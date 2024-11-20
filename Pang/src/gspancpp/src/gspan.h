#ifndef GSPANDT_H
#define GSPANDT_H
/*
    gSpanSubs 2011/03/08 implemented by Lei Zhao based on:

    gspan.h, v 1.0 2009/12/21 by Marisa Thoma

    is a modification based on :

    $Id: gspan.h,v 1.6 2004/05/21 05:50:13 taku-ku Exp $;

   Copyright (C) 2004 Taku Kudo, All rights reserved.
     This is free software with ABSOLUTELY NO WARRANTY.

   This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 2 of the License, or
     (at your option) any later version.

   This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
     along with this program; if not, write to the Free Software
     Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
     02111-1307, USA

    Modified by Marisa Thoma (extracted from matlab; extended to CORK pruning)
*/

#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <list>
#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>
#include <igraph.h> // For gSpanSubs 

namespace GSPAN {

template <class T> inline void _swap (T &x, T &y) {
    T z = x;
    x = y;
    y = z;
}

struct Edge {
    int from;
    int to;
    int elabel;
    unsigned int id;
    Edge(): from(0), to(0), elabel(0), id(0) {};
};

class Vertex
{
public:
    typedef std::vector<Edge>::const_iterator const_edge_iterator;
    typedef std::vector<Edge>::iterator edge_iterator;

    int label;
    std::vector<Edge> edge;

    void push (int from, int to, int elabel)
    {
        edge.resize (edge.size()+1);
        edge[edge.size()-1].from = from;
        edge[edge.size()-1].to = to;
        edge[edge.size()-1].elabel = elabel;
        return;
    }
    Vertex() {
        label = -1;
    }
    Vertex(const Vertex & v) {
        label = v.label;
        edge.resize(v.edge.size());
        for (unsigned int i = 0; i < v.edge.size(); ++i) {
            edge[i].from = v.edge[i].from;
            edge[i].to = v.edge[i].to;
            edge[i].elabel = v.edge[i].elabel;
            edge[i].id = v.edge[i].id;
        }
    }
    Vertex& operator=(const Vertex & v) {
        if (&v != this) {
            label = v.label;
            edge.resize(v.edge.size());
            for (unsigned int i = 0; i < v.edge.size(); ++i) {
                edge[i].from = v.edge[i].from;
                edge[i].to = v.edge[i].to;
                edge[i].elabel = v.edge[i].elabel;
                edge[i].id = v.edge[i].id;
            }
        }
        return *this;
    }
};

/* For gSpanSubs: igraph structures */
struct igraphGraph { 
    igraph_t graph;
    igraph_vector_long_t vertex_labels;
    igraph_vector_long_t edge_labels;
};
typedef struct igraphGraph IGraph;
typedef std::map <unsigned int, std::map <unsigned int, std::vector <IGraph> > > igraph_map2; 

class Graph: public std::vector<Vertex> {
private:
    unsigned int edge_size_;
public:
    typedef std::vector<Vertex>::const_iterator const_vertex_iterator;
    typedef std::vector<Vertex>::iterator vertex_iterator;

    Graph (bool _directed)
    {
        directed = _directed;
    };
    Graph(const Graph & graph) {
        directed = graph.directed;
        edge_size_ = graph.edge_size_;
        resize(graph.size());
        std::copy(graph.begin(), graph.end(), begin());
    }
    Graph& operator=(const Graph & graph) {
        if (this != &graph) {
            directed = graph.directed;
            edge_size_ = graph.edge_size_;
            resize(graph.size());
            std::copy(graph.begin(), graph.end(), begin());
        }
        return *this;
    }
    bool directed;

    unsigned int edge_size ()  const {
        return edge_size_;
    }
    unsigned int vertex_size () const {
        return (unsigned int)size();    // wrapper
    }
    void buildEdge ();
    std::istream &read (std::istream &, std::string * gName = NULL); // read
    std::ostream &write (std::ostream &os) const; // write
    void check (void);

    Graph* copyGraph() const {
        Graph * graph = new Graph(directed);
        graph->edge_size_ = edge_size_;
        unsigned int i;
        for (i = 0; i < vertex_size(); ++i) {
            graph->push_back(at(i));
        }
        return graph;
    }
    IGraph get_IGraph(); // For gSpanSubs

    Graph(): edge_size_(0), directed(false) {};
};

struct PDFS {
    unsigned int id;	// ID of the original input graph
    Edge        *edge;
    PDFS        *prev;
    PDFS(): id(0), edge(0), prev(0) {};
};

class Projected: public std::vector<PDFS> {
public:
    void push (int id, Edge *edge, PDFS *prev) // graph id, edge ref, predecessor
    {
        resize (size() + 1);
        PDFS &d = (*this)[size()-1]; // == init
        d.id = id;
        d.edge = edge;
        d.prev = prev;
    }
};

class DFS {
public:
    int from;
    int to;
    int fromlabel;
    int elabel;
    int tolabel;
    unsigned int cover_subs_num;
    /* For gSpanSubs: the number of cores this DFS contains
     * also the id of next core will be checked
     */
    friend bool operator == (const DFS &d1, const DFS &d2)
    {
        return (d1.from == d2.from && d1.to == d2.to && d1.fromlabel == d2.fromlabel
                && d1.elabel == d2.elabel && d1.tolabel == d2.tolabel);
    }
    friend bool operator != (const DFS &d1, const DFS &d2) {
        return (! (d1 == d2));
    }
    DFS(): from(0), to(0), fromlabel(0), elabel(0), tolabel(0), cover_subs_num(0) {};
    friend bool operator < (const DFS &d1, const DFS &d2)
    // For gSpanSubs: map compare operator
    {
        return ( d1.from < d2.from ||
                 (d1.from == d2.from && d1.to < d2.to) ||
                 (d1.from == d2.from && d1.to == d2.to && d1.elabel < d2.elabel) ||
                 (d1.from == d2.from && d1.to == d2.to && d1.elabel == d2.elabel && d1.tolabel < d2.tolabel));
    }

};

typedef std::vector<int> RMPath;

struct DFSCode: public std::vector <DFS> {
private:
    RMPath rmpath;
public:
    /* Backtrack the right-most path of this DFS code */
    const RMPath& buildRMPath ();

    /* Convert current DFS code into a graph. */
    bool toGraph (Graph &) const;

    /* Clear current DFS code and build code from the given graph. */
    void fromGraph (Graph &g);

    /* Clear current DFS code and build code from given compressed code. */
    void fromCompressedDFSCode(std::vector<int> & compressed);

    /* Return number of nodes in the graph. */
    unsigned int nodeCount (void) const;
    unsigned int edgeCount (void) const; // For gSpanSubs

    void push (int from, int to, int fromlabel, int elabel, int tolabel)
    {
        resize (size() + 1);
        DFS &d = (*this)[size()-1];

        d.from = from;
        d.to = to;
        d.fromlabel = fromlabel;
        d.elabel = elabel;
        d.tolabel = tolabel;
    }
    void pop () {
        resize (size()-1);
    }
    std::ostream &write (std::ostream &); // write
    std::vector<int> getCompressedDFSCode() const;

    /* For gSpanSubs */  
    bool is_min ();  // the minimality test
    // function is originally in the scope of gSpan
    bool project_is_min (Projected &, Graph &, DFSCode &); 
    // function is originally in the scope of gSpan
    DFSCode get_min_DFS();
    // compute minimum DFS code for core
    bool project_min_DFS(Projected &, Graph &, DFSCode &);
    // recursively compute minimum DFS code for core
    IGraph get_IGraph(); // get igraph structure
    bool cover_all_subs(std::map< unsigned int, IGraph > &); 
    // in multiple constraints, check if contains all cores
};

/* For gSpanSubs: tree of all discovered patterns' minimum DFS codes 
 * used in using minimum DFS code to eliminate duplicate patterns
 */
struct DFSCodeTree
{
    DFS curDFS;
    bool isPattern;
    std::map<DFS, DFSCodeTree *> children;
public:
    bool insert(DFSCode);
    void show(unsigned int indent);
    DFSCodeTree(): isPattern(false) {};
};

class History: public std::vector<Edge*>
{
private:
    std::vector<int> edge;
    std::vector<int> vertex;
    std::map<int, int> vertex_rank; // For gSpanSubs: vertex -> id
    std::vector<int> sorted_vertex; // For gSpanSubs: id -> vertex

public:
    Edge * lastFwdEdge; // For gSpanSubs 
    bool hasEdge   (unsigned int id) {
        return (bool)edge[id];
    }
    bool hasVertex (unsigned int id) {
        return (bool)vertex[id];
    }
    void build     (Graph &, PDFS *);
    void buildSubsInfo (); // For gSpanSubs: extract information for gSpanSubs pattern growth
    History() {};
    History (Graph& g, PDFS *p) {
        build (g, p);
    }

    unsigned int get_vertex_rank(int id)
    // For gSpanSubs 
    {
        return vertex_rank[id];
    }

    unsigned int get_vertex(int rank)
    // For gSpanSubs 
    {
        return sorted_vertex[rank];
    }

    unsigned int nodeCount()
    // For gSpanSubs 
    {
        return sorted_vertex.size();
    }
};

typedef std::vector <Edge*> EdgeList;

bool  get_forward_pure   (Graph&, Edge *,  int,    History&, EdgeList &);
bool  get_forward_rmpath (Graph&, Edge *,  int,    History&, EdgeList &);
bool  get_forward_root   (Graph&, Vertex&, EdgeList &);
Edge *get_backward       (Graph&, Edge *,  Edge *, History&);

/* For gSpanSubs*/
bool  get_forward_root_subs (Graph&, Vertex&, EdgeList &, History); 
bool  get_backward_root_subs (Graph&, Vertex&, EdgeList &, History); 
bool get_forward_subs (Graph&, Edge *, unsigned int, History&, EdgeList &); 
bool get_backward_subs (Graph&, Edge *, unsigned int, History&, EdgeList &); 

  /* Strucutre for maintaining correspondence counts and possible extensions.
   */
struct Correspondence {
    int father; // fathering correspondence id
    std::vector<unsigned int> counts; // # matching instances for each class (C_x)
    std::vector<unsigned int> ext; // matching extension (for a super-pattern) instances for each class (C_x_1)
};

/* For gSpanSubs: these definitions are originally in the scope of gSpan */
typedef std::map<int, std::map <int, std::map <int, Projected> > >           Projected_map3; // from (vlab) -> via (elab) -> to (vlab)
typedef std::map<int, std::map <int, Projected> >                            Projected_map2;
typedef std::map<int, Projected>                                             Projected_map1;
typedef std::map<int, std::map <int, std::map <int, Projected> > >::iterator Projected_iterator3;
typedef std::map<int, std::map <int, Projected> >::iterator                  Projected_iterator2;
typedef std::map<int, Projected>::iterator                                   Projected_iterator1;
typedef std::map<int, std::map <int, std::map <int, Projected> > >::reverse_iterator Projected_riterator3;

typedef std::pair<int, int> swapLabelPair;
// For gSpanSubs

class gSpan {

private:
    std::vector < Graph >       TRANS; // the graph dataset
    DFSCode                     DFS_CODE; // the currently examined DFS code
    bool verbose;
    unsigned int ID;  // the currently mined subgraph's id
    unsigned int minsup;  // the minimum frequency bound
    unsigned int maxpat_min;	// lower bound on node count
    unsigned int maxpat_max;	// upper bound on node count
    unsigned int where; // traceback option:
    // 0 => no traceback, 1 => traceback with graph ids,
    // 2 => traceback with graph ids and frequencies
    bool enc;  // encode output graphs as DFS codes?
    bool xml;  // display output as xml?
    bool directed;  // use directed edges?
    std::ostream* os;  // output stream
    std::ostream* infoStream;  // info stream for statistical output

    /* Singular vertex handling stuff; not supported by gSpan_CORK
     * [graph][vertexlabel] = count.
     */
    std::map<unsigned int, std::map<unsigned int, unsigned int> > singleVertex;
    std::map<unsigned int, unsigned int> singleVertexLabel;
    void report_single (Graph &g, std::map<unsigned int, unsigned int>& ncount);

    // transforms the matches in projected to support counts
    std::map<unsigned int, unsigned int> support_counts (Projected &projected);
    unsigned int support (Projected &);
    void project         (Projected &);
    void report          (Projected &, unsigned int);
    // read TRANS from is, optionally using an id filter
    std::istream &read (std::istream & is, std::set<unsigned int> const * filter = NULL);
    void run_intern (void);
    /* Nested Feature Selection variables */
    std::string fs_option;
    bool fs; // perform CORK-filtering?
    bool do_le; // accept new subgraphs even if they do not improve CORK
    // needed for reaching a fixed number of selected features
    unsigned int number_of_classes;
    /*    graph vectors (corresponding to TRANS)     */
    std::vector<unsigned int> class_labels;
    std::vector<unsigned int> equ_classes; // equivalence classes
    std::vector<bool> possible_changes; // equ_classes may be split for these graphs (= subgraph(s) found)
    std::vector<unsigned int> originalInstanceLabels; // initialized if a graph file has not been read as a whole, but filtered
    std::vector<unsigned int> originalInstanceRanks;
    bool delete_resolved; // delete graphs (and their corresponding vectors) if they are no longer part of a correspondence
    // graph name => class label  (for creating the class_labels mapping)
    std::map<unsigned int,unsigned int> instances2classLabels;

    // all unresolved correspondences classes for the current feature choice
    std::vector<Correspondence> correspondence_classes;
    DFSCode current_best_subgraph;
    std::map<unsigned int, unsigned int> current_best_support_counts;
    std::vector<DFSCode> selectedSubgraphs; // DFS codes of selected subgraphs
    unsigned int correspondences; // for the current set of selected subgraphs
    // statistics for a noisy output:
    unsigned int tested_subgraphs; // counts the number of subgraphs tested in one iteration of CORK
    unsigned int frequent_subgraphs; // counts the number of tested frequent subgraphs in one iteration of CORK
    unsigned int minimal_subgraphs; // counts the number of tested minimal subgraphs in one iteration of CORK
    unsigned int winner_subgraphs; // counts the number of subgraphs replacing the current best subgraph in one iteration of CORK
    unsigned int prunedFS_subgraphs; // counts the number of subgraphs pruned by CORK in one iteration of CORK
    // these statistics are streamed to infoStream

    // translates the graph ids to their original ids
    std::map<unsigned int,unsigned int> translateSupportCounts(std::map<unsigned int,unsigned int> const & _instances2classLabels);
    // create the class label mapping for TRANS
    std::vector<unsigned int> assign_class_labels(std::string class_label_file, std::set<unsigned int> const * graphFilter = NULL);
    // output subgraph toReport using the given support counts
    void report(DFSCode & toReport, std::map<unsigned int, unsigned int> & supp_counts);

    // reset the extension fields in correspondence_classes
    void resetCORK_Extensions();
    // extend correspondence_classes using their extension fields
    void extendCORK();
    // get the CORK value of the current subgraph and decide whether or not its children must be examined
    bool getAskCORK(Projected const * projected, unsigned int * get_cork_value = NULL, bool use_conserved = false);

    /* For gSpanSubs */
    clock_t mining_time, edge_grow_time, filter_sup_time, filter_dup_time, filter_subs_time, relabel_time, start_time, gspan_filter_sup_time, gspan_filter_dup_time, gspan_mining_time, gspan_edge_grow_time; 
    bool report_result; // for performance test
    unsigned int pattern_num; 
    Graph BASE_SUBS_GRAPH; // the base for multiple constraints
    DFSCode SUBS_MINDFS; // minimum DFS code for base or core
    std::vector <Graph> SUBS_SET; // set of all cores
    std::map <unsigned int, IGraph > SUBS_IGRAPH_SET; // igraph structures for all cores
    Projected SUBS_APPEARANCE; // appearance of core or base
    DFSCodeTree minDFSTree; // minimum DFS code tree structure for all discovered pattern
    unsigned int BASE_SUBS_GRAPH_ID;
    igraph_map2 FREQ_PATTERNS; // frequent patterns with the same edge/vertex size
    void read_SUBS_GRAPH(std::string); 
    void get_BASE_SUBS_MINDFS(void); 
    bool project_SUBS_MINDFS (Projected &, Projected &); 
    void project_subs(Projected &); 
    void project_subs_get_root(Projected &); 
    bool insert_pattern(IGraph); // insert pattern to all discovered patterns
    Projected ask_all_subs(Projected);
    // filter from projected those which does not contain all subs
    void relabel(swapLabelPair &, swapLabelPair &);

public:
    gSpan (void);

    // run gSpan only
    void run (std::istream &is, std::ostream &_os,
              unsigned int _minsup,
              unsigned int _maxpat_min, unsigned int _maxpat_max,
              bool _enc,
              unsigned int _where,
              bool _xml,
              bool _directed, std::set<unsigned int> const * graphFilter,
              bool _report_result = true, // For gSpanSubs: performance
              bool _verbose = false);

    // run gSpan_CORK giving a class label file
    std::vector<DFSCode> run_gSpan_CORK (std::istream &is, std::ostream &_os,
                                         unsigned int _minsup,
                                         unsigned int _maxpat_min, unsigned int _maxpat_max,
                                         bool _enc, unsigned int _where,
                                         bool _xml, bool _directed,
                                         std::set<unsigned int> const * graphFilter,
                                         std::string _fs_option = "",
                                         std::string class_label_file = "",
                                         std::vector<std::map<unsigned int, unsigned int> > * supportMap = NULL,
                                         std::string subs_file = "", //  
                                         bool _report_result = true, // For gSpanSubs
                                         bool _verbose = false);

    // run gSpan_CORK giving a class label map
    std::vector<DFSCode> run_gSpan_CORK (std::istream &is, std::ostream &_os,
                                         unsigned int _minsup,
                                         unsigned int _maxpat_min, unsigned int _maxpat_max,
                                         bool _enc, unsigned int _where,
                                         bool _xml, bool _directed,
                                         std::set<unsigned int> const * graphFilter,
                                         std::string _fs_option = "",
                                         std::map<unsigned int,unsigned int> const * _instances2classLabels = NULL,
                                         std::vector<std::map<unsigned int, unsigned int> > * supportMap = NULL,
                                         std::string subs_file = "", // For gSpanSubs
                                         bool _report_result = true, // For gSpanSubs
                                         bool _verbose = false);

    void setInfoStream(std::ostream &infoS);
};
};
#endif // GSPANDT_H
