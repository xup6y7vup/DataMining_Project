#include "gengraph_header.h"
#include "gengraph_graph_molloy_optimized.h"
#include "gengraph_graph_molloy_hash.h"
#include "gengraph_degree_sequence.h"
#include "gengraph_random.h"

#include "igraph.h"

namespace gengraph {

// return negative number if program should exit
int parse_options(int &argc, char** &argv);

// options
static bool MONITOR_TIME = false;
static int  SHUFFLE_TYPE = FINAL_HEURISTICS;
static bool RAW_DEGREES  = false;
static FILE *Fdeg = stdin;

//_________________________________________________________________________
// int main(int argc, char** argv) {

//   // options
//   SET_VERBOSE(VERBOSE_NONE);
//   if(parse_options(argc, argv) < 0) return -1;

//   //Read degree distribution
//   degree_sequence dd(Fdeg, !RAW_DEGREES);
  
//   //Allocate memory
//   if(VERBOSE()) fprintf(stderr,"Allocate memory for graph...");
//   graph_molloy_opt g(dd);
//   dd.~degree_sequence();
//   //Realize degree sequence
//   if(VERBOSE()) fprintf(stderr,"done\nRealize degree sequence...");
//   bool FAILED = !g.havelhakimi();
//   if(VERBOSE()) fprintf(stderr," %s\n", FAILED ? "Failed" : "Success");
//   if(FAILED) return 2;
//   //Merge connected components together
//   if(VERBOSE()) fprintf(stderr,"Connecting...");
//   FAILED = !g.make_connected();
//   if(VERBOSE()) fprintf(stderr," %s\n", FAILED ? "Failed" : "Success");
//   if(FAILED) return 3;
//   //Convert graph_molloy_opt to graph_molloy_hash
//   if(VERBOSE()) fprintf(stderr,"Convert adjacency lists into hash tables...");
//   int *hc = g.hard_copy();
//   g.~graph_molloy_opt();
//   graph_molloy_hash gh(hc);
//   delete[] hc;
//   if(VERBOSE()) fprintf(stderr,"Done\n");
//   //Shuffle
//   gh.shuffle(5*gh.nbarcs(), SHUFFLE_TYPE);
//   //Output
//   gh.print();
//   if(MONITOR_TIME) {
//     double t = double(clock()) / double(CLOCKS_PER_SEC);
//     fprintf(stderr,"Time used: %f\n", t);
//   }
//   return 0;
// }

//_________________________________________________________________________
int parse_options(int &argc, char** &argv) {
  bool HELP = false;
  int argc0 = argc;
  argc = 1;
  for(int a=1; a<argc0; a++) {
    if(strcmp(argv[a],"-v")==0) SET_VERBOSE(VERBOSE_SOME);
    else if(strcmp(argv[a],"-vv")==0) SET_VERBOSE(VERBOSE_LOTS);
    else if(strcmp(argv[a],"-s")==0) my_srandom(0);
    else if(strcmp(argv[a],"-?")==0 || strcmp(argv[1],"--help")==0 || strcmp(argv[1],"/?")==0) HELP = true;
    else if(strcmp(argv[a],"-t")==0) MONITOR_TIME = true;
    else if(strcmp(argv[a],"-g")==0) SHUFFLE_TYPE = GKAN_HEURISTICS;
    else if(strcmp(argv[a],"-b")==0) SHUFFLE_TYPE = BRUTE_FORCE_HEURISTICS;
    else if(strcmp(argv[a],"-f")==0) SHUFFLE_TYPE = FAB_HEURISTICS;
    else if(strcmp(argv[a],"-o")==0) SHUFFLE_TYPE = OPTIMAL_HEURISTICS;
    else if(strcmp(argv[a],"-raw")==0) RAW_DEGREES=true;
    else // No option present
      argv[argc++] = argv[a];
  }
  if(!HELP && argc==2) {
    Fdeg = fopen(argv[1],"r");
    if(Fdeg==NULL) {
      fprintf(stderr,"Error : couldn't open file \"%s\" for reading\n",argv[1]);
      return -1;
    }
    argv[1]=argv[0];
    argv++;
    argc--;
  }
  if(HELP || argc!=1) {
    fprintf(stderr,"Usage : %s [options] [file containing degree distribution]\n",argv[0]);
    fprintf(stderr," -> %s returns a graph in its standard output\n",argv[0]);
    fprintf(stderr,"    If no file is given, %s reads its standard input\n",argv[0]);
    fprintf(stderr,"    [-v] and [-vv] options causes extra verbose.\n");
    fprintf(stderr,"    [-g] option uses the Gkantsidis heuristics.\n");
    fprintf(stderr,"    [-b] option uses the Brute Force heuristics.\n");
    fprintf(stderr,"    [-f] option uses the Modified Gkantsidis heuristics.\n");
    fprintf(stderr,"    [-o] option uses the Optimal Gkantsidis heuristics.\n");
    fprintf(stderr,"    [-t] option monitors computation time\n");
    fprintf(stderr,"    [-s] does a srandom(0) to get a constant random graph\n");
    fprintf(stderr,"    [-raw] is to take raw degree sequences as input\n");
    return -1;
  }
  return 0;
}


} // namespace gengraph

using namespace gengraph;

extern "C" {

int igraph_degree_sequence_game_vl(igraph_t *graph,
				   const igraph_vector_t *out_seq,
				   const igraph_vector_t *in_seq) {

  RNG_BEGIN();

  if (in_seq && igraph_vector_size(in_seq) != 0) {
    RNG_END();
    IGRAPH_ERROR("This generator works with undirected graphs only", IGRAPH_EINVAL);
  }
  degree_sequence *dd = new degree_sequence(out_seq);
  
  graph_molloy_opt *g = new graph_molloy_opt(*dd);
  delete dd;

  if (!g->havelhakimi()) {
    delete g;
    RNG_END();
    IGRAPH_ERROR("Cannot realize the given degree sequence as an undirected, simple graph",
		 IGRAPH_EINVAL);
  }
  
  if (!g->make_connected()) {
    delete g;
    RNG_END();
    IGRAPH_ERROR("Cannot make a connected graph from the given degree sequence", 
		 IGRAPH_EINVAL);
  }
  
  int *hc = g->hard_copy();
  delete g;
  graph_molloy_hash *gh = new graph_molloy_hash(hc);
  delete [] hc;

  gh->shuffle(5*gh->nbarcs(), SHUFFLE_TYPE);
  
  IGRAPH_CHECK(gh->print(graph));
  delete gh;

  RNG_END();
  
  return 0;
}

}
