/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2006  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "igraph.h"
#include "memory.h"
#include "random.h"

#include <string.h>

/**
 * \function igraph_motifs_randesu
 * \brief Count the number of motifs in a graph
 * 
 * </para><para>
 * Motifs are small subgraphs of a given structure in a graph. It is
 * argued that the motif profile (ie. the number of different motifs
 * in the graph) is characteristic for different types of networks and
 * network function is related to the motifs in the graph. 
 * 
 * </para><para>
 * This function is able to find the different motifs of size three
 * and four (ie. the number of different subgraphs with three and four
 * vertices) in the network. (This limitation is the result of the
 * lack of code to decide graph isomorphism for larger graphs.)
 * 
 * </para><para>
 * In a big network the total number of motifs can be very large, so
 * it takes a lot of time to find all of them, a sampling method can
 * be used. This function is capable of doing sampling via the
 * \c cut_prob argument. This argument gives the probability that
 * a branch of the motif search tree will not be explored. See
 * S. Wernicke and F. Rasche: FANMOD: a tool for fast network motif
 * detection, Bioinformatics 22(9), 1152--1153, 2006 for details.
 * 
 * </para><para>
 * Set the \c cut_prob argument to a zero vector for finding all
 * motifs. 
 * 
 * </para><para> 
 * Directed motifs will be counted in directed graphs and undirected
 * motifs in undirected graphs.
 *
 * \param graph The graph to find the motifs in.
 * \param hist The result of the computation, it gives the number of
 *        motifs found for each isomorphism class. See
 *        \ref igraph_isoclass() for help about isomorphism classes.
 * \param size The size of the motifs to search for. Only three and
 *        four are implemented currently. The limitation is not in the
 *        motif finding code, but the graph isomorphism code.
 * \param cut_prob Vector of probabilities for cutting the search tree
 *        at a given level. The first element is the first level, etc.
 *        Supply all zeros here (of length \c size) to find all motifs 
 *        in a graph.
 * \return Error code.
 * \sa \ref igraph_motifs_randesu_estimate() for estimating the number
 * of motifs in a graph, this can help to set the \c cut_prob
 * parameter; \ref igraph_motifs_randesu_no() to calculate the total
 * number of motifs of a given size in a graph.
 * 
 * Time complexity: TODO.
 */

int igraph_motifs_randesu(const igraph_t *graph, igraph_vector_t *hist, 
			  int size, const igraph_vector_t *cut_prob) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_adjlist_t allneis, alloutneis;
  igraph_vector_t *neis;
  long int father;
  long int i, j, s;
  long int motifs=0;

  igraph_vector_t vids;		/* this is G */
  igraph_vector_t adjverts;	/* this is V_E */
  igraph_stack_t stack;		/* this is S */
  long int *added;
  char *subg;
  
  long int histlen;
  unsigned int *arr_idx, *arr_code;
  int code=0;
  unsigned char mul, idx;

  if (size != 3 && size != 4) {
    IGRAPH_ERROR("Only 3 and 4 vertex motifs are implemented",
		 IGRAPH_EINVAL);
  }
  if (size==3) {
    mul=3;
    if (igraph_is_directed(graph)) {
      histlen=16;
      arr_idx=igraph_i_isoclass_3_idx;
      arr_code=igraph_i_isoclass2_3;
    } else {
      histlen=4;
      arr_idx=igraph_i_isoclass_3u_idx;
      arr_code=igraph_i_isoclass2_3u;
    }
  } else {
    mul=4;
    if (igraph_is_directed(graph)) {
      histlen=218;
      arr_idx=igraph_i_isoclass_4_idx;
      arr_code=igraph_i_isoclass2_4;
    } else {
      histlen=11;
      arr_idx=igraph_i_isoclass_4u_idx;
      arr_code=igraph_i_isoclass2_4u;
    }
  }

  IGRAPH_CHECK(igraph_vector_resize(hist, histlen));
  igraph_vector_null(hist);
  
  added=igraph_Calloc(no_of_nodes, long int);
  if (added==0) {
    IGRAPH_ERROR("Cannot find motifs", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, added);

  subg=igraph_Calloc(no_of_nodes, char);
  if (subg==0) {
    IGRAPH_ERROR("Cannot find motifs", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, subg);

  IGRAPH_CHECK(igraph_adjlist_init(graph, &allneis, IGRAPH_ALL));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &allneis);  
  IGRAPH_CHECK(igraph_adjlist_init(graph, &alloutneis, IGRAPH_OUT));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &alloutneis);  

  IGRAPH_VECTOR_INIT_FINALLY(&vids, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&adjverts, 0);
  IGRAPH_CHECK(igraph_stack_init(&stack, 0));
  IGRAPH_FINALLY(igraph_stack_destroy, &stack);

  RNG_BEGIN();

  for (father=0; father<no_of_nodes; father++) {
    long int level;

    IGRAPH_ALLOW_INTERRUPTION();
    
    /* init G */
    igraph_vector_clear(&vids); level=0;
    IGRAPH_CHECK(igraph_vector_push_back(&vids, father));
    subg[father]=1; added[father] += 1; level += 1;
    
    /* init V_E */
    igraph_vector_clear(&adjverts);
    neis=igraph_adjlist_get(&allneis, father);
    s=igraph_vector_size(neis);
    for (i=0; i<s; i++) {
      long int nei=VECTOR(*neis)[i];
      if (!added[nei] && nei > father) {
	IGRAPH_CHECK(igraph_vector_push_back(&adjverts, nei));
	IGRAPH_CHECK(igraph_vector_push_back(&adjverts, father));
      }
      added[nei] += 1;
    }
    
    /* init S */
    igraph_stack_clear(&stack);

    while (level > 1 || !igraph_vector_empty(&adjverts)) {
      igraph_real_t cp=VECTOR(*cut_prob)[level];

      if (level==size-1) {
	s=igraph_vector_size(&adjverts)/2;
	for (i=0; i<s; i++) {
	  long int k, s2;
	  long int last;

	  if (cp!=0 && RNG_UNIF01() < cp) { continue; }
	  motifs+=1;
	  
	  last=VECTOR(adjverts)[2*i];
	  IGRAPH_CHECK(igraph_vector_push_back(&vids, last));
	  subg[last]=size;

	  code=0; idx=0;
	  for (k=0; k<size; k++) {
	    long int from=VECTOR(vids)[k];
 	    neis=igraph_adjlist_get(&alloutneis, from);
	    s2=igraph_vector_size(neis);
	    for (j=0; j<s2; j++) {
	      long int nei=VECTOR(*neis)[j];
	      if (subg[nei] && k != subg[nei]-1) {
		idx=mul*k+(subg[nei]-1);
		code |= arr_idx[idx];
	      }
	    }
	  }

	  igraph_vector_pop_back(&vids);
	  subg[last]=0;
	  VECTOR(*hist)[arr_code[code]] += 1;
	}
      }

      /* can we step down? */
      if (level < size-1 && 
	  !igraph_vector_empty(&adjverts) && 
	  (cp==0 || RNG_UNIF01() > cp)) {
	/* yes, step down */
	long int neifather=igraph_vector_pop_back(&adjverts);
	long int nei=igraph_vector_pop_back(&adjverts);
	IGRAPH_CHECK(igraph_vector_push_back(&vids, nei));
	subg[nei] = level+1; added[nei] += 1; level += 1;

	IGRAPH_CHECK(igraph_stack_push(&stack, neifather));
	IGRAPH_CHECK(igraph_stack_push(&stack, nei));
	IGRAPH_CHECK(igraph_stack_push(&stack, level));
	
	neis=igraph_adjlist_get(&allneis, nei);
	s=igraph_vector_size(neis);
	for (i=0; i<s; i++) {
	  long int nei2=VECTOR(*neis)[i];
	  if (!added[nei2] && nei2 > father) {
	    IGRAPH_CHECK(igraph_vector_push_back(&adjverts, nei2));
	    IGRAPH_CHECK(igraph_vector_push_back(&adjverts, nei));
	  }
	  added[nei2] += 1;
	}
      } else {
	/* no, step back */
	long int nei, neifather;
	while (!igraph_stack_empty(&stack) &&
	       level==igraph_stack_top(&stack)-1) {
	  igraph_stack_pop(&stack);
	  nei=igraph_stack_pop(&stack);
	  neifather=igraph_stack_pop(&stack);
	  igraph_vector_push_back(&adjverts, nei);
	  igraph_vector_push_back(&adjverts, neifather);
	}

	nei=igraph_vector_pop_back(&vids);
	subg[nei]=0; added[nei] -= 1; level -= 1;
	neis=igraph_adjlist_get(&allneis, nei);
	s=igraph_vector_size(neis);
	for (i=0; i<s; i++) {
	  added[ (long int) VECTOR(*neis)[i] ] -= 1;
	}
	while (!igraph_vector_empty(&adjverts) && 
	       igraph_vector_tail(&adjverts)==nei) {
	  igraph_vector_pop_back(&adjverts);
	  igraph_vector_pop_back(&adjverts);
	}
      }
      
    } /* while */

    /* clear the added vector */
    added[father] -= 1;
    subg[father] = 0;
    neis=igraph_adjlist_get(&allneis, father);
    s=igraph_vector_size(neis);
    for (i=0; i<s; i++) {
      added[ (long int) VECTOR(*neis)[i] ] -= 1;
    }

  } /* for father */

  RNG_END();

  igraph_Free(added);
  igraph_Free(subg);
  igraph_vector_destroy(&vids);
  igraph_vector_destroy(&adjverts);
  igraph_adjlist_destroy(&alloutneis);
  igraph_adjlist_destroy(&allneis);
  igraph_stack_destroy(&stack);
  IGRAPH_FINALLY_CLEAN(7);
  return 0;
}

/**
 * \function igraph_motifs_randesu_estimate
 * \brief Estimate the total number of motifs in a graph
 * 
 * </para><para>
 * This function is useful for large graphs for which it is not
 * feasible to count all the different motifs, because there is very
 * many of them.
 *
 * </para><para>
 * The total number of motifs is estimated by taking a sample of
 * vertices and counts all motifs in which these vertices are
 * included. (There is also a \c cut_prob parameter which gives the
 * probabilities to cut a branch of the search tree.)
 *
 * </para><para> 
 * Directed motifs will be counted in directed graphs and undirected
 * motifs in undirected graphs.
 *
 * \param graph The graph object to study.
 * \param est Pointer to an integer type, the result will be stored
 *        here.
 * \param size The size of the motif to look for.
 * \param cut_prob Vector giving the probabilities to cut a branch of
 *        the search tree and omit counting the motifs in that branch.
 *        It contains a probability for each level. Supply \c size
 *        zeros here to count all the motifs in the sample.
 * \param sample_size The number of vertices to use as the
 *        sample. This parameter is only used if the \c parsample
 *        argument is a null pointer.
 * \param parsample Either pointer to an initialized vector or a null
 *        pointer. If a vector then the vertex ids in the vector are
 *        used as a sample. If a null pointer then the \c sample_size
 *        argument is used to create a sample of vertices drawn with
 *        uniform probability.
 * \return Error code.
 * \sa \ref igraph_motifs_randesu(), \ref igraph_motifs_randesu_no().
 * 
 * Time complexity: TODO.
 */

int igraph_motifs_randesu_estimate(const igraph_t *graph, igraph_integer_t *est,
				   int size, const igraph_vector_t *cut_prob, 
				   igraph_integer_t sample_size, 
				   const igraph_vector_t *parsample) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_t neis;
    
  igraph_vector_t vids;		/* this is G */
  igraph_vector_t adjverts;	/* this is V_E */
  igraph_stack_t stack;		/* this is S */
  long int *added;
  igraph_vector_t *sample;
  long int sam;
  long int i;

  added=igraph_Calloc(no_of_nodes, long int);
  if (added==0) {
    IGRAPH_ERROR("Cannot find motifs", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, added);

  IGRAPH_VECTOR_INIT_FINALLY(&vids, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&adjverts, 0);
  IGRAPH_CHECK(igraph_stack_init(&stack, 0));
  IGRAPH_FINALLY(igraph_stack_destroy, &stack);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);

  if (parsample==0) {
    sample=igraph_Calloc(1, igraph_vector_t);
    if (sample==0) {
      IGRAPH_ERROR("Cannot estimate motifs", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, sample);
    IGRAPH_VECTOR_INIT_FINALLY(sample, 0);  
    IGRAPH_CHECK(igraph_random_sample(sample, 0, no_of_nodes-1, sample_size));
  } else {
    sample=(igraph_vector_t*)parsample;
    sample_size=igraph_vector_size(sample);
  }

  *est=0;

  RNG_BEGIN();

  for (sam=0; sam<sample_size; sam++) {
    long int father=VECTOR(*sample)[sam];
    long int level, s;

    IGRAPH_ALLOW_INTERRUPTION();

    /* init G */
    igraph_vector_clear(&vids); level=0;
    IGRAPH_CHECK(igraph_vector_push_back(&vids, father));
    added[father] += 1; level += 1;
    
    /* init V_E */
    igraph_vector_clear(&adjverts);
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, father, IGRAPH_ALL));
    s=igraph_vector_size(&neis);
    for (i=0; i<s; i++) {
      long int nei=VECTOR(neis)[i];
      if (!added[nei]) {
	IGRAPH_CHECK(igraph_vector_push_back(&adjverts, nei));
	IGRAPH_CHECK(igraph_vector_push_back(&adjverts, father));
      }
      added[nei] += 1;
    }
    
    /* init S */
    igraph_stack_clear(&stack);

    while (level > 1 || !igraph_vector_empty(&adjverts)) {
      igraph_real_t cp=VECTOR(*cut_prob)[level];
      
      if (level==size-1) {
	s=igraph_vector_size(&adjverts)/2;
	for (i=0; i<s; i++) {
	  if (cp!=0 && RNG_UNIF01() < cp) { continue; }
	  (*est) += 1;
	}
      }

      if (level < size-1 && 
	  !igraph_vector_empty(&adjverts) && 
	  (cp==0 || RNG_UNIF01() > cp)) {
	/* yes, step down */
	long int neifather=igraph_vector_pop_back(&adjverts);
	long int nei=igraph_vector_pop_back(&adjverts);
	IGRAPH_CHECK(igraph_vector_push_back(&vids, nei));
	added[nei] += 1; level += 1;

	IGRAPH_CHECK(igraph_stack_push(&stack, neifather));
	IGRAPH_CHECK(igraph_stack_push(&stack, nei));
	IGRAPH_CHECK(igraph_stack_push(&stack, level));

	IGRAPH_CHECK(igraph_neighbors(graph, &neis, nei, IGRAPH_ALL));
	s=igraph_vector_size(&neis);
	for (i=0; i<s; i++) {
	  long int nei2=VECTOR(neis)[i];
	  if (!added[nei2]) {
	    IGRAPH_CHECK(igraph_vector_push_back(&adjverts, nei2));
	    IGRAPH_CHECK(igraph_vector_push_back(&adjverts, nei));
	  }
	  added[nei2] += 1;
	}
      } else {
	/* no, step back */
	long int nei, neifather;
	while (!igraph_stack_empty(&stack) &&
	       level==igraph_stack_top(&stack)-1) {
	  igraph_stack_pop(&stack);
	  nei=igraph_stack_pop(&stack);
	  neifather=igraph_stack_pop(&stack);
	  igraph_vector_push_back(&adjverts, nei);
	  igraph_vector_push_back(&adjverts, neifather);
	}

	nei=igraph_vector_pop_back(&vids);
	added[nei] -= 1; level -= 1;
	IGRAPH_CHECK(igraph_neighbors(graph, &neis, nei, IGRAPH_ALL));
	s=igraph_vector_size(&neis);
	for (i=0; i<s; i++) {
	  added[ (long int) VECTOR(neis)[i] ] -= 1;
	}
	while (!igraph_vector_empty(&adjverts) && 
	       igraph_vector_tail(&adjverts)==nei) {
	  igraph_vector_pop_back(&adjverts);
	  igraph_vector_pop_back(&adjverts);
	}
      }
      
    } /* while */

    /* clear the added vector */
    added[father] -= 1;
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, father, IGRAPH_ALL));
    s=igraph_vector_size(&neis);
    for (i=0; i<s; i++) {
      added[ (long int) VECTOR(neis)[i] ] -= 1;
    }

  } /* for father */

  RNG_END();

  (*est) *= ((double)no_of_nodes/sample_size);
  (*est) /= size;

  if (parsample==0) {
    igraph_vector_destroy(sample);
    igraph_Free(sample);
    IGRAPH_FINALLY_CLEAN(2);
  }

  igraph_Free(added);
  igraph_vector_destroy(&vids);
  igraph_vector_destroy(&adjverts);
  igraph_stack_destroy(&stack);
  igraph_vector_destroy(&neis);
  IGRAPH_FINALLY_CLEAN(5);
  return 0;
}

/**
 * \function igraph_motifs_randesu_no
 * \brief Count the total number of motifs in a graph
 *
 * </para><para> 
 * This function counts the total number of motifs in a graph without
 * assigning isomorphism classes to them. 
 *
 * </para><para> 
 * Directed motifs will be counted in directed graphs and undirected
 * motifs in undirected graphs.
 *
 * \param graph The graph object to study.
 * \param no Pointer to an integer type, the result will be stored
 *        here. 
 * \param size The size of the motifs to count.
 * \param cut_prob Vector giving the probabilities that a branch of
 *        the search tree will be cut at a given level.
 * \return Error code.
 * \sa \ref igraph_motifs_randesu(), \ref
 *     igraph_motifs_randesu_estimate(). 
 * 
 * Time complexity: TODO.
 */

int igraph_motifs_randesu_no(const igraph_t *graph, igraph_integer_t *no,
			     int size, const igraph_vector_t *cut_prob) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_t neis;
    
  igraph_vector_t vids;		/* this is G */
  igraph_vector_t adjverts;	/* this is V_E */
  igraph_stack_t stack;		/* this is S */
  long int *added;
  long int father;
  long int i;

  added=igraph_Calloc(no_of_nodes, long int);
  if (added==0) {
    IGRAPH_ERROR("Cannot find motifs", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, added);

  IGRAPH_VECTOR_INIT_FINALLY(&vids, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&adjverts, 0);
  IGRAPH_CHECK(igraph_stack_init(&stack, 0));
  IGRAPH_FINALLY(igraph_stack_destroy, &stack);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);

  *no=0;

  RNG_BEGIN();

  for (father=0; father<no_of_nodes; father++) {
    long int level, s;

    IGRAPH_ALLOW_INTERRUPTION();

    /* init G */
    igraph_vector_clear(&vids); level=0;
    IGRAPH_CHECK(igraph_vector_push_back(&vids, father));
    added[father] += 1; level += 1;
    
    /* init V_E */
    igraph_vector_clear(&adjverts);
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, father, IGRAPH_ALL));
    s=igraph_vector_size(&neis);
    for (i=0; i<s; i++) {
      long int nei=VECTOR(neis)[i];
      if (!added[nei]) {
	IGRAPH_CHECK(igraph_vector_push_back(&adjverts, nei));
	IGRAPH_CHECK(igraph_vector_push_back(&adjverts, father));
      }
      added[nei] += 1;
    }
    
    /* init S */
    igraph_stack_clear(&stack);

    while (level > 1 || !igraph_vector_empty(&adjverts)) {
      igraph_real_t cp=VECTOR(*cut_prob)[level];
      
      if (level==size-1) {
	s=igraph_vector_size(&adjverts)/2;
	for (i=0; i<s; i++) {
	  if (cp!=0 && RNG_UNIF01() < cp) { continue; }
	  (*no) += 1;
	}
      }

      if (level < size-1 && 
	  !igraph_vector_empty(&adjverts) && 
	  (cp==0 || RNG_UNIF01() > cp)) {
	/* yes, step down */
	long int neifather=igraph_vector_pop_back(&adjverts);
	long int nei=igraph_vector_pop_back(&adjverts);
	IGRAPH_CHECK(igraph_vector_push_back(&vids, nei));
	added[nei] += 1; level += 1;

	IGRAPH_CHECK(igraph_stack_push(&stack, neifather));
	IGRAPH_CHECK(igraph_stack_push(&stack, nei));
	IGRAPH_CHECK(igraph_stack_push(&stack, level));

	IGRAPH_CHECK(igraph_neighbors(graph, &neis, nei, IGRAPH_ALL));
	s=igraph_vector_size(&neis);
	for (i=0; i<s; i++) {
	  long int nei2=VECTOR(neis)[i];
	  if (!added[nei2]) {
	    IGRAPH_CHECK(igraph_vector_push_back(&adjverts, nei2));
	    IGRAPH_CHECK(igraph_vector_push_back(&adjverts, nei));
	  }
	  added[nei2] += 1;
	}
      } else {
	/* no, step back */
	long int nei, neifather;
	while (!igraph_stack_empty(&stack) &&
	       level==igraph_stack_top(&stack)-1) {
	  igraph_stack_pop(&stack);
	  nei=igraph_stack_pop(&stack);
	  neifather=igraph_stack_pop(&stack);
	  igraph_vector_push_back(&adjverts, nei);
	  igraph_vector_push_back(&adjverts, neifather);
	}

	nei=igraph_vector_pop_back(&vids);
	added[nei] -= 1; level -= 1;
	IGRAPH_CHECK(igraph_neighbors(graph, &neis, nei, IGRAPH_ALL));
	s=igraph_vector_size(&neis);
	for (i=0; i<s; i++) {
	  added[ (long int) VECTOR(neis)[i] ] -= 1;
	}
	while (!igraph_vector_empty(&adjverts) && 
	       igraph_vector_tail(&adjverts)==nei) {
	  igraph_vector_pop_back(&adjverts);
	  igraph_vector_pop_back(&adjverts);
	}
      }
      
    } /* while */

    /* clear the added vector */
    added[father] -= 1;
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, father, IGRAPH_ALL));
    s=igraph_vector_size(&neis);
    for (i=0; i<s; i++) {
      added[ (long int) VECTOR(neis)[i] ] -= 1;
    }

  } /* for father */

  RNG_END();

  *no /= size;

  igraph_Free(added);
  igraph_vector_destroy(&vids);
  igraph_vector_destroy(&adjverts);
  igraph_stack_destroy(&stack);
  igraph_vector_destroy(&neis);
  IGRAPH_FINALLY_CLEAN(5);
  return 0;
}

/**
 * \function igraph_dyad_census
 * \brief Calculating the dyad census as defined by Holland and Leinhardt
 * 
 * </para><para>
 * Dyad census means classifying each pair of vertices of a directed
 * graph into three categories: mutual, there is an edge from \c a to
 * \c b and also from \c b to \c a; asymmetric, there is an edge
 * either from \c a to \c b or from \c b to \c a but not the other way
 * and null, no edges between \c a and \c b.
 * 
 * </para><para>
 * Holland, P.W. and Leinhardt, S.  (1970).  A Method for Detecting
 * Structure in Sociometric Data.  American Journal of Sociology,
 * 70, 492-513. 
 * \param graph The input graph, a warning is given if undirected as 
 *    the results are undefined for undirected graphs.
 * \param mut Pointer to an integer, the number of mutual dyads is
 *    stored here.
 * \param asym Pointer to an integer, the number of asymmetric dyads
 *    is stored here.
 * \param null Pointer to an integer, the number of null dyads is
 *    stored here.
 * \return Error code.
 *
 * \sa \ref igraph_reciprocity(), \ref igraph_triad_census().
 * 
 * Time complexity: O(|V|+|E|), the number of vertices plus the number
 * of edges.
 */

int igraph_dyad_census(const igraph_t *graph, igraph_integer_t *mut,
		       igraph_integer_t *asym, igraph_integer_t *null) {
  
  igraph_integer_t vc=igraph_vcount(graph);
  igraph_integer_t ec=igraph_ecount(graph);
  igraph_real_t rec;

  if (!igraph_is_directed(graph)) {
    IGRAPH_WARNING("Dyad census called on undirected graph");
  }

  IGRAPH_CHECK(igraph_reciprocity(graph, &rec, 1));
  *mut = rec*ec/(rec+1);
  *asym = ec-2*(*mut);
  *null = vc*(vc-1)/2-(*mut)-(*asym);
  
  return 0;
}

/**
 * \function igraph_triad_census_24
 * TODO
 */

int igraph_triad_census_24(const igraph_t *graph, igraph_integer_t *res2,
			   igraph_integer_t *res4) {
  
  long int vc=igraph_vcount(graph);
  igraph_vector_long_t seen;
  igraph_vector_t *neis, *neis2;
  long int i, j, k, s, neilen, neilen2, ign;
  igraph_adjlist_t adjlist;
  
  IGRAPH_CHECK(igraph_vector_long_init(&seen, vc));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &seen);
  IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
  *res2=*res4=0;

  for (i=0; i<vc; i++) {
    IGRAPH_ALLOW_INTERRUPTION();
    
    neis=igraph_adjlist_get(&adjlist, i);
    neilen=igraph_vector_size(neis);
    /* mark neighbors of i & i itself */
    VECTOR(seen)[i]=i+1;
    ign=0;
    for (j=0; j<neilen; j++) {
      long int nei=VECTOR(*neis)[j];
      if (VECTOR(seen)[nei]==i+1 || VECTOR(seen)[nei]==-(i+1)) {
	/* multiple edges or loop edge */
	VECTOR(seen)[nei]=-(i+1);
	ign++;
      } else {
	VECTOR(seen)[nei]=i+1;
      }
    }
    
    for (j=0; j<neilen; j++) {
      long int nei=VECTOR(*neis)[j];
      if (nei<=i || (j>0 && nei==VECTOR(*neis)[j-1])) { continue; }
      neis2=igraph_adjlist_get(&adjlist, nei);
      neilen2=igraph_vector_size(neis2);
      s=0;
      for (k=0; k<neilen2; k++) {
	long int nei2=VECTOR(*neis2)[k];
	if (k>0 && nei2==VECTOR(*neis2)[k-1]) { continue; }	
	if (VECTOR(seen)[nei2] != i+1 && VECTOR(seen)[nei2] != -(i+1)) {
	  s++;
	}
      }
      if (VECTOR(seen)[nei] > 0) {
	*res2 += vc-s-neilen+ign-1;
      } else {
	*res4 += vc-s-neilen+ign-1;
      }
    }
  }

  igraph_adjlist_destroy(&adjlist);
  igraph_vector_long_destroy(&seen);
  IGRAPH_FINALLY_CLEAN(2);

  return 0;
}

/**
 * \function igraph_triad_census
 * \brief Triad census, as defined by Davis and Leinhardt
 * 
 * </para><para>
 * Calculating the triad census means classifying every triple of
 * vertices in a directed graph. A triple can be in one of 16 states:
 * \clist
 * \cli 003 
 *      A, B, C, the empty graph.
 * \cli 012 
 *      A->B, C, a graph with a single directed edge.
 * \cli 102 
 *      A&lt;->B, C, a graph with a mutual connection between two vertices.
 * \cli 021D 
 *      A&lt;-B->C, the binary out-tree.
 * \cli 021U 
 *      A->B&lt;-C, the binary in-tree.
 * \cli 021C 
 *      A->B->C, the directed line.
 * \cli 111D 
 *      A&lt;->B&lt;-C.
 * \cli 111U 
 *      A&lt;->B->C.
 * \cli 030T 
 *      A->B&lt;-C, A->C.
 * \cli 030C 
 *      A&lt;-B&lt;-C, A->C.
 * \cli 201 
 *      A&lt;->B&lt;->C.
 * \cli 120D 
 *      A&lt;-B->C, A&lt;->C.
 * \cli 120U 
 *      A->B&lt;-C, A&lt;->C.
 * \cli 120C 
 *      A->B->C, A&lt;->C.
 * \cli 210 
 *      A->B&lt;->C, A&lt;->C.
 * \cli 300 
 *      A&lt;->B&lt;->C, A&lt;->C, the complete graph.
 * \endclist
 *
 * </para><para>
 * See also Davis, J.A. and Leinhardt, S.  (1972).  The Structure of
 * Positive Interpersonal Relations in Small Groups.  In J. Berger
 * (Ed.), Sociological Theories in Progress, Volume 2, 218-251. 
 * Boston: Houghton Mifflin.
 * 
 * </para><para>
 * This function calls \ref igraph_motifs_randesu() which is an
 * implementation of the FANMOD motif finder tool, see \ref
 * igraph_motifs_randesu() for details. Note that the order of the
 * triads is not the same for \ref igraph_triad_census() and \ref
 * igraph_motifs_randesu().
 * 
 * \param graph The input graph. A warning is given for undirected
 *   graphs, as the result is undefined for those.
 * \param res Pointer to an initialized vector, the result is stored
 *   here in the same order as given in the list above. Note that this
 *   order is different than the one used by \ref igraph_motifs_randesu().
 * \return Error code.
 * 
 * \sa \ref igraph_motifs_randesu(), \ref igraph_dyad_census().
 * 
 * Time complexity: TODO.
 */

int igraph_triad_census(const igraph_t *graph, igraph_vector_t *res) {

  igraph_vector_t cut_prob;
  igraph_integer_t m2, m4;
  igraph_vector_t tmp;
  igraph_integer_t vc=igraph_vcount(graph);

  if (!igraph_is_directed(graph)) {
    IGRAPH_WARNING("Triad census called on an undirected graph");
  }

  IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&cut_prob, 3); /* all zeros */
  IGRAPH_CHECK(igraph_vector_resize(res, 16));
  IGRAPH_CHECK(igraph_motifs_randesu(graph, &tmp, 3, &cut_prob));
  
  IGRAPH_CHECK(igraph_triad_census_24(graph, &m2, &m4));
  VECTOR(tmp)[1]=m2;
  VECTOR(tmp)[3]=m4;
  VECTOR(tmp)[0]=vc*(vc-1)*(vc-2)/6 - igraph_vector_sum(&tmp);
  
  /* Reorder */
  VECTOR(*res)[0] = VECTOR(tmp)[0];
  VECTOR(*res)[1] = VECTOR(tmp)[1];
  VECTOR(*res)[2] = VECTOR(tmp)[3];
  VECTOR(*res)[3] = VECTOR(tmp)[6];
  VECTOR(*res)[4] = VECTOR(tmp)[2];
  VECTOR(*res)[5] = VECTOR(tmp)[4];
  VECTOR(*res)[6] = VECTOR(tmp)[5];
  VECTOR(*res)[7] = VECTOR(tmp)[9];
  VECTOR(*res)[8] = VECTOR(tmp)[7];
  VECTOR(*res)[9] = VECTOR(tmp)[11];
  VECTOR(*res)[10] = VECTOR(tmp)[10];
  VECTOR(*res)[11] = VECTOR(tmp)[8];
  VECTOR(*res)[12] = VECTOR(tmp)[13];
  VECTOR(*res)[13] = VECTOR(tmp)[12];
  VECTOR(*res)[14] = VECTOR(tmp)[14];
  VECTOR(*res)[15] = VECTOR(tmp)[15];
  
  igraph_vector_destroy(&cut_prob);
  igraph_vector_destroy(&tmp);
  IGRAPH_FINALLY_CLEAN(2);  
  
  return 0;
}

