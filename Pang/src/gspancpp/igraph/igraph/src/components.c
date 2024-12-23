/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2003, 2004, 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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
#include <string.h>
#include <limits.h>

int igraph_clusters_weak(const igraph_t *graph, igraph_vector_t *membership,
			 igraph_vector_t *csize, igraph_integer_t *no);

int igraph_clusters_strong(const igraph_t *graph, igraph_vector_t *membership,
			   igraph_vector_t *csize, igraph_integer_t *no);

/**
 * \ingroup structural
 * \function igraph_clusters
 * \brief Calculates the (weakly or strongly) connected components in a graph. 
 *
 * \param graph The graph object to analyze.
 * \param membership First half of the result will be stored here. For
 *        every vertex the id of its component is given. The vector
 *        has to be preinitialized and will be resized. Alternatively 
 *        this argument can be \c NULL, in which case it is ignored.
 * \param csize The second half of the result. For every component it
 *        gives its size, the order is defined by the component ids.
 *        The vector has to be preinitialized and will be resized.
 *        Alternatively this argument can be \c NULL, in which
 *        case it is ignored.
 * \param no Pointer to an integer, if not \c NULL then the number of 
 *        clusters will be stored here.
 * \param mode For directed graph this specifies whether to calculate
 *        weakly or strongly connected components. Possible values: 
 *        \c IGRAPH_WEAK,
 *        \c IGRAPH_STRONG. This argument is 
 *        ignored for undirected graphs.
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid mode argument.
 * 
 * Time complexity: O(|V|+|E|),
 * |V| and 
 * |E| are the number of vertices and
 * edges in the graph. 
 */

int igraph_clusters(const igraph_t *graph, igraph_vector_t *membership, 
		    igraph_vector_t *csize, igraph_integer_t *no,
		    igraph_connectedness_t mode) {
  if (mode==IGRAPH_WEAK || !igraph_is_directed(graph)) {
    return igraph_clusters_weak(graph, membership, csize, no);
  } else if (mode==IGRAPH_STRONG) {
    return igraph_clusters_strong(graph, membership, csize, no);
  } else {
    IGRAPH_ERROR("Cannot calculate clusters", IGRAPH_EINVAL);
  }
  
  return 1;
}

int igraph_clusters_weak(const igraph_t *graph, igraph_vector_t *membership,
			 igraph_vector_t *csize, igraph_integer_t *no) {

  long int no_of_nodes=igraph_vcount(graph);
  char *already_added;
  long int first_node, act_cluster_size=0, no_of_clusters=1;
  
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  
  long int i;
  igraph_vector_t neis=IGRAPH_VECTOR_NULL;

  already_added=igraph_Calloc(no_of_nodes,char);
  if (already_added==0) {
    IGRAPH_ERROR("Cannot calculate clusters", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, already_added);

  IGRAPH_DQUEUE_INIT_FINALLY(&q, no_of_nodes > 100000 ? 10000 : no_of_nodes/10);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);

  /* Memory for result, csize is dynamically allocated */
  if (membership) { 
    IGRAPH_CHECK(igraph_vector_resize(membership, no_of_nodes));
  }
  if (csize) { 
    igraph_vector_clear(csize);
  }

  /* The algorithm */

  for (first_node=0; first_node < no_of_nodes; ++first_node) {
    if (already_added[first_node]==1) continue;
    IGRAPH_ALLOW_INTERRUPTION();

    already_added[first_node]=1;
    act_cluster_size=1;
    if (membership) {
      VECTOR(*membership)[first_node]=no_of_clusters-1;
    }
    IGRAPH_CHECK(igraph_dqueue_push(&q, first_node));
    
    while ( !igraph_dqueue_empty(&q) ) {
      long int act_node=igraph_dqueue_pop(&q);
      IGRAPH_CHECK(igraph_neighbors(graph, &neis, act_node, IGRAPH_ALL));
      for (i=0; i<igraph_vector_size(&neis); i++) {
	long int neighbor=VECTOR(neis)[i];
	if (already_added[neighbor]==1) { continue; }
	IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
	already_added[neighbor]=1;
	act_cluster_size++;
	if (membership) {
	  VECTOR(*membership)[neighbor]=no_of_clusters-1;
	}
      }
    }
    no_of_clusters++;
    if (csize) {
      IGRAPH_CHECK(igraph_vector_push_back(csize, act_cluster_size));
    }
  }
  
  /* Cleaning up */
  
  if (no) { *no = no_of_clusters-1; }
  
  igraph_Free(already_added);
  igraph_dqueue_destroy(&q);
  igraph_vector_destroy(&neis);
  IGRAPH_FINALLY_CLEAN(3);
  
  return 0;
}

int igraph_clusters_strong(const igraph_t *graph, igraph_vector_t *membership,
			   igraph_vector_t *csize, igraph_integer_t *no) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_t next_nei=IGRAPH_VECTOR_NULL;
  
  long int i;
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  
  long int no_of_clusters=1;
  long int act_cluster_size;

  igraph_vector_t out=IGRAPH_VECTOR_NULL;
  igraph_vector_t tmp=IGRAPH_VECTOR_NULL;

  /* The result */

  IGRAPH_VECTOR_INIT_FINALLY(&next_nei, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&out, 0);
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
  IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);

  if (membership) {
    IGRAPH_CHECK(igraph_vector_resize(membership, no_of_nodes));
  }
  IGRAPH_CHECK(igraph_vector_reserve(&out, no_of_nodes));

  igraph_vector_null(&out);
  if (csize) {
    igraph_vector_clear(csize);
  }
  
  for (i=0; i<no_of_nodes; i++) {
    IGRAPH_ALLOW_INTERRUPTION();
    IGRAPH_CHECK(igraph_neighbors(graph, &tmp, i, IGRAPH_OUT));
    if (VECTOR(next_nei)[i] > igraph_vector_size(&tmp)) { continue; }
    
    IGRAPH_CHECK(igraph_dqueue_push(&q, i));
    while (!igraph_dqueue_empty(&q)) {
      long int act_node=igraph_dqueue_back(&q);
      IGRAPH_CHECK(igraph_neighbors(graph, &tmp, act_node, IGRAPH_OUT));
      if (VECTOR(next_nei)[act_node]==0) {
	/* this is the first time we've met this vertex */
	VECTOR(next_nei)[act_node]++;
      } else if (VECTOR(next_nei)[act_node] <= igraph_vector_size(&tmp)) {
	/* we've already met this vertex but it has more children */
	long int neighbor=VECTOR(tmp)[(long int)VECTOR(next_nei)[act_node]-1];
	if (VECTOR(next_nei)[neighbor] == 0) {
	  IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
	}
	VECTOR(next_nei)[act_node]++;
      } else {
	/* we've met this vertex and it has no more children */
	IGRAPH_CHECK(igraph_vector_push_back(&out, act_node));
	igraph_dqueue_pop_back(&q);
      }
    } /* while q */
  }  /* for */

  /* OK, we've the 'out' values for the nodes, let's use them in
     descreasing order with the help of a heap */

  igraph_vector_null(&next_nei);                            /* mark already
							added vertices */
  while (!igraph_vector_empty(&out)) {
    long int grandfather=igraph_vector_pop_back(&out);
    IGRAPH_ALLOW_INTERRUPTION();
    if (VECTOR(next_nei)[grandfather] != 0) { continue; }
    VECTOR(next_nei)[grandfather]=1;
    act_cluster_size=1;
    if (membership) {
      VECTOR(*membership)[grandfather]=no_of_clusters-1;
    }
    IGRAPH_CHECK(igraph_dqueue_push(&q, grandfather));
    
    while (!igraph_dqueue_empty(&q)) {
      long int act_node=igraph_dqueue_pop_back(&q);
      IGRAPH_CHECK(igraph_neighbors(graph, &tmp, act_node, IGRAPH_IN));
      for (i=0; i<igraph_vector_size(&tmp); i++) {
	long int neighbor=VECTOR(tmp)[i];
	if (VECTOR(next_nei)[neighbor] != 0) { continue; }
	IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
	VECTOR(next_nei)[neighbor]=1;
	act_cluster_size++;
	if (membership) {
	  VECTOR(*membership)[neighbor]=no_of_clusters-1;
	}
      }
    }
    no_of_clusters++;
    if (csize) {
      IGRAPH_CHECK(igraph_vector_push_back(csize, act_cluster_size));
    }
  }
  
  if (no) { *no=no_of_clusters-1; }

  /* Clean up, return */

  igraph_vector_destroy(&out);
  igraph_vector_destroy(&tmp);
  igraph_dqueue_destroy(&q);
  igraph_vector_destroy(&next_nei);
  IGRAPH_FINALLY_CLEAN(4);

  return 0;
}

int igraph_is_connected_weak(const igraph_t *graph, igraph_bool_t *res);

/**
 * \ingroup structural
 * \function igraph_is_connected
 * \brief Decides whether the graph is (weakly or strongly) connected.
 * 
 * \param graph The graph object to analyze.
 * \param res Pointer to a logical variable, the result will be stored
 *        here. 
 * \param mode For directed graph this specifies whether to calculate
 *        weak or strong connectedness. Possible values: 
 *        \c IGRAPH_WEAK,
 *        \c IGRAPH_STRONG. This argument is 
 *        igrored for undirected graphs.
 * \return Error code:
 *        \c IGRAPH_EINVAL: invalid mode argument.
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices 
 * plus the number of edges in the graph.
 */


int igraph_is_connected(const igraph_t *graph, igraph_bool_t *res, 
			igraph_connectedness_t mode) {
  if (mode==IGRAPH_WEAK || !igraph_is_directed(graph)) {
    return igraph_is_connected_weak(graph, res);
  } else if (mode==IGRAPH_STRONG) {
    int retval;
    igraph_integer_t no;
    retval = igraph_clusters_strong(graph, 0, 0, &no);
    *res = (no==1);
    return retval;
  } else {
    IGRAPH_ERROR("mode argument", IGRAPH_EINVAL);
  }
  return 0;
}

int igraph_is_connected_weak(const igraph_t *graph, igraph_bool_t *res) {

  long int no_of_nodes=igraph_vcount(graph);
  char *already_added;
  igraph_vector_t neis=IGRAPH_VECTOR_NULL;
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  
  long int i, j;

  already_added=igraph_Calloc(no_of_nodes, char);
  if (already_added==0) {
    IGRAPH_ERROR("is connected (weak) failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, already_added); /* TODO: hack */

  IGRAPH_DQUEUE_INIT_FINALLY(&q, 10);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  
  /* Try to find at least two clusters */
  already_added[0]=1;
  IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
  
  j=1;
  while ( !igraph_dqueue_empty(&q)) {
    long int actnode=igraph_dqueue_pop(&q);
    IGRAPH_ALLOW_INTERRUPTION();
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, actnode, IGRAPH_ALL));
    for (i=0; i <igraph_vector_size(&neis); i++) {
      long int neighbor=VECTOR(neis)[i];
      if (already_added[neighbor] != 0) { continue; }
      IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
      j++;
      already_added[neighbor]++;
    }
  }
  
  /* Connected? */
  *res = (j == no_of_nodes);

  igraph_Free(already_added);
  igraph_dqueue_destroy(&q);
  igraph_vector_destroy(&neis);
  IGRAPH_FINALLY_CLEAN(3);

  return 0;
}

void igraph_i_decompose_free(igraph_vector_ptr_t *complist) {
  long int i;
  for (i=0; i<igraph_vector_ptr_size(complist); i++) {
    if (VECTOR(*complist)[i] != 0) {
      igraph_destroy(VECTOR(*complist)[i]);
      igraph_free(VECTOR(*complist)[i]);
    }
  }
}

/**
 * \function igraph_decompose
 * \brief Decompose a graph into connected components.
 * 
 * Create separate graph for each component of a graph. Note that the
 * vertex ids in the new graphs will be different than in the original 
 * graph. (Except if there is only one component in the original graph.)
 * 
 * \param graph The original graph.
 * \param components This pointer vector will contain pointers to the
 *   subcomponent graphs. It should be initialized before calling this
 *   function and will be resized to hold the graphs. Don't forget to 
 *   call \ref igraph_destroy() and igraph_free() on the elements of
 *   this pointer vector to free unneeded memory.
 * \param mode Either \c IGRAPH_WEAK or \c IGRAPH_STRONG for weakly
 *    and strongly connected components respectively. Right now only
 *    the former is implemented.
 * \param maxcompno The maximum number of components to return. The
 *    first \p maxcompno components will be returned (which hold at
 *    least \p minelements vertices, see the next parameter), the
 *    others will be ignored. Supply -1 here if you don't want to limit
 *    the number of components.
 * \param minelements The minimum number of vertices a component
 *    should contain in order to place it in the \p components
 *    vector. Eg. supply 2 here to ignore isolate vertices.
 * \return Error code, \c IGRAPH_ENOMEM if there is not enough memory
 *   to perform the operation.
 *
 * Added in version 0.2.</para><para>
 * 
 * Time complexity: O(|V|+|E|), the number of vertices plus the number
 * of edges. 
 */

int igraph_decompose(const igraph_t *graph, igraph_vector_ptr_t *components, 
		     igraph_connectedness_t mode,
		     long int maxcompno, long int minelements) {

  long int actstart;
  long int no_of_nodes=igraph_vcount(graph);
  long int resco=0;		/* number of graphs created so far */ 
  char *already_added;
  igraph_dqueue_t q;
  igraph_vector_t verts;
  igraph_vector_t neis;
  long int i;
  igraph_t *newg;

  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_WEAK;
  }

  if (mode != IGRAPH_WEAK) {
    IGRAPH_ERROR("only 'IGRAPH_WEAK' is implemented", IGRAPH_EINVAL);
  }

  if (maxcompno<0) {
    maxcompno=LONG_MAX;
  }

  already_added=igraph_Calloc(no_of_nodes, char);
  if (already_added==0) {
    IGRAPH_ERROR("Cannot decompose graph", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, already_added);

  IGRAPH_CHECK(igraph_dqueue_init(&q, 100));
  IGRAPH_VECTOR_INIT_FINALLY(&verts, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  igraph_vector_ptr_clear(components);
  IGRAPH_FINALLY(igraph_i_decompose_free, components);
  
  for(actstart=0; resco<maxcompno && actstart < no_of_nodes; actstart++) {
    
    if (already_added[actstart]) { continue; }
    IGRAPH_ALLOW_INTERRUPTION();
    
    igraph_vector_clear(&verts);
    already_added[actstart]=1;
    IGRAPH_CHECK(igraph_vector_push_back(&verts, actstart));
    IGRAPH_CHECK(igraph_dqueue_push(&q, actstart));
    
    while (!igraph_dqueue_empty(&q) ) {
      long int actvert=igraph_dqueue_pop(&q);
      IGRAPH_CHECK(igraph_neighbors(graph, &neis, actvert, IGRAPH_ALL));
      for (i=0; i<igraph_vector_size(&neis); i++) {
	long int neighbor=VECTOR(neis)[i];
	if (already_added[neighbor]==1) { continue; }
	IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
	IGRAPH_CHECK(igraph_vector_push_back(&verts, neighbor));
	already_added[neighbor]=1;
      }
    }
    
    /* ok, we have a component */
    if (igraph_vector_size(&verts)<minelements) { continue; }

    newg=igraph_Calloc(1, igraph_t);
    if (newg==0) {
      IGRAPH_ERROR("Cannot decompose graph", IGRAPH_ENOMEM);
    }
    IGRAPH_CHECK(igraph_vector_ptr_push_back(components, newg));
    IGRAPH_FINALLY(igraph_destroy, newg);
    IGRAPH_CHECK(igraph_subgraph(graph, newg, 
				 igraph_vss_vector(&verts)));
    IGRAPH_FINALLY_CLEAN(1);
    resco++;
    
  } /* for actstart++ */

  igraph_dqueue_destroy(&q);
  igraph_free(already_added);
  igraph_vector_destroy(&verts);
  igraph_vector_destroy(&neis);
  
  IGRAPH_FINALLY_CLEAN(4);
  return 0;
}

/** 
 * \function igraph_articulation_points
 * Find the articulation points in a graph.
 * 
 * A vertex is an articulation point if its removal increases 
 * the number of connected components in the graph.
 * \param graph The input graph.
 * \param res Pointer to an initialized vector, the 
 *    articulation points will be stored here.
 * \return Error code.
 * 
 * Time complexity: O(|V|+|E|), linear in the number of vertices and edges.
 * 
 * \sa \ref igraph_biconnected_components(), \ref igraph_clusters()
 */

int igraph_articulation_points(const igraph_t *graph,
			       igraph_vector_t *res) {

  igraph_integer_t no;
  return igraph_biconnected_components(graph, &no, 0, res);
}

/**
 * \function igraph_biconnected_components
 * Calculate biconnected components
 * 
 * A graph is biconnected if the removal of any single vertex (and
 * its adjacent edges) does not disconnect it.
 * 
 * </para><para>
 * A biconnected component of a graph is a maximal biconnected
 * subgraph of it. The biconnected components of a graph can be given
 * by the partition of its edges: every edge is a member of exactly
 * one biconnected component. Note that this is not true for
 * vertices: the same vertex can be part of many biconnected
 * components.
 * \param graph The input graph
 * \param no The number of biconnected components will be stored here.
 * \param components If not a NULL points, then the found components
 *     are stored here, in a list of vectors. Every vector in the list
 *     is a biconnected component, represented by its edges. More precisely, 
 *     a spanning tree of the biconnected component is returned.
 *     Note you'll have to 
 *     destroy each vector first by calling \ref igraph_vector_destroy()
 *     and then <code>free()</code> on it, plus you need to call 
 *     \ref igraph_vector_ptr_destroy() on the list to regain all 
 *     allocated memory.
 * \param articulation_points If not a NULL pointer, then the 
 *     articulation points of the graph are stored in this vector.
 *     A vertex is an articulation point if its removal increases the 
 *     number of (weakly) connected components in the graph.
 * \return Error code.
 * 
 * Time complexity: O(|V|+|E|), linear in the number of vertices and edges.
 * 
 * \sa \ref igraph_articulation_points(), \ref igraph_clusters().
 */

int igraph_biconnected_components(const igraph_t *graph,
				  igraph_integer_t *no,
				  igraph_vector_ptr_t *components,
				  igraph_vector_t *articulation_points) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_long_t nextptr;
  igraph_vector_long_t num, low;
  igraph_vector_bool_t found;
  igraph_vector_t *adjedges;
  igraph_stack_t path;
  igraph_vector_t edgestack;
  igraph_adjedgelist_t adjedgelist;
  long int i, counter, rootdfs=0;  

  IGRAPH_CHECK(igraph_vector_long_init(&nextptr, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &nextptr);
  IGRAPH_CHECK(igraph_vector_long_init(&num, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &num);
  IGRAPH_CHECK(igraph_vector_long_init(&low, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &low);
  IGRAPH_CHECK(igraph_vector_bool_init(&found, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_bool_destroy, &found);

  IGRAPH_CHECK(igraph_stack_init(&path, 100));
  IGRAPH_FINALLY(igraph_stack_destroy, &path);
  IGRAPH_VECTOR_INIT_FINALLY(&edgestack, 0);
  IGRAPH_CHECK(igraph_vector_reserve(&edgestack, 100));

  IGRAPH_CHECK(igraph_adjedgelist_init(graph, &adjedgelist, IGRAPH_ALL));
  IGRAPH_FINALLY(igraph_adjedgelist_destroy, &adjedgelist);

  if (no) {
    *no=0;
  }
  if (components) {
    igraph_vector_ptr_clear(components);
  }
  if (articulation_points) {
    igraph_vector_clear(articulation_points);
  }

  for (i=0; i<no_of_nodes; i++) {
    
    if (VECTOR(low)[i] != 0) { continue; } /* already visited */

    IGRAPH_ALLOW_INTERRUPTION();

    IGRAPH_CHECK(igraph_stack_push(&path, i));
    counter=1; 
    rootdfs=0;
    VECTOR(low)[i]=VECTOR(num)[i]=counter++;
    while (!igraph_stack_empty(&path)) {
      long int n;
      long int act=igraph_stack_top(&path);
      long int actnext=VECTOR(nextptr)[act];
      
      adjedges=igraph_adjedgelist_get(&adjedgelist, act);
      n=igraph_vector_size(adjedges);
      if (actnext < n) {
	/* Step down (maybe) */
	long int edge=VECTOR(*adjedges)[actnext];
	long int nei=IGRAPH_OTHER(graph, edge, act);
	if (VECTOR(low)[nei] == 0) {
	  if (act==i) { rootdfs++; }
	  IGRAPH_CHECK(igraph_vector_push_back(&edgestack, edge));
	  IGRAPH_CHECK(igraph_stack_push(&path, nei));
	  VECTOR(low)[nei] = VECTOR(num)[nei]=counter++;
	} else {
	  /* Update low value if needed */
	  if (VECTOR(num)[nei] < VECTOR(low)[act]) {
	    VECTOR(low)[act]=VECTOR(num)[nei];
	  }
	}
	VECTOR(nextptr)[act] += 1;
      } else {
	/* Step up */
	igraph_stack_pop(&path);
	if (!igraph_stack_empty(&path)) {
	  long int prev=igraph_stack_top(&path);
	  /* Update LOW value if needed */
	  if (VECTOR(low)[act] < VECTOR(low)[prev]) {
	    VECTOR(low)[prev] = VECTOR(low)[act];
	  }
	  /* Check for articulation point */
	  if (VECTOR(low)[act] >= VECTOR(num)[prev]) {
	    if (articulation_points && !VECTOR(found)[prev] 
		&& prev != i /* the root */) {
	      IGRAPH_CHECK(igraph_vector_push_back(articulation_points, prev));
	      VECTOR(found)[prev] = 1;
	    }
	    if (no) { *no += 1; }
	    if (components) {
	      igraph_vector_t *v=igraph_Calloc(1, igraph_vector_t);
	      IGRAPH_CHECK(igraph_vector_init(v, 0));
	      while (!igraph_vector_empty(&edgestack)) {
		long int e=igraph_vector_pop_back(&edgestack);
		IGRAPH_CHECK(igraph_vector_push_back(v, e));
		if (IGRAPH_FROM(graph,e)==prev || IGRAPH_TO(graph,e)==prev) {
		  break;
		}
	      }
	      IGRAPH_CHECK(igraph_vector_ptr_push_back(components, v));
	    }
	  }
	} /* !igraph_stack_empty(&path) */
      }
      
    } /* !igraph_stack_empty(&path) */
    
    if (articulation_points && rootdfs >= 2) {
      IGRAPH_CHECK(igraph_vector_push_back(articulation_points, i));
    }

  } /* i < no_of_nodes */

  igraph_adjedgelist_destroy(&adjedgelist);
  igraph_vector_destroy(&edgestack);
  igraph_stack_destroy(&path);
  igraph_vector_bool_destroy(&found);
  igraph_vector_long_destroy(&low);
  igraph_vector_long_destroy(&num);
  igraph_vector_long_destroy(&nextptr);
  IGRAPH_FINALLY_CLEAN(7);

  return 0;
}

