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

/* The original version of this file was written by Joerg Reichardt 
   The original copyright notice follows here */

/***************************************************************************
                          main.cpp  -  description
                             -------------------
    begin                : Tue Jul 13 11:26:47 CEST 2004
    copyright            : (C) 2004 by 
    email                : 
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "NetDataTypes.h"
#include "NetRoutines.h"
#include "pottsmodel_2.h"

#include "igraph.h"
#include "error.h"
#include "random.h"

/**
 * \function igraph_community_spinglass
 * \brief Community detection based on statistical mechanics
 * 
 * This function implements the community structure detection
 * algorithm proposed by Joerg Reichardt and Stefan Bornholdt. 
 * The algorithm is described in their paper: Statistical Mechanics of 
 * Community Detection, http://arxiv.org/abs/cond-mat/0603718.
 * \param graph The input graph, it may be directed but the direction
 *     of the edge is not used in the algorithm.
 * \param weights The vector giving the edge weights, it may be \c NULL, 
 *     in which case all edges are weighted equally.
 * \param modularity Pointer to a real number, if not \c NULL then the
 *     modularity score of the solution will be stored here, see
 *     M. E. J. Newman and M. Girvan, Phys. Rev. E 69, 026113 (2004)
 *     for details. 
 * \param temperature Pointer to a real number, if not \c NULL then
 *     the temperature at the end of the algorithm will be stored
 *     here.
 * \param membership Pointer to an initialized vector or \c NULL. If
 *     not \c NULL then the result of the clustering will be stored
 *     here, for each vertex the number of its cluster is given, the 
 *     first cluster is numbered zero. The vector will be resized as
 *     needed. 
 * \param csize Pointer to an initialized vector or \c NULL. If not \c
 *     NULL then the sizes of the clusters will stored here in cluster
 *     number order. The vector will be resized as needed.
 * \param spins Integer giving the number of spins, ie. the maximum
 *     number of clusters. Usually it is not a program to give a high
 *     number here, the default was 25 in the original code. Even if
 *     the number of spins is high the number of clusters in the
 *     result might small.
 * \param parupdate A logical constant, whether to update all spins in
 *     parallel. The default for this argument was \c FALSE (ie. 0) in
 *     the original code.
 * \param starttemp Real number, the temperature at the start. The
 *     value of this argument was 1.0 in the original code.
 * \param stoptemp Real number, the algorithm stops at this
 *     temperature. The default was 0.01 in the original code.
 * \param coolfact Real number, the coolinf factor for the simulated
 *     annealing. The default was 0.99 in the original code.
 * \param update_rule The type of the update rule. Possible values: \c
 *     IGRAPH_SPINCOMM_UPDATE_SIMPLE and \c
 *     IGRAPH_SPINCOMM_UPDATE_CONFIG. Basically this parameter defined
 *     the null model based on which the actual clustering is done. If
 *     this is \c IGRAPH_SPINCOMM_UPDATE_SIMPLE then the random graph
 *     (ie. G(n,p)), if it is \c IGRAPH_SPINCOMM_UPDATE then the
 *     configuration model is used. The configuration means that the
 *     baseline for the clustering is a random graph with the same
 *     degree distribution as the input graph.
 * \param gamma Real number. The gamma parameter of the
 *     algorithm. This defined the weight of the missing and existing
 *     links in the quality function for the clustering. The default
 *     value in the original code was 1.0, which is equal weight to 
 *     missing and existing edges. Smaller values make the existing
 *     links contibute more to the energy function which is minimized
 *     in the algorithm. Bigger values make the missing links more
 *     important. (If my understanding is correct.)
 * \return Error code.
 * 
 * \sa igraph_community_spinglass_single() for calculating the community
 * of a single vertex.
 * 
 * Time complexity: TODO.
 */

int igraph_community_spinglass(const igraph_t *graph,
			       const igraph_vector_t *weights,
			       igraph_real_t *modularity,
			       igraph_real_t *temperature,
			       igraph_vector_t *membership, 
			       igraph_vector_t *csize, 
			       igraph_integer_t spins,
			       igraph_bool_t parupdate,
			       igraph_real_t starttemp,
			       igraph_real_t stoptemp,
			       igraph_real_t coolfact,
			       igraph_spincomm_update_t update_rule,
			       igraph_real_t gamma) {

  unsigned long changes, runs;
  igraph_bool_t use_weights=0;
  bool zeroT;
  double Q, kT, acc, prob;
  ClusterList<NNode*> *cl_cur;
  network *net;
  PottsModel *pm;

  /* Check arguments */

  if (spins < 2 || spins > 500) {
    IGRAPH_ERROR("Invalid number of spins", IGRAPH_EINVAL);
  }
  if (update_rule != IGRAPH_SPINCOMM_UPDATE_SIMPLE &&
      update_rule != IGRAPH_SPINCOMM_UPDATE_CONFIG) {
    IGRAPH_ERROR("Invalid update rule", IGRAPH_EINVAL);
  }
  if (weights) {
    if (igraph_vector_size(weights) != igraph_ecount(graph)) {
      IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
    }
    use_weights=1;
  }
  if (coolfact < 0 || coolfact>=1.0) {
    IGRAPH_ERROR("Invalid cooling factor", IGRAPH_EINVAL);
  }
  if (gamma < 0.0) {
    IGRAPH_ERROR("Invalid gamme value", IGRAPH_EINVAL);
  }
  if (starttemp/stoptemp<1.0) {
    IGRAPH_ERROR("starttemp should be larger in absolute value than stoptemp",
		 IGRAPH_EINVAL);
  }
  
  /* Check whether we have a single component */
  igraph_bool_t conn;
  IGRAPH_CHECK(igraph_is_connected(graph, &conn, IGRAPH_WEAK));
  if (!conn) {
    IGRAPH_ERROR("Cannot work with unconnected graph", IGRAPH_EINVAL);
  }

  net = new network;
  net->node_list   =new DL_Indexed_List<NNode*>();
  net->link_list   =new DL_Indexed_List<NLink*>();
  net->cluster_list=new DL_Indexed_List<ClusterList<NNode*>*>();

  /* Transform the igraph_t */
  IGRAPH_CHECK(igraph_i_read_network(graph, weights,
				     net, 0.0, use_weights, 0));

  prob=2.0*net->sum_weights/double(net->node_list->Size())
    /double(net->node_list->Size()-1);

  pm=new PottsModel(net,(unsigned int)spins,update_rule);

  /* initialize the random number generator */
  RNG_BEGIN();
  
  if ((stoptemp==0.0) && (starttemp==0.0)) zeroT=true; else zeroT=false;
  if (!zeroT) kT=pm->FindStartTemp(gamma, prob, starttemp); else kT=stoptemp;
  /* assign random initial configuration */
  pm->assign_initial_conf(-1);
  Q=pm->initialize_Qmatrix();
  runs=0;
  changes=1;

  while (changes>0 && (kT/stoptemp>1.0 || (zeroT && runs<150))) {

    IGRAPH_ALLOW_INTERRUPTION(); /* This is not clean.... */
    
    runs++;
    if (!zeroT) {
      kT*=coolfact;
      if (parupdate) { 
	changes=pm->HeatBathParallelLookup(gamma, prob, kT, 50);
      } else {
	acc=pm->HeatBathLookup(gamma, prob, kT, 50);
	if (acc<(1.0-1.0/double(spins))*0.01) {
	  changes=0; 
	} else { 
	  changes=1;
	}
      }
    } else {
      if (parupdate) { 
	changes=pm->HeatBathParallelLookupZeroTemp(gamma, prob, 50);
      } else {
	acc=pm->HeatBathLookupZeroTemp(gamma, prob, 50);
	/* less than 1 percent acceptance ratio */
	if (acc<(1.0-1.0/double(spins))*0.01) {
	  changes=0; 
	} else { 
	  changes=1;
	}
      }
    }
  } /* while loop */

  pm->WriteClusters(modularity, temperature, csize, membership, kT);

  while (net->link_list->Size()) delete net->link_list->Pop();
  while (net->node_list->Size()) delete net->node_list->Pop();
  while (net->cluster_list->Size())
    {
      cl_cur=net->cluster_list->Pop();
      while (cl_cur->Size()) cl_cur->Pop();
      delete cl_cur;
    }
  delete net->link_list;
  delete net->node_list;
  delete net->cluster_list;
  
  RNG_END();
  
  delete net;
  delete pm;

  return 0;
}

int igraph_spinglass_community(const igraph_t *graph,
			       const igraph_vector_t *weights,
			       igraph_real_t *modularity,
			       igraph_real_t *temperature,
			       igraph_vector_t *membership, 
			       igraph_vector_t *csize, 
			       igraph_integer_t spins,
			       igraph_bool_t parupdate,
			       igraph_real_t starttemp,
			       igraph_real_t stoptemp,
			       igraph_real_t coolfact,
			       igraph_spincomm_update_t update_rule,
			       igraph_real_t gamma) {

  IGRAPH_WARNING("This function was renamed to igraph_community_spinglass");
  return igraph_community_spinglass(graph, weights, modularity, temperature,
				    membership, csize, spins, parupdate,
				    starttemp, stoptemp, coolfact, update_rule,
				    gamma);
}

/**
 * \function igraph_community_spinglass_single
 * \brief Community of a single node based on statistical mechanics
 * 
 * This function implements the community structure detection
 * algorithm proposed by Joerg Reichardt and Stefan Bornholdt. It is
 * described in their paper: Statistical Mechanics of 
 * Community Detection, http://arxiv.org/abs/cond-mat/0603718.
 * 
 * </para><para>
 * This function calculates the community of a single vertex without
 * calculating all the communities in the graph.
 * 
 * \param graph The input graph, it may be directed but the direction
 *    of the edges is not used in the algorithm.
 * \param weights Pointer to a vector with the weights of the edges.
 *    Alternatively \c NULL can be supplied to have the same weight
 *    for every edge.
 * \param vertex The vertex id of the vertex of which ths community is 
 *    calculated.
 * \param community Pointer to an initialized vector, the result, the
 *    ids of the vertices in the community of the input vertex will be
 *    stored here. The vector will be resized as needed.
 * \param cohesion Pointer to a real variable, if not \c NULL the
 *     cohesion index of the community will be stored here.
 * \param adhesion Pointer to a real variable, if not \c NULL the
 *     adhesion index of the community will be stored here.
 * \param inner_links Pointer to an integer, if not \c NULL the 
 *     number of edges within the community is stored here.
 * \param outer_links Pointer to an integer, if not \c NULL the 
 *     number of edges between the community and the rest of the graph  
 *     will be stored here.
 * \param spins The number of spins to use, this can be higher than
 *    the actual number of clusters in the network, in which case some
 *    clusters will contain zero vertices.
 * \param update_rule The type of the update rule. Possible values: \c
 *     IGRAPH_SPINCOMM_UPDATE_SIMPLE and \c
 *     IGRAPH_SPINCOMM_UPDATE_CONFIG. Basically this parameter defined
 *     the null model based on which the actual clustering is done. If
 *     this is \c IGRAPH_SPINCOMM_UPDATE_SIMPLE then the random graph
 *     (ie. G(n,p)), if it is \c IGRAPH_SPINCOMM_UPDATE then the
 *     configuration model is used. The configuration means that the
 *     baseline for the clustering is a random graph with the same
 *     degree distribution as the input graph.
 * \param gamma Real number. The gamma parameter of the
 *     algorithm. This defined the weight of the missing and existing
 *     links in the quality function for the clustering. The default
 *     value in the original code was 1.0, which is equal weight to 
 *     missing and existing edges. Smaller values make the existing
 *     links contibute more to the energy function which is minimized
 *     in the algorithm. Bigger values make the missing links more
 *     important. (If my understanding is correct.)
 * \return Error code.
 * 
 * \sa igraph_community_spinglass() for the traditional version of the 
 * algorithm.
 * 
 * Time complexity: TODO.
 */

int igraph_community_spinglass_single(const igraph_t *graph,
				      const igraph_vector_t *weights,
				      igraph_integer_t vertex,
				      igraph_vector_t *community,
				      igraph_real_t *cohesion,
				      igraph_real_t *adhesion,
				      igraph_integer_t *inner_links,
				      igraph_integer_t *outer_links,
				      igraph_integer_t spins,
				      igraph_spincomm_update_t update_rule,
				      igraph_real_t gamma) {

  igraph_bool_t use_weights=0;
  double prob;
  ClusterList<NNode*> *cl_cur;
  network *net;
  PottsModel *pm;
  char startnode[255];

  /* Check arguments */

  if (spins < 2 || spins > 500) {
    IGRAPH_ERROR("Invalid number of spins", IGRAPH_EINVAL);
  }
  if (update_rule != IGRAPH_SPINCOMM_UPDATE_SIMPLE &&
      update_rule != IGRAPH_SPINCOMM_UPDATE_CONFIG) {
    IGRAPH_ERROR("Invalid update rule", IGRAPH_EINVAL);
  }
  if (weights) {
    if (igraph_vector_size(weights) != igraph_ecount(graph)) {
      IGRAPH_ERROR("Invalid weight vector length", IGRAPH_EINVAL);
    }
    use_weights=1;
  }
  if (gamma < 0.0) {
    IGRAPH_ERROR("Invalid gamme value", IGRAPH_EINVAL);
  }
  if (vertex < 0 || vertex > igraph_vcount(graph)) {
    IGRAPH_ERROR("Invalid vertex id", IGRAPH_EINVAL);
  }
  
  /* Check whether we have a single component */
  igraph_bool_t conn;
  IGRAPH_CHECK(igraph_is_connected(graph, &conn, IGRAPH_WEAK));
  if (!conn) {
    IGRAPH_ERROR("Cannot work with unconnected graph", IGRAPH_EINVAL);
  }

  net = new network;
  net->node_list   =new DL_Indexed_List<NNode*>();
  net->link_list   =new DL_Indexed_List<NLink*>();
  net->cluster_list=new DL_Indexed_List<ClusterList<NNode*>*>();

  /* Transform the igraph_t */
  IGRAPH_CHECK(igraph_i_read_network(graph, weights,
				     net, 0.0, use_weights, 0));

  prob=2.0*net->sum_weights/double(net->node_list->Size())
    /double(net->node_list->Size()-1);

  pm=new PottsModel(net,(unsigned int)spins,update_rule);

  /* initialize the random number generator */
  RNG_BEGIN();

  /* to be exected, if we want to find the community around a particular node*/
  /* the initial conf is needed, because otherwise, 
     the degree of the nodes is not in the weight property, stupid!!! */
  pm->assign_initial_conf(-1);
  snprintf(startnode, 255, "%li", (long int)vertex+1);
  pm->FindCommunityFromStart(gamma, prob, startnode, community,
			     cohesion, adhesion, inner_links, outer_links);
  
  while (net->link_list->Size()) delete net->link_list->Pop();
  while (net->node_list->Size()) delete net->node_list->Pop();
  while (net->cluster_list->Size())
    {
      cl_cur=net->cluster_list->Pop();
      while (cl_cur->Size()) cl_cur->Pop();
      delete cl_cur;
    }
  delete net->link_list;
  delete net->node_list;
  delete net->cluster_list;
  
  RNG_END();

  delete net;
  delete pm;

  return 0;
}

int igraph_spinglass_my_community(const igraph_t *graph,
				  const igraph_vector_t *weights,
				  igraph_integer_t vertex,
				  igraph_vector_t *community,
				  igraph_real_t *cohesion,
				  igraph_real_t *adhesion,
				  igraph_integer_t *inner_links,
				  igraph_integer_t *outer_links,
				  igraph_integer_t spins,
				  igraph_spincomm_update_t update_rule,
				  igraph_real_t gamma) {

  IGRAPH_WARNING("this function was renamed to `igraph_community_spinglass_single");
  return igraph_community_spinglass_single(graph, weights, vertex, community,
					   cohesion, adhesion, inner_links, 
					   outer_links, spins, update_rule,
					   gamma);
}
