/*
    gSpanSubs 2011/03/08 implemented by Lei Zhao based on:

    $Id: misc.cpp,v 1.6 2004/05/21 05:50:13 taku-ku Exp $;

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
*/

#include "gspan.h"

#include <assert.h>

namespace GSPAN {

const RMPath &DFSCode::buildRMPath ()
{
    rmpath.clear ();

    int old_from = -1;

    for (int i = size() - 1 ; i >= 0 ; --i) {
        if ((*this)[i].from < (*this)[i].to && // forward
                (rmpath.empty() || old_from == (*this)[i].to))
        {
            rmpath.push_back (i);
            old_from = (*this)[i].from;
        }
    }
    assert(rmpath.size() <= size());
    return rmpath;
}

void History::buildSubsInfo()
// For gSpanSubs: add more information in history for subs
{
    sorted_vertex.clear();
    vertex_rank.clear(); 
    int vertex_counter = 0;

    for (unsigned int i = 0; i < size(); i ++)
    // sort the vertex as the sequence they are discovered
    {
        if (vertex_rank.count((*this)[i] -> from) == 0)
        {
            vertex_rank[(*this)[i] -> from] = vertex_counter;
            sorted_vertex.push_back((*this)[i] -> from);
            vertex_counter ++;
        }
        if (vertex_rank.count((*this)[i] -> to) == 0)
        {
            vertex_rank[(*this)[i] -> to] = vertex_counter;
            sorted_vertex.push_back((*this)[i] -> to);
            vertex_counter ++;
        }
    }

    int last_vertex = sorted_vertex.back();
    for (unsigned int i = size() - 1; i >= 0; i --)
    // start from the end of history, find the last forward edge
        if ((*this)[i] -> to == last_vertex &&
                (*this)[i] -> from != last_vertex)
        {
            lastFwdEdge = (*this)[i];
            break;
        }
}

void History::build (Graph &graph, PDFS *e)
{
    // first build history
    clear ();
    edge.clear ();
    edge.resize (graph.edge_size());
    vertex.clear ();
    vertex.resize (graph.size());

    if (e) {
        push_back (e->edge);
        edge[e->edge->id] = vertex[e->edge->from] = vertex[e->edge->to] = 1;

        for (PDFS *p = e->prev ; p ; p = p->prev) {
            push_back (p->edge);	// this line eats 8% of overall instructions(!)
            edge[p->edge->id] = vertex[p->edge->from] = vertex[p->edge->to] = 1;
        }
        std::reverse (begin(), end());
    }
}

/* get_forward_pure ()
   e1 (from1, elabel1, to1)
   from から繋がる edge e2(from2, elabel2, to2) を返す.

   minlabel <= elabel2,
   (elabel1 < elabel2 ||
   (elabel == elabel2 && tolabel1 < tolabel2) の条件をみたす.
   (elabel1, to1) のほうが先に探索されるべき
   また, いままで見た vertex には逝かない (backward のやくめ)
*/
bool get_forward_rmpath (Graph &graph, Edge *e, int minlabel, History& history, EdgeList &result)
{
    result.clear ();
    assert (e->to >= 0 && e->to < static_cast<int>(graph.size ()));
    assert (e->from >= 0 && e->from < static_cast<int>(graph.size ()));
    int tolabel = graph[e->to].label;

    for (Vertex::edge_iterator it = graph[e->from].edge.begin() ;
            it != graph[e->from].edge.end() ; ++it)
    {
        int tolabel2 = graph[it->to].label;
        if (e->to == it->to || minlabel > tolabel2 || history.hasVertex (it->to))
            continue; // mixture of test and pre-scanning?

        if (e->elabel < it->elabel || (e->elabel == it->elabel && tolabel <= tolabel2))
            result.push_back (&(*it)); // Edge not yet tested => add
    }

    return (! result.empty());
}

/* get_forward_pure ()
   e (from, elabel, to)
   to から繋がる edge を返す
   ただし, minlabel より大きいものにしかいかない (DFSの制約)
   また, いままで見た vertex には逝かない (backward のやくめ)
 */
bool get_forward_pure (Graph &graph, Edge *e, int minlabel, History& history, EdgeList &result)
{
    result.clear ();

    assert (e->to >= 0 && e->to < static_cast<int>(graph.size ()));

    /* Walk all edges leaving from vertex e->to.
     */
    for (Vertex::edge_iterator it = graph[e->to].edge.begin() ;
            it != graph[e->to].edge.end() ; ++it)
    {
        /* -e-> [e->to] -it-> [it->to]
         */
        assert (it->to >= 0 && it->to < static_cast<int>(graph.size ()));
        if (minlabel > graph[it->to].label || history.hasVertex (it->to))
            continue;

        result.push_back (&(*it));
    }

    return (! result.empty());
}

/* graph の vertex からはえる edge を探す
   ただし, fromlabel <= tolabel の性質を満たす.
*/
bool get_forward_root (Graph &g, Vertex &v, EdgeList &result)
{
    result.clear ();
    for (Vertex::edge_iterator it = v.edge.begin(); it != v.edge.end(); ++it) {
        assert (it->to >= 0 && it->to < static_cast<int>(g.size ()));
        if (v.label <= g[it->to].label)
            result.push_back (&(*it));
    }

    return (! result.empty());
}

/* get_backward (graph, e1, e2, history);
   e1 (from1, elabel1, to1)
   e2 (from2, elabel2, to2)
   to2 -> from1 に繋がるかどうかしらべる.

   (elabel1 < elabel2 ||
   (elabel == elabel2 && tolabel1 < tolabel2) の条件をみたす. (elabel1, to1) のほうが先に探索されるべき
 */
Edge *get_backward (Graph &graph, Edge* e1, Edge* e2, History& history)
{
    if (e1 == e2)
        return 0;

    assert (e1->from >= 0 && e1->from < static_cast<int>(graph.size ()));
    assert (e1->to >= 0 && e1->to < static_cast<int>(graph.size ()));
    assert (e2->to >= 0 && e2->to < static_cast<int>(graph.size ()));

    for (Vertex::edge_iterator it = graph[e2->to].edge.begin() ;
            it != graph[e2->to].edge.end() ; ++it)
    {
        if (history.hasEdge (it->id))
            continue;

        if ( (it->to == e1->from) &&
                ( (e1->elabel < it->elabel) ||
                  ( (e1->elabel == it->elabel) &&
                    (graph[e1->to].label <= graph[e2->to].label) )
                ) )
        {
            return &(*it);
        }
    }

    return 0;
}
bool get_forward_root_subs (Graph &g, Vertex &v, EdgeList &result, History history)
// For gSpanSubs: one-node grow from subs
{
    result.clear ();
    for (Vertex::edge_iterator it = v.edge.begin(); it != v.edge.end(); ++it)
        if (!history.hasEdge(it -> id) &&
                !history.hasVertex(it -> to))
            result.push_back (&(*it));
    return (! result.empty());
}

bool get_backward_root_subs (Graph &g, Vertex &v, EdgeList &result, History history)
// For gSpanSubs: one-edge grow within subs
{
    result.clear ();
    for (Vertex::edge_iterator it = v.edge.begin(); it != v.edge.end(); ++it)
        if (!history.hasEdge(it -> id) &&
                history.hasVertex(it -> to) &&
                history.get_vertex_rank(it -> to) <= history.get_vertex_rank(it -> from))
            result.push_back (&(*it));
    return (! result.empty());
}

bool get_backward_subs (Graph &graph, Edge* e, unsigned int subsNodeCount, History& history, EdgeList &result)
// For gSpanSubs: get BFS-based backward edge
{
    result.clear();
    unsigned int efrom = history.get_vertex_rank(e -> from);
    unsigned int eto = history.get_vertex_rank(e -> to);
    if (efrom < eto) // forward edge
    {
        /* if current edge is (3, 8)
         * add (8, 4), (8, 5), (8, 6), (8, 7), (8, 8) to backward edges
         */
        for (Vertex::edge_iterator it = graph[e->to].edge.begin() ;
                it != graph[e->to].edge.end() ; ++it)
            if (!history.hasEdge (it->id) &&
                    history.hasVertex(it->to) &&
                    history.get_vertex_rank(it->to) > efrom &&
                    history.get_vertex_rank(it->to) <= eto)
                result.push_back (&(*it));
    }
    else // backward edge
    {
        /* if current pattern has 10 vertices, and current edge is (8, 5)
         * first add (8, 6), (8, 7), (8, 8) to backward edges
         */
        for (Vertex::edge_iterator it = graph[e->from].edge.begin();
                it != graph[e->from].edge.end(); it ++)
        {
            if (!history.hasEdge(it->id) &&
                    history.hasVertex(it->to) &&
                    history.get_vertex_rank(it->to) > eto &&
                    history.get_vertex_rank(it->to) <= efrom)
                result.push_back(&(*it));
        }
        if (efrom < subsNodeCount) // backward edge is in subs
        {
            /* if subs has 10 vertices, and current edge is (8, 5)
             * then add (9, 0), (9, 1),...(9, 9) to backward edges
             */
            for(unsigned int i = efrom + 1; i < subsNodeCount; i ++)
                for (Vertex::edge_iterator it = graph[history.get_vertex(i)].edge.begin();
                        it != graph[history.get_vertex(i)].edge.end(); it ++)
                    if (!history.hasEdge(it -> id) &&
                            history.hasVertex(it -> to) &&
                            history.get_vertex_rank(it -> to) <= i)
                        result.push_back (&(*it));
        }
    }
    return (! result.empty());
}

bool get_forward_subs (Graph &graph, Edge *e, unsigned int subsNodeCount, History& history, EdgeList &result)
// For gSpanSubs: get BFS-based forward edge
/* Edge e here is the last forward edge in the pattern */
{
    result.clear ();
    unsigned int start_rank;

    /* If the last forward edge of the current pattern is an edge in subs
     * then the current pattern is got by adding some backward edges
     * among vertices in subs, forward edge should be all edges going
     * from a vertex in subs to a vertex outside subs
     */
    if (history.get_vertex_rank(e->to) == subsNodeCount - 1)
        start_rank = 0;
    else
    {
        /* If the current pattern has more vertices than subs:
         * if the current edge is a backward edge,
         * it must from the last discoverted vertex,
         * which is the to vertex of the last forward edge.
         * Say the current edge is (8, 5), last forward edge is (3, 8)
         * we should add (3, 9), (4, 9), ... (8, 9) to forward edges
         * If the current edge is a forward edge, say (5, 8)
         * we should add (5, 9), (6, 9), (7, 9), (8, 9) to forward edges
         */
        for (Vertex::edge_iterator it = graph[e->from].edge.begin() ;
                it != graph[e->from].edge.end() ; ++it)
            /* for all the edges starting from the same vertex as this edge
             * forward edges are the ones that have larger edge label
             * or the ones that have the same edge label, but larger destination vertex label
             */
        {
            if (!history.hasEdge(it -> id) &&
                    !history.hasVertex(it -> to) &&
                    (it -> elabel > e -> elabel ||
                     (it -> elabel == e -> elabel &&
                      graph[it -> to].label >= graph[e -> to].label)))
                result.push_back (&(*it));
        }
        start_rank = history.get_vertex_rank(e -> from) + 1;
    }

    for(unsigned int i = start_rank; i < history.nodeCount(); i ++)
        /* for all the nodes added to DFS_CODE after the from vertex
         * all non-included edges are forward edges
         */
        for(Vertex::edge_iterator it = graph[history.get_vertex(i)].edge.begin();
                it != graph[history.get_vertex(i)].edge.end(); it ++)
        {
            if (!history.hasEdge(it -> id) &&
                    !history.hasVertex(it -> to))
                result.push_back (&(*it));
        }
    return (! result.empty());
}
}

