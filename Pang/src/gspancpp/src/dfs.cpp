/*
    gSpanSubs 2011/03/08 implemented by Lei Zhao based on:

    dfs.cpp, v 1.0 2009/12/21 by Marisa Thoma

    is a modification based on :

    $Id: dfs.cpp,v 1.3 2004/05/21 05:50:13 taku-ku Exp $;

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

    Modifications: added some format converters
*/

#include "gspan.h"

#include <cstring>
#include <assert.h>

#include <string>
#include <iterator>
#include <set>

namespace GSPAN {

/* Build a DFS code from a given graph.
 */
void
DFSCode::fromGraph (Graph &g)
{
    clear ();

    EdgeList edges;
    for (unsigned int from = 0 ; from < g.size () ; ++from) {
        if (get_forward_root (g, g[from], edges) == false)
            continue;

        for (EdgeList::iterator it = edges.begin () ; it != edges.end () ; ++it)
            push (from, (*it)->to, g[(*it)->from].label, (*it)->elabel, g[(*it)->to].label);
    }
}

/* Clear current DFS code and build code from given compressed code.
 */
void DFSCode::fromCompressedDFSCode(std::vector<int> & compressed) {
    clear ();

    for (unsigned int from = 0 ; from < (compressed.size()-1)/3 ; ++from) {
        std::cerr<<"error: sorry, DFSCode::fromCompressedDFSCode not implemented yet"<<std::endl;
        exit(1);
    }
}


bool DFSCode::toGraph (Graph &g) const
{
    g.clear ();

    for (DFSCode::const_iterator it = begin(); it != end(); ++it) {
        if (it -> fromlabel != -1 || it -> tolabel != -1)
            /* For gSpanSubs: after support subs, it is possible that the current from/to are not the last discovered vertex
             */
            g.resize (std::max (it->from, it->to) + 1);

        if (it->fromlabel != -1) {
            assert(g[it->from].label == -1 || g[it->from].label == it->fromlabel);
            g[it->from].label = it->fromlabel;
        }
        if (it->tolabel != -1) {
            assert(g[it->to].label == -1 || g[it->to].label == it->tolabel);
            g[it->to].label = it->tolabel;
        }

        g[it->from].push (it->from, it->to, it->elabel);
        if (g.directed == false)
            g[it->to].push (it->to, it->from, it->elabel);
    }

    g.buildEdge ();

    return (true);
}

unsigned int
DFSCode::nodeCount (void) const
{
    unsigned int nodecount = 0;

    for (DFSCode::const_iterator it = begin() ; it != end() ; ++it)
        nodecount = std::max (nodecount, (unsigned int) (std::max (it->from, it->to) + 1));

    return (nodecount);
}


std::ostream &DFSCode::write (std::ostream &os)
{
    if (size() == 0) return os;
    for (unsigned int i = 0; i < size(); ++i)
        os << '(' << (*this)[i].from << ' ' << (*this)[i].to << ' ' << (*this)[i].fromlabel << ' ' << (*this)[i].elabel << ' ' << (*this)[i].tolabel << ')' << " -> "; // For gSpanSubs
    os << std::endl;
    return os;
}

unsigned int DFSCode::edgeCount (void) const {
// For gSpanSubs
    unsigned int edgeCount = 0;
    for (DFSCode::const_iterator it = begin(); it != end(); ++it)
        edgeCount ++;
    return edgeCount;
}

std::vector<int> DFSCode::getCompressedDFSCode() const {
    std::vector<int> cDFSCode(1+size()*3);
    if (size() == 0) {
        cDFSCode.clear();
        return cDFSCode;
    }
    cDFSCode[0] = at(0).fromlabel;

    for (unsigned int i = 0; i < size(); i++) {
        DFS dfstemp = at(i);
        if (dfstemp.from < dfstemp.to) { // forward edge
            cDFSCode[i*3+2] = dfstemp.from;
            cDFSCode[i*3+3] = dfstemp.tolabel;
        } else { // backward edge
            cDFSCode[i*3+2] = -1;
            cDFSCode[i*3+3] = dfstemp.to;
        }
        cDFSCode[i*3+1] = dfstemp.elabel;
    }
    return cDFSCode;
}

bool DFSCode::is_min ()
// For gSpanSubs: function originally in the scope of gSpan
{
    if (size() == 1)
    {
        if (back().fromlabel <= back().tolabel)
            return (true);
        else
            return (false);
    }

    DFSCode DFSCodeIsMin;
    Graph g;
    toGraph (g);

    Projected_map3 root;
    EdgeList       edges;

    for (unsigned int from = 0; from < g.size() ; ++from)
        if (get_forward_root (g, g[from], edges))
            for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it)
                root[g[from].label][(*it)->elabel][g[(*it)->to].label].push (0, *it, 0);

    Projected_iterator3 fromlabel = root.begin();
    Projected_iterator2 elabel    = fromlabel->second.begin();
    Projected_iterator1 tolabel   = elabel->second.begin();

    DFSCodeIsMin.push (0, 1, fromlabel->first, elabel->first, tolabel->first);

    return (project_is_min (tolabel->second, g, DFSCodeIsMin));
}

bool DFSCode::project_is_min (Projected &projected, Graph &g, DFSCode &DFSCodeIsMin)
// For gSpanSubs: function originally in the scope of gSpan
{
    const RMPath& rmpath = DFSCodeIsMin.buildRMPath ();
    int minlabel         = DFSCodeIsMin[0].fromlabel;
    int maxtoc           = DFSCodeIsMin[rmpath[0]].to;

    {
        Projected_map1 root;
        bool flg = false;
        int newto = 0;

        for (int i = rmpath.size()-1; ! flg  && i >= 1; --i)
        {
            for (unsigned int n = 0; n < projected.size(); ++n)
            {
                PDFS *cur = &projected[n];
                History history (g, cur);
                Edge *e = get_backward (g, history[rmpath[i]], history[rmpath[0]], history);
                if (e)
                {
                    root[e->elabel].push (0, e, cur);
                    newto = DFSCodeIsMin[rmpath[i]].from;
                    flg = true;
                }
            }
        }

        if (flg)
        {
            Projected_iterator1 elabel = root.begin();
            DFSCodeIsMin.push (maxtoc, newto, -1, elabel->first, -1);
            /* The equality test only depends on the start and end vertex positions as well as the edge label */
            unsigned int back_position = DFSCodeIsMin.size()-1;
            if ((*this)[back_position].from != DFSCodeIsMin[back_position].from
                    || (*this)[back_position].to != DFSCodeIsMin[back_position].to
                    || (*this)[back_position].elabel != DFSCodeIsMin[back_position].elabel)
                return false;
            return project_is_min (elabel->second, g, DFSCodeIsMin);
        }
    }

    {
        bool flg = false;
        int newfrom = 0;
        Projected_map2 root;
        EdgeList edges;

        for (unsigned int n = 0; n < projected.size(); ++n)
        {
            PDFS *cur = &projected[n];
            History history (g, cur);
            if (get_forward_pure (g, history[rmpath[0]], minlabel, history, edges))
            {
                flg = true;
                newfrom = maxtoc;
                for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it)
                    root[(*it)->elabel][g[(*it)->to].label].push (0, *it, cur);
            }
        }

        for (int i = 0; ! flg && i < (int)rmpath.size(); ++i)
        {
            for (unsigned int n = 0; n < projected.size(); ++n)
            {
                PDFS *cur = &projected[n];
                History history (g, cur);
                if (get_forward_rmpath (g, history[rmpath[i]], minlabel, history, edges))
                {
                    flg = true;
                    newfrom = DFSCodeIsMin[rmpath[i]].from;
                    for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it)
                        root[(*it)->elabel][g[(*it)->to].label].push (0, *it, cur);
                }
            }
        }

        if (flg)
        {
            Projected_iterator2 elabel  = root.begin();
            Projected_iterator1 tolabel = elabel->second.begin();
            DFSCodeIsMin.push (newfrom, maxtoc + 1, -1, elabel->first, tolabel->first);
            /* The equality test only depends on the start and end vertex positions as well as the edge and to label */
            unsigned int back_position = DFSCodeIsMin.size()-1;
            if ((*this)[back_position].from != DFSCodeIsMin[back_position].from
                    || (*this)[back_position].to != DFSCodeIsMin[back_position].to
                    || (*this)[back_position].elabel != DFSCodeIsMin[back_position].elabel
                    || (*this)[back_position].tolabel != DFSCodeIsMin[back_position].tolabel)
                return false;
            return project_is_min (tolabel->second, g, DFSCodeIsMin);
        }
    }

    return true;
}

DFSCode DFSCode::get_min_DFS()
// For gSpanSubs: compute minimum DFS Code
{
    Graph g;
    toGraph(g);
    EdgeList edges;
    Projected_map3 root;
    DFSCode minDFSCode;

    for (unsigned int from = 0; from < g.size() ; ++from)
        if (get_forward_root (g, g[from], edges))
            for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it)
            {
                root[g[from].label][(*it)->elabel][g[(*it)->to].label].push (0, *it, 0);
            }

    for (Projected_iterator3 fromlabel = root.begin() ;
            fromlabel != root.end() ; ++fromlabel)
        for (Projected_iterator2 elabel = fromlabel->second.begin() ;
                elabel != fromlabel->second.end() ; ++elabel)
            for (Projected_iterator1 tolabel = elabel->second.begin();
                    tolabel != elabel->second.end(); ++tolabel)
            {
                /* Build the initial two-node graph.  It will be grown
                *  recursively within the project.
                */
                minDFSCode.push (0, 1, fromlabel->first, elabel->first, tolabel->first);
                if (project_min_DFS (tolabel->second, g, minDFSCode))
                    return minDFSCode;
                minDFSCode.pop ();
            }
    std::cerr<< "Compute MinDFS Error (gSpan::get_min_DFS:"<< __LINE__<<")"<<std::endl;
    exit(1);
}

bool DFSCode::project_min_DFS(Projected &projected, Graph &g, DFSCode &minDFSCode)
/* For gSpanSubs: recursively compute min DFS Code
 * similar as project
 */
{
    if (minDFSCode.is_min()== false)
        return false;

    if (minDFSCode.edgeCount() == g.edge_size())
        return true;

    const RMPath &rmpath = minDFSCode.buildRMPath ();
    int minlabel = minDFSCode[0].fromlabel;
    int maxtoc = minDFSCode[rmpath[0]].to;

    Projected_map3 fwd_root;
    Projected_map2 bck_root;
    EdgeList edges;

    // Enumerate all possible one edge extensions for subs
    for (unsigned int n = 0; n < projected.size(); ++n) {
        // for all appearances in the graph
        PDFS *cur = &projected[n];
        History history (g, cur);

        for (int i = (int)rmpath.size()-1; i >= 1; --i)
        {
            Edge *e = get_backward (g, history[rmpath[i]], history[rmpath[0]], history);
            if (e)
                bck_root[minDFSCode[rmpath[i]].from][e->elabel].push (0, e, cur);
        }

        if (get_forward_pure (g, history[rmpath[0]], minlabel, history, edges))
            for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
                fwd_root[maxtoc][(*it)->elabel][g[(*it)->to].label].push (0, *it, cur);

        for (int i = 0; i < (int)rmpath.size(); ++i)
            if (get_forward_rmpath (g, history[rmpath[i]], minlabel, history, edges))
                for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it)
                    fwd_root[minDFSCode[rmpath[i]].from][(*it)->elabel][g[(*it)->to].label].push (0, *it, cur);
    }

    /* Test all extended substructures. */
    // backward
    for (Projected_iterator2 to = bck_root.begin(); to != bck_root.end(); ++to)
    {
        for (Projected_iterator1 elabel = to->second.begin(); elabel != to->second.end(); ++elabel)
        {
            minDFSCode.push (maxtoc, to->first, -1, elabel->first, -1);
            if(project_min_DFS (elabel->second, g, minDFSCode))
                return true;
            minDFSCode.pop();
        }
    }

    // forward
    for (Projected_riterator3 from = fwd_root.rbegin() ;
            from != fwd_root.rend() ; ++from)
        for (Projected_iterator2 elabel = from->second.begin() ;
                elabel != from->second.end() ; ++elabel)
            for (Projected_iterator1 tolabel = elabel->second.begin();
                    tolabel != elabel->second.end(); ++tolabel)
            {
                minDFSCode.push (from->first, maxtoc+1, -1, elabel->first, tolabel->first);
                if(project_min_DFS (tolabel->second, g, minDFSCode))
                    return true;
                minDFSCode.pop ();
            }
    return false;
}

bool DFSCodeTree::insert(DFSCode dfscode)
/* For gSpanSubs: add a new dfs code to the dfs code tree
 * if it is already in the tree, return false
 */
{
    DFSCodeTree *current, *newchild;
    current = this;
    for (unsigned int i = 0; i < dfscode.size(); i ++)
    {
        if (current -> children.count(dfscode[i]) != 0)
            current = current -> children[dfscode[i]];
        else
        {
            newchild = new DFSCodeTree();
            newchild -> curDFS = dfscode[i];
            newchild -> children.clear();
            current -> children[dfscode[i]] = newchild;
            current = newchild;
        }
    }
    if (current -> isPattern == true)
        return false;
    else
    {
        current -> isPattern = true;
        return true;
    }
}

void DFSCodeTree::show(unsigned int indent)
// For gSpanSubs: display the dfs code tree
{
    std::cout << std::string(indent, ' ');
    std::cout << '(' << (*this).curDFS.from << ' ' << (*this).curDFS.to << ' ' << (*this).curDFS.fromlabel << ' ' << (*this).curDFS.elabel << ' ' << (*this).curDFS.tolabel << ')' << std::endl;
    for(std::map<DFS, DFSCodeTree *>::iterator it = (*this).children.begin(); it != (*this).children.end(); ++it)
        (*(it -> second)).show(indent + 2);
}

IGraph DFSCode::get_IGraph()
// For gSpanSubs: get igraph structure for dfs code
{
    IGraph igraph;
    Graph g;
    (*this).toGraph(g);
    igraph_empty(&igraph.graph, 0, false);
    igraph_vector_long_init(&igraph.vertex_labels, 0);
    igraph_vector_long_init(&igraph.edge_labels, 0);
    for(unsigned int i = 0; i < g.size(); i ++) {
        igraph_add_vertices(&igraph.graph, 1, NULL);
        igraph_vector_long_push_back(&igraph.vertex_labels, (long int)g[i].label);
    }
    for(unsigned int i = 0; i < g.size(); i ++)
        for(Vertex::edge_iterator it = g[i].edge.begin();
                it != g[i].edge.end(); it ++) {
            igraph_add_edge(&igraph.graph, (long int)(it -> from), (long int)(it -> to));
            igraph_vector_long_push_back(&igraph.edge_labels, (long int)(it -> elabel));
        }
    return igraph;
}

bool DFSCode::cover_all_subs(std::map< unsigned int, IGraph > &subs_set)
// For gSpanSubs: test if current pattern cover all cores
{
    IGraph ig = get_IGraph();
    igraph_bool_t result;
    for (unsigned int i = back().cover_subs_num; i < subs_set.size(); i ++)
    {
        igraph_subisomorphic_vf2(
                &ig.graph,
                &(subs_set[i].graph),
                &(ig.vertex_labels),
                &(subs_set[i].vertex_labels),
                &(ig.edge_labels),
                &(subs_set[i].edge_labels),
                &result, NULL, NULL);
        if (result)
            back().cover_subs_num ++;
        else
            break;
    }
    if (back().cover_subs_num == subs_set.size())
        return true;
    else
        return false;
}
}
