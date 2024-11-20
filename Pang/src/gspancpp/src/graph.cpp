/*
    gSpanSubs 2011/03/08 implemented by Lei Zhao based on:

    graph.cpp, v 1.0 2009/12/21 by Marisa Thoma

    is a modification based on :

    $Id: graph.cpp,v 1.4 2004/05/21 05:50:13 taku-ku Exp $;

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

     Modifications: excluded from matlab, handled some IO issues
*/

#include "gspan.h"

#include <cstring>
#include <assert.h>

#include <string>
#include <iterator>
#include <strstream>
#include <set>

namespace GSPAN {

template <class T, class Iterator>
void tokenize (const char *str, Iterator iterator)
{
    std::istrstream is (str, std::strlen(str));
    std::copy (std::istream_iterator <T> (is), std::istream_iterator <T> (), iterator);
}

void Graph::buildEdge ()
{
    char buf[512];
    std::map <std::string, unsigned int> tmp;

    unsigned int id = 0;
    for (int from = 0; from < (int)size (); ++from) {
        for (Vertex::edge_iterator it = (*this)[from].edge.begin ();
                it != (*this)[from].edge.end (); ++it)
        {
            if (directed || from <= it->to)
                std::sprintf (buf, "%d %d %d", from, it->to, it->elabel);
            else
                std::sprintf (buf, "%d %d %d", it->to, from, it->elabel);

            // Assign unique id's for the edges.
            if (tmp.find (buf) == tmp.end()) {
                it->id = id;
                tmp[buf] = id;
                ++id;
            } else {
                it->id = tmp[buf];
            }
        }
    }

    edge_size_ = id;
}

std::istream &Graph::read (std::istream &is, std::string * gName)
{
    std::vector <std::string> result;
    char line[1024];

    clear ();

    while (is.peek() != 't') {

        if (is.peek() == EOF)
            return is;

        is.getline (line, 1024);
    }
    if (! is.getline (line, 1024)) {
        std::cerr << "error: no graph entry here at all Graph::read:" <<__LINE__  << std::endl;
        exit(1);
    }
    if (line[1] != ' ' || line[2] != '#') {
        std::cerr << "error: illegal graph identifier \""<<line<<"\""  << std::endl;
        std::cerr << "       needs DIMACS format: \"t # id\" (Graph::read:" <<__LINE__<<")"  << std::endl;
        exit(1);
    }
    std::string graph_name(line+4);

    if (gName != NULL)
        *gName = graph_name;

    while (is.peek() != 't') {


        if (! is.getline (line, 1024))
            break;

        result.clear();
        tokenize<std::string>(line, std::back_inserter (result));

        if (result.empty()) {
            // do nothing
        } else if (result[0] == "t") {

            std::cerr << "error: should not arrive here: Graph::read:"<<__LINE__ << std::endl;
            exit(1);

            if (! empty()) { // use as delimiter

            } else {
                /*
                 * y = atoi (result[3].c_str());
                 */
            }
        } else if (result[0] == "v" && result.size() >= 3) {
            unsigned int id    = atoi (result[1].c_str());
            this->resize (id + 1);
            (*this)[id].label = atoi (result[2].c_str());
        } else if (result[0] == "e" && result.size() >= 4) {
            int from   = atoi (result[1].c_str());
            int to     = atoi (result[2].c_str());
            int elabel = atoi (result[3].c_str());

            if ((int)size () <= from || (int)size () <= to) {
                std::cerr << "Format Error:  define vertex lists before edges (graph '"<<graph_name<<"', ["<<from<<", "<<to<<", "<<elabel<<"]" << std::endl;
                exit (-1);
            }

            (*this)[from].push (from, to, elabel);
            if (directed == false)
                (*this)[to].push (to, from, elabel);
        }
    }

    buildEdge ();

    return is;
}

std::ostream &Graph::write (std::ostream &os) const
{
    char buf[512];
    std::set <std::string> tmp;

    for (int from = 0; from < (int)size (); ++from) {
        os << "v " << from << " " << (*this)[from].label << std::endl;

        Vertex::const_edge_iterator theEnd = (*this)[from].edge.end();

        for (Vertex::const_edge_iterator it = (*this)[from].edge.begin ();
                it != theEnd; ++it) {

            if (directed || from <= it->to) {
                std::sprintf (buf, "%d %d %d", from, it->to,   it->elabel);
            } else {
                std::sprintf (buf, "%d %d %d", it->to,   from, it->elabel);
            }
            tmp.insert (buf);
        }
    }

    for (std::set<std::string>::iterator it = tmp.begin(); it != tmp.end(); ++it) {
        os << "e " << *it << std::endl;
    }

    return os;
}

void Graph::check (void)
{
    /* Check all indices */
    for (int from = 0 ; from < (int)size () ; ++from) {

        for (Vertex::edge_iterator it = (*this)[from].edge.begin ();
                it != (*this)[from].edge.end (); ++it)
        {
            assert (it->from >= 0 && it->from < static_cast<int>(size ()));
            assert (it->to >= 0 && it->to < static_cast<int>( size ()));
        }
    }
}

IGraph Graph::get_IGraph()
// For gSpanSubs: return igraph structure for graph
{
    IGraph igraph;
    igraph_empty(&igraph.graph, 0, false);
    igraph_vector_long_init(&igraph.vertex_labels, 0);
    igraph_vector_long_init(&igraph.edge_labels, 0);
    for(unsigned int i = 0; i < size(); i ++) {
        igraph_add_vertices(&igraph.graph, 1, NULL);
        igraph_vector_long_push_back(&igraph.vertex_labels, (long int)(*this)[i].label);
    }
    for(unsigned int i = 0; i < size(); i ++)
        for(Vertex::edge_iterator it = (*this)[i].edge.begin();
                it != (*this)[i].edge.end(); it ++) {
            igraph_add_edge(&igraph.graph, (long int)(it -> from), (long int)(it -> to));
            igraph_vector_long_push_back(&igraph.edge_labels, (long int)(it -> elabel));
        }
    return igraph;
}

}

