/*
    gSpanSubs 2011/03/08 implemented by Lei Zhao based on:

    gspan.cpp, v 1.0 2009/12/21 by Marisa Thoma

    is a modification based on :

    $Id: gspan.cpp,v 1.8 2004/05/21 09:27:17 taku-ku Exp $;

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

   Mmodified by Marisa Thoma (extracted from matlab; extended to CORK pruning)

*/
#include "gspan.h"

#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <limits>

#include <iterator>
#include <iostream>
#include <cassert>

namespace GSPAN {

gSpan::gSpan (void) {
    infoStream = &std::cerr;
}

std::istream &gSpan::read (std::istream &is,
                           std::set<unsigned int> const * filter) {

    Graph g(directed);
    std::string gName;
    unsigned int gNameInt;
    while (true) {
        g.read (is, &gName);
        gNameInt = strtoul(gName.c_str(),NULL,10);
        if (g.empty()) break;
        if (filter != NULL) { // select for elements in filter
            if (filter->find(gNameInt) == filter->end()) // graph not selected
                continue;
        }
        TRANS.push_back (g);
        // assign class labels if possible
        if (instances2classLabels.size() != 0) {
            std::map<unsigned int, unsigned int>::const_iterator clIter = instances2classLabels.find(gNameInt);
            if (clIter == instances2classLabels.end()) {
                std::cerr <<"error: instance \""<<gNameInt<<"\" is not covered in class label map (gSpan::read:"<<__LINE__<<")"<<std::endl;
                exit(1);
            }
            class_labels.push_back(clIter->second);
            originalInstanceLabels.push_back(gNameInt);
        }
    }

    return is;
}

std::map<unsigned int, unsigned int> gSpan::support_counts (Projected &projected) {
    std::map<unsigned int, unsigned int> counts;

    for (Projected::iterator cur = projected.begin() ;
            cur != projected.end() ; ++cur)
    {
        counts[cur->id] ++;
    }
    return counts;
}

unsigned int gSpan::support (Projected &projected) {
    unsigned int oid = 0xffffffff;
    unsigned int size = 0;

    for (Projected::iterator cur = projected.begin(); cur != projected.end(); ++cur) {
        if (oid != cur->id) {
            ++size;
        }
        oid = cur->id;
    }

    return size;
}

/* Special report function for single node graphs. -- NOT supported by FS
 */
void gSpan::report_single (Graph &g,
                           std::map<unsigned int, unsigned int>& ncount) {
    unsigned int sup = 0;
    for (std::map<unsigned int, unsigned int>::iterator it = ncount.begin () ;
            it != ncount.end () ; ++it) {
        sup += (*it).second;
    }

    if (maxpat_max > maxpat_min && g.size () > maxpat_max)
        return;
    if (maxpat_min > 0 && g.size () < maxpat_min)
        return;
    if (enc == false) {
        *os << "t # " << ID << " * " << sup;
        *os << '\n';

        g.write (*os);
        *os << '\n';
    } else {
        std::cerr << "report_single not implemented for non-Matlab calls" << std::endl;
    }

}

void gSpan::report (Projected &projected, unsigned int sup) {

    /* Filter too small / too large graphs.
     */
    if (maxpat_max > maxpat_min && DFS_CODE.nodeCount () > maxpat_max)
        return;
    if (maxpat_min > 0 && DFS_CODE.nodeCount () < maxpat_min)
        return;

    if (xml && where) {
        *os << "<pattern>\n";
        *os << "<id>" << ID << "</id>\n";
        *os << "<support>" << sup << "</support>\n";
        *os << "<what>";
    }

    if (! enc) {
        Graph g(directed);
        DFS_CODE.toGraph (g);

        if (! xml && where) // separate
            *os << std::endl;

        *os << "t # " << ID << " * " << sup;
        *os << '\n';
        g.write (*os);
    } else {
        if (! xml || ! where)
            *os << '<' << ID << ">    " << sup << " [";

        DFS_CODE.write (*os);
        if (! xml || ! where)
            *os << ']';
    }

    if ((bool)where) { // list graph ids for pattern
        if (xml) // close tag
            *os << "</what>\n<where>";
        else
            *os << " {";

        unsigned int oid = 0xffffffff;
        if (where == 1) { // => only report graph ids as traceback
            for (Projected::iterator cur = projected.begin(); cur != projected.end(); ++cur) {
                if (oid != cur->id) {
                    if (cur != projected.begin())
                        *os << ' ';
                    if (originalInstanceLabels.size() == 0)
                        *os << cur->id;
                    else
                        *os << originalInstanceLabels.at(cur->id);
                }
                oid = cur->id;
            }
        } else { // must be 2 => also report hit frequencies
            unsigned int freqCount = 0;
            for (Projected::iterator cur = projected.begin(); cur != projected.end(); ++cur) {
                if (oid != cur->id) {
                    if (cur != projected.begin()) { // freqCount == 0
                        *os << ':' << freqCount << ' ';
                        freqCount = 0;
                    }
                    if (originalInstanceLabels.size() == 0)
                        *os << cur->id;
                    else
                        *os << originalInstanceLabels.at(cur->id);

                    oid = cur->id;
                }
                freqCount++;
            }
            *os << ':' << freqCount; // last frequency
        }
        if (xml)
            *os << "</where>\n</pattern>";
        else
            *os << '}';
    }

    *os << '\n';
    ++ID;
}

void gSpan::report (DFSCode & toReport,
                    std::map<unsigned int, unsigned int> & supp_counts) {

    /* Filter too small / too large graphs.
     */
    if (maxpat_max > maxpat_min && toReport.nodeCount () > maxpat_max)
        return;
    if (maxpat_min > 0 && toReport.nodeCount () < maxpat_min)
        return;

    if (xml && where) {
        *os << "<pattern>\n";
        *os << "<id>" << ID << "</id>\n";
        *os << "<support>" << supp_counts.size() << "</support>\n";
        *os << "<what>";
    }

    if (! enc) {
        Graph g(directed);
        toReport.toGraph (g);

        if (! xml && where) // separate
            *os << std::endl;

        *os << "t # " << ID << " * " << supp_counts.size();
        *os << std::endl;
        g.write (*os);
    } else {
        if (! xml || ! where)
            *os << '<' << ID << ">    " << supp_counts.size() << " [";

        toReport.write (*os);
        if (! xml || ! where)
            *os << ']';
    }

    if ((bool)where) { // list graph ids for pattern
        if (xml) // close tag
            *os << "</what>\n<where>";
        else
            *os << " {";

        if (where == 1)
            for (std::map<unsigned int, unsigned int>::const_iterator cur_iter = supp_counts.begin(); cur_iter != supp_counts.end(); ++cur_iter) {
                assert( originalInstanceLabels.size() == 0
                        || cur_iter->first < originalInstanceLabels.size() );

                if (cur_iter != supp_counts.begin())
                    *os << ' ';
                *os << (originalInstanceLabels.size() == 0 ?
                        cur_iter->first :
                        originalInstanceLabels.at(cur_iter->first));
            }
        else { // must be 2
            for (std::map<unsigned int, unsigned int>::const_iterator cur_iter = supp_counts.begin(); cur_iter != supp_counts.end(); ++cur_iter) {
                assert( originalInstanceLabels.size() == 0
                        || cur_iter->first < originalInstanceLabels.size() );

                if (cur_iter != supp_counts.begin())
                    *os << ' ';
                *os <<  (originalInstanceLabels.size() == 0 ?
                         cur_iter->first :
                         originalInstanceLabels.at(cur_iter->first)) <<
                    ':' << cur_iter->second;
            }
        }
        if (xml)
            *os << "</where>\n</pattern>";
        else
            *os << '}';
    }

    *os << '\n';
    ++ID;
}

/**
  * Get the maximally possible improvement achieved by the extensions stored in
  * correspondence class cor.
  */
unsigned int getPartialImprovment(const Correspondence * cor) 
{
    // improvement = max { 0,     A1 * (B1 - B0),     B1 * (A1 - B0) } (for all classes; 1-vs.-rest)
    unsigned int tempImprovement = 0, s = (cor->counts).size(), subImprovement;
    int a_bdiff, b_adiff;
    unsigned int total_counts = 0, total_ext = 0;
    for (unsigned int i = 0; i < s-1; i++) {
        total_counts += (cor->counts)[i];
        total_ext += (cor->ext)[i];
    }
    for (unsigned int i = 0; i < s-1; i++) { // class A  ---   B is the rest of the dataset
        subImprovement = 0; // no improvement
        // A1 * (B1 - B0) = A1 * (B1 - (B - B1)) = A1 * (2 B1 - B)
        a_bdiff = (cor->ext)[i] * (2*(total_ext - (cor->ext)[i]) - (total_counts - (cor->counts)[i]));
        if (a_bdiff > 0)
            subImprovement = static_cast<unsigned int>( a_bdiff );
        // B1 * (A1 - A0) = B1 * (A1 - (A - A1)) = B1 * (2 A1 - A)
        b_adiff = (total_ext - (cor->ext)[i]) * (2*(cor->ext)[i] - (cor->counts)[i]);
        if (b_adiff > static_cast<int>(subImprovement))
            subImprovement = static_cast<unsigned int>( b_adiff );
        tempImprovement += subImprovement;
    }
    return tempImprovement;
}

/**
  * Multi-class correspondences are:
  * sum_{class a} #correspondences_{"a" vs. "not a"}
  */
unsigned int oneAgainstRestCorresondences(std::vector<unsigned int> const & classes) {
    unsigned int corrs = 0, numMatches = 0;
    for (std::vector<unsigned int>::const_iterator cIt = classes.begin(); cIt != classes.end(); cIt++)
        numMatches += *cIt;
    for (std::vector<unsigned int>::const_iterator cIt = classes.begin(); cIt != classes.end(); cIt++) {
        corrs += *cIt * (numMatches - *cIt);
    }
    return corrs;
}

/*
 * Returns true if any child graph of the current subgraph represented in
 *   "projected" (or, if "use_conserved", in "current_best_support_counts")
 *   can improve the current CORK-value by more than the current
 *   maximal CORK value - if "get_cork_value" != NULL, it gets assigned the
 *   new CORK value of the current subgraph.
 */
bool gSpan::getAskCORK(Projected const * projected,
                       unsigned int * get_cork_value,
                       bool use_conserved) {
    /* identify the effects of all matching graphs on the correspondence
         * equivalence classes */
    unsigned int oid = 0xffffffff;
    Correspondence * cor;
    if (use_conserved) {
        for (std::map<unsigned int, unsigned int>::const_iterator c_iter = current_best_support_counts.begin(); c_iter != current_best_support_counts.end(); c_iter++) {
            assert(equ_classes.size() > c_iter->first && class_labels.size() > c_iter->first && equ_classes.size() > c_iter->first && possible_changes.size() > c_iter->first && correspondence_classes.size() > equ_classes[c_iter->first]);
            cor = &(correspondence_classes[equ_classes[c_iter->first]]);
            (cor->ext).resize((cor->counts).size()); // ensure ext is initialized
            (cor->ext)[class_labels[c_iter->first]]++;
            possible_changes[c_iter->first] = true; // prepare extend() step
        }
    } else {
        for (Projected::const_iterator cur = projected->begin(); cur != projected->end(); ++cur) {
            if (oid != cur->id) {
                // std::cout<<"ec.s()="<<equ_classes.size()<<", cc.s()="<<correspondence_classes.size()<<", pc.s()="<<possible_changes.size()<<", cur->id="<<cur->id<<", equ_classes[cur->id]="<<equ_classes[cur->id]<<", lab="<<class_labels[cur->id]<<std::endl;
                assert(equ_classes.size() > cur->id && class_labels.size() > cur->id && equ_classes.size() > cur->id && possible_changes.size() > cur->id && correspondence_classes.size() > equ_classes[cur->id]);
                cor = &(correspondence_classes[equ_classes[cur->id]]);
                // std::cout << "cor: "<<cor->father<<", cs="<<(cor->counts).size()<<", ex="<<(cor->ext).size() << std::endl;
                (cor->ext).resize((cor->counts).size()); // ensure ext is initialized
                (cor->ext)[class_labels[cur->id]]++;
                oid = cur->id;
            }
        }
    }

    if (get_cork_value == NULL) // only re-freshing extension information - already passed CORK test
        return true;

    /* calculate current CORK value and possible further gain */

    *get_cork_value = 0;
    unsigned int maxImprovement = 0, tempImprovement = 0;
    for (std::vector<Correspondence>::iterator cor_iter = correspondence_classes.begin(); cor_iter != correspondence_classes.end(); cor_iter++) {
        //  compute CORK
        if ((cor_iter->ext).size() == 0) {
            (cor_iter->ext).resize((cor_iter->counts).size());
        }
        assert((cor_iter->ext).size() == (cor_iter->counts).size());

        std::vector<unsigned int> misses = std::vector<unsigned int>(cor_iter->counts);
        std::vector<unsigned int>::const_iterator hIt = (cor_iter->ext).begin();
        for (std::vector<unsigned int>::iterator mIt = misses.begin(); mIt != misses.end(); mIt++) {
            *mIt -= *hIt;
            hIt++;
        }
        *get_cork_value += oneAgainstRestCorresondences(cor_iter->ext);
        *get_cork_value += oneAgainstRestCorresondences(misses);

        tempImprovement = getPartialImprovment(&(*cor_iter));

        maxImprovement += tempImprovement;
    }

    /* test CORK pruning threshold */
    assert(maxImprovement <= *get_cork_value);
    if (*get_cork_value - maxImprovement < correspondences)
        return true;
    return false;
}

void gSpan::resetCORK_Extensions() {
    for (std::vector<Correspondence>::iterator cor_iter = correspondence_classes.begin(); cor_iter != correspondence_classes.end(); cor_iter++) {
        (cor_iter->ext).assign((cor_iter->counts).size(),0);
    }
}

void gSpan::extendCORK() {
    std::vector<unsigned int> correspondence_children(correspondence_classes.size());

    unsigned int numCs = correspondence_classes.size();
    for (unsigned int i = 0; i < numCs; i++) {
        assert(i < correspondence_classes.size());

        unsigned int num0 = 0, extNot0 = 0, numEqMax = 0;
        for(unsigned int j = 0; j < correspondence_classes[i].counts.size(); j++) {
            if (correspondence_classes[i].ext[j] == correspondence_classes[i].counts[j])
                numEqMax++;
            if (correspondence_classes[i].ext[j] != 0)
                extNot0++;
            if (correspondence_classes[i].counts[j] == 0)
                num0++;
        }

        if (extNot0 > 0) {
            if (numEqMax == correspondence_classes[i].counts.size()) {
                correspondence_classes[i].ext.assign(correspondence_classes[i].counts.size(), 0);
                continue; // not a really new equivalence class
            }

            /* possibility of avoiding meaningless correspondence classes: */
            if (delete_resolved)
                if (num0 >= correspondence_classes[i].counts.size()-1) {
                    // split up only slows us down; however this means that we will not know the correct number of patterns generated by the selected subgraphs
                    correspondence_classes[i].ext.assign(correspondence_classes[i].counts.size(), 0);
                    continue;
                }
            /*****/

            correspondence_children[i] = correspondence_classes.size();

            // build and fill new correspondence
            correspondence_classes.resize(correspondence_classes.size()+1);
            correspondence_classes[correspondence_classes.size()-1].counts = correspondence_classes[i].ext;
            correspondence_classes[correspondence_classes.size()-1].father = i;

            // reduce "old" correspondence
            for (unsigned int j = 0; j < correspondence_classes[i].counts.size(); j++)
                correspondence_classes[i].counts[j] -= correspondence_classes[i].ext[j];
            correspondence_classes[i].ext.assign(correspondence_classes[i].counts.size(), 0);
        }
    }

    // assign new correspondence equivalence classes
    for (unsigned int i = 0; i < equ_classes.size(); i++) {
        if (possible_changes[i]) {
            if (correspondence_children[equ_classes[i]] != 0)
                equ_classes[i] = correspondence_children[equ_classes[i]];
            possible_changes[i] = false;
        }
    }

    // delete all graphs which are not further participating in feature selection
    if (delete_resolved) {
        std::vector<unsigned int>::iterator clIter = class_labels.begin();
        std::vector<unsigned int>::iterator ecIter = equ_classes.begin(), oilIter = originalInstanceLabels.begin(), oirIter = originalInstanceRanks.begin();
        std::vector < Graph >::iterator tIter = TRANS.begin();
        while (ecIter != equ_classes.end()) {
            unsigned int num0 = 0;
            for(unsigned int j = 0; j < correspondence_classes[*ecIter].counts.size(); j++) {
                if (correspondence_classes[*ecIter].counts[j] == 0)
                    num0++;
            }
            if (num0 >= correspondence_classes[*ecIter].counts.size()-1) { // resolved
                ecIter = equ_classes.erase(ecIter);
                clIter = class_labels.erase(clIter);
                tIter = TRANS.erase(tIter);
                if (oilIter != originalInstanceLabels.end())
                    oilIter = originalInstanceLabels.erase(oilIter);
                if (oirIter != originalInstanceRanks.end())
                    oirIter = originalInstanceRanks.erase(oirIter);
            } else {
                ecIter++;
                clIter++;
                tIter++;
                if (oilIter != originalInstanceLabels.end())
                    oilIter++;
                if (oirIter != originalInstanceRanks.end())
                    oirIter++;
            }
        }
        possible_changes.resize(TRANS.size());
    }

}

/* Recursive subgraph mining function (similar to subprocedure 1
 * Subgraph_Mining in [Yan2002]).
 */
void gSpan::project (Projected &projected) {

    tested_subgraphs++;
    bool tempresult; // For gSpanSubs: performance

    /* Check if the pattern is frequent enough.
     */
    start_time = clock();
    unsigned int sup = support (projected);
    gspan_filter_sup_time = gspan_filter_sup_time + clock() - start_time;
    if (sup < minsup)
        return;

    frequent_subgraphs++;

    if (verbose)
        std::cout<<"frequent ("<<sup<<")"<<std::endl;


    /* The minimal DFS code check is more expensive than the support check,
     *  hence it is done now, after checking the support.
     */
    start_time = clock();
    tempresult = DFS_CODE.is_min();
    gspan_filter_dup_time = gspan_filter_dup_time + clock() - start_time;

    if (tempresult == false) // For gSpanSubs
        return;

    if (verbose)
        std::cout<<"minimal"<<std::endl;

    minimal_subgraphs++;

    if (fs) { // exploit feature selection pruning

        assert(equ_classes.size() == class_labels.size() && class_labels.size() == possible_changes.size());

        if (1) { // enable further selection options

            unsigned int corrs = 0;
            bool testCORK = getAskCORK(&projected, &corrs);

            if (verbose) {
                std::cout<<corrs<<" correspondences"<<std::endl;
            }

            if ( corrs < correspondences ) { // found better subgraph
                correspondences = corrs;
                if (verbose)
                    std::cout<<"new winner"<<std::endl;

                winner_subgraphs++;

                // TODO: more efficiency via explicit vector copy calls?
                current_best_support_counts = support_counts(projected);
                current_best_subgraph = DFS_CODE;
                // might be avoided ...
            }

            // TODO: make more efficient - using projected
            resetCORK_Extensions(); // reset to zero-extensions

            if (verbose)
                std::cout<<"reset"<<std::endl;

            if (!testCORK) { // failed CORK-test
                if (verbose)
                    std::cout<<"failed test"<<std::endl;
                prunedFS_subgraphs++;
                return; // do not extend
            }

            if (verbose)
                std::cout<<"passed test"<<std::endl;

        }
    }
    else // output is delayed or omitted
    {
        if (report_result) // For gSpanSubs
            report (projected, sup); // Output the frequent substructure
        pattern_num ++;
    }

    /* In case we have a valid upper bound and our graph already exceeds it,
     * return.  Note: we do not check for equality as the DFS exploration may
     * still add edges within an existing subgraph, without increasing the
     * number of nodes.
     */
    if (maxpat_max > maxpat_min && DFS_CODE.nodeCount () > maxpat_max)
        return;

    /* We just outputted a frequent subgraph.  As it is frequent enough, so
     * might be its (n+1)-extension-graphs, hence we enumerate them all.
     */
    start_time = clock(); // For gSpanSubs
    const RMPath &rmpath = DFS_CODE.buildRMPath ();
    int minlabel = DFS_CODE[0].fromlabel;
    int maxtoc = DFS_CODE[rmpath[0]].to;

    Projected_map3 new_fwd_root;
    Projected_map2 new_bck_root;
    EdgeList edges;
    gspan_edge_grow_time = gspan_edge_grow_time + clock() - start_time; // For gSpanSubs
    /* Enumerate all possible one edge extensions of the current substructure.
     */
    for (unsigned int n = 0; n < projected.size(); ++n) { // for all fitting graphs

        unsigned int id = projected[n].id;
        PDFS *cur = &projected[n];
        History history (TRANS[id], cur);

        // XXX: do we have to change something here for directed edges?
        
        start_time = clock();// For gSpanSubs
        // backward
        for (int i = (int)rmpath.size()-1; i >= 1; --i) {
            Edge *e = get_backward (TRANS[id], history[rmpath[i]], history[rmpath[0]], history);
            if (e)
                new_bck_root[DFS_CODE[rmpath[i]].from][e->elabel].push (id, e, cur);
        }

        // pure forward
        // FIXME: here we pass a too large e->to (== history[rmpath[0]]->to
        // into get_forward_pure, such that the assertion fails.
        //
        // The problem is:
        // history[rmpath[0]]->to > TRANS[id].size()
        if (get_forward_pure (TRANS[id], history[rmpath[0]], minlabel, history, edges))
            for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
                new_fwd_root[maxtoc][(*it)->elabel][TRANS[id][(*it)->to].label].push (id, *it, cur);

        // backtracked forward
        for (int i = 0; i < (int)rmpath.size(); ++i)
            if (get_forward_rmpath (TRANS[id], history[rmpath[i]], minlabel, history, edges))
                for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it)
                    new_fwd_root[DFS_CODE[rmpath[i]].from][(*it)->elabel][TRANS[id][(*it)->to].label].push (id, *it, cur);
        gspan_edge_grow_time = gspan_edge_grow_time + clock() - start_time;
    }

    /* Test all extended substructures.
     */
    // backward
    for (Projected_iterator2 to = new_bck_root.begin(); to != new_bck_root.end(); ++to) {
        for (Projected_iterator1 elabel = to->second.begin(); elabel != to->second.end(); ++elabel) { // mis-used parameter names ...
            DFS_CODE.push (maxtoc, to->first, -1, elabel->first, -1);
            project (elabel->second);
            DFS_CODE.pop();
        }
    }

    // forward
    for (Projected_riterator3 from = new_fwd_root.rbegin() ;
            from != new_fwd_root.rend() ; ++from)	{
        for (Projected_iterator2 elabel = from->second.begin() ;
                elabel != from->second.end() ; ++elabel) {
            for (Projected_iterator1 tolabel = elabel->second.begin();
                    tolabel != elabel->second.end(); ++tolabel)	{
                DFS_CODE.push (from->first, maxtoc+1, -1, elabel->first, tolabel->first);
                project (tolabel->second);
                DFS_CODE.pop ();
            }
        }
    }

    return;
}


std::vector<unsigned int> gSpan::assign_class_labels(std::string class_label_file, std::set<unsigned int> const * graphFilter) {
    if (TRANS.size() == 0) {
        std::cerr<< "error: no graphs read in before assigning graph class labels (gSpan::assign_class_labels:"<< __LINE__<<")"<<std::endl;
        exit(1);
    }
    std::vector<unsigned int> class_counts;
    if (class_label_file=="") {
        // class labels must already have been set;
        // only assert dimensions and update class counts
        if (class_labels.size() != TRANS.size()) {
            std::cerr<< "error: # mapped class labels ("<<class_labels.size()<<") must be equal zu #graphs ("<<TRANS.size()<<") (gSpan::assign_class_labels:"<< __LINE__<<")"<<std::endl;
            exit(1);
        }
        std::set<unsigned int> disjunctClasses;
        unsigned int max_label = 0;
        for (unsigned int i = 0; i < class_labels.size(); i++) {
            if (class_counts.size() <= class_labels[i])
                class_counts.resize(class_labels[i]+1);
            class_counts[class_labels[i]]++;
            disjunctClasses.insert(class_labels[i]);
            if (max_label < class_labels[i])
                max_label = class_labels[i];
        }
        number_of_classes = disjunctClasses.size();
        if (max_label != number_of_classes-1) {
            std::cerr<< "error: require a class labelling from 0 to (num_classes-1) " << max_label <<"!=" << (number_of_classes-1) <<"(gSpan::assign_class_labels:"<< __LINE__<<")"<<std::endl;
            exit(1);
        }
        return class_counts;
    }

    class_labels.resize(TRANS.size());

    std::ifstream classLabelFile(class_label_file.c_str());
    if (!classLabelFile) {
        std::cerr<< "error: cannot open file containing graph labels \""<<class_label_file<<"\" (gSpan::assign_class_labels:"<< __LINE__<<")"<<std::endl;
        exit(1);
    }
    std::string str, tempStr;

    getline(classLabelFile, str);
    std::istringstream iss(str);
    unsigned int graph_id = 0, /* accepted */ graph_index = 0 /* encountered */;
    if (!(iss>>tempStr)) { // take first word of the line
        tempStr=str; // one entry only
        if (tempStr=="") {
            std::cerr << "error: graph id mapping file \"" << class_label_file << "\" does not contain" << std::endl;
            std::cerr << "       a class label for graph " << graph_id << " (gSpan::assign_class_labels:" <<  __LINE__ << ")" << std::endl;
            exit(1);
        }
    }
    if (tempStr[0] == '-') {
        std::cerr << "error: graph id mapping file \""<<class_label_file<<"\" has a negative"<<std::endl;
        std::cerr << "       class label \""<<tempStr<<"\" for graph "<<graph_index<<" (gSpan::assign_class_labels:"<< __LINE__<<")"<<std::endl;
        exit(1);
    }
    unsigned int classLabel = strtoul(tempStr.c_str(),NULL,10);
    if (graphFilter == NULL || graphFilter->find(graph_index) != graphFilter->end()) {
        class_labels[graph_id] = classLabel;
        class_counts.resize(classLabel+1);
        class_counts[classLabel]++;
        graph_id++;
    }
    graph_index++;

    std::set<unsigned int> disjunctClasses;
    unsigned int max_label = 0;

    while(getline(classLabelFile, str)) {
        if (graphFilter != NULL && graphFilter->find(graph_index) == graphFilter->end()) {
            graph_index++; // skip this graph
            continue;
        }
        iss.str(str);
        if (!(iss>>tempStr)) { // take first word of the line
            tempStr=str; // one entry only
            if (tempStr=="") {
                std::cerr<< "error: graph id mapping file \""<<class_label_file<<"\" does not contain a class label for graph "<<graph_index<<" (gSpan::assign_class_labels:"<< __LINE__<<")"<<std::endl;
                std::cerr<< "       '"<<str<<"' to '"<<tempStr<<"' from '"<<iss.str()<<"'"<<std::endl;
                exit(1);
            }
        }
        if (tempStr[0] == '-') {
            std::cerr << "error: graph id mapping file \""<<class_label_file<<"\" has a negative"<<std::endl;
            std::cerr << "       class label \""<<tempStr<<"\" for graph "<<graph_index<<" (gSpan::assign_class_labels:"<< __LINE__<<")"<<std::endl;
            exit(1);
        }
        classLabel = strtoul(tempStr.c_str(), NULL, 10);
        class_labels[graph_id++] = classLabel;
        if (class_counts.size() <= classLabel)
            class_counts.resize(classLabel+1);
        class_counts[classLabel]++;
        disjunctClasses.insert(classLabel);
        if (max_label < classLabel)
            max_label = classLabel;
        graph_index++;
    }
    classLabelFile.close();

    number_of_classes = disjunctClasses.size();
    if (max_label != disjunctClasses.size()-1) {
        // must re-label the class labels (empty class slots cost time)
        unsigned int new_labels[max_label+1];
        unsigned int class_index = 0;
        for (std::set<unsigned int>::const_iterator clIt = disjunctClasses.begin(); clIt != disjunctClasses.end(); clIt++) {
            class_counts[class_index] = class_counts[*clIt];
            new_labels[*clIt] = class_index++;
        }
        for (std::vector<unsigned int>::iterator clIt = class_labels.begin(); clIt != class_labels.end(); clIt++) {
            class_index = new_labels[*clIt];
            *clIt = class_index;
        }
        class_counts.resize(disjunctClasses.size());
    }
    return class_counts;
}

std::map<unsigned int,unsigned int> gSpan::translateSupportCounts(std::map<unsigned int,unsigned int> const & _instances2classLabels) {
    if (originalInstanceRanks.size() == 0) {
        unsigned int i = 0;
        originalInstanceRanks.resize(originalInstanceLabels.size());
        std::map<unsigned int,unsigned int> lab2rank;
        for (std::map<unsigned int,unsigned int>::const_iterator i2cIt = instances2classLabels.begin(); i2cIt != instances2classLabels.end(); i2cIt++)
            lab2rank[i2cIt->first] = i++;
        for (i = 0; i < originalInstanceLabels.size(); i++)
            originalInstanceRanks.at(i) = lab2rank[originalInstanceLabels[i]];
    }
    std::map<unsigned int,unsigned int> supports_translated;
    for (std::map<unsigned int,unsigned int>::const_iterator i2cIt = _instances2classLabels.begin(); i2cIt != _instances2classLabels.end(); i2cIt++)
        supports_translated[originalInstanceRanks.at(i2cIt->first)] = i2cIt->second;

    return supports_translated;
}

void gSpan::setInfoStream(std::ostream &infoS) {
    infoStream = &infoS;
}

void gSpan::run (std::istream &is, std::ostream &_os,
                 unsigned int _minsup,
                 unsigned int _maxpat_min, unsigned int _maxpat_max,
                 bool _enc,
                 unsigned int _where,
                 bool _xml,
                 bool _directed, std::set<unsigned int> const * graphFilter,
                 bool _report_result, // For gSpanSubs
                 bool _verbose) {

    os = &_os;
    ID = 0;
    minsup = _minsup;
    maxpat_min = _maxpat_min;
    maxpat_max = _maxpat_max;
    enc = _enc;
    xml = _xml;
    where = _where;
    directed = _directed;
    report_result = _report_result; // For gSpanSubs
    verbose = _verbose;

    read (is, graphFilter);
    run_intern ();

}

std::vector<DFSCode> gSpan::run_gSpan_CORK (std::istream &is, std::ostream &_os,
        unsigned int _minsup,
        unsigned int _maxpat_min, unsigned int _maxpat_max,
        bool _enc, unsigned int _where,
        bool _xml, bool _directed,
        std::set<unsigned int> const * graphFilter,
        std::string _fs_option, std::string class_label_file,
        std::vector<std::map<unsigned int, unsigned int> > * supportMap,
        std::string subs_file, // For gSpanSubs
        bool _report_result,    // For gSpanSubs
        bool _verbose) {

    os = &_os;

    fs = false;
    ID = 0;
    minsup = _minsup;
    maxpat_min = _maxpat_min;
    maxpat_max = _maxpat_max;
    enc = _enc;
    xml = _xml;
    where = _where;
    directed = _directed;
    verbose = _verbose;
    delete_resolved = false;
    report_result = _report_result; // For gSpanSubs
    pattern_num = 0; // For gSpanSubs

    time_t time_now = time(NULL);
    start_time = clock();
    
    std::cerr << "reading graphs ......" <<std::endl;
    read (is, graphFilter); // neither is this

    std::cerr << "graph size: "<<TRANS.size()
              <<", time: "<< (float)(clock() - start_time) / CLOCKS_PER_SEC 
              << " s" << std::endl; // For gSpanSubs

    double minsupRatio = static_cast<double>(minsup) / TRANS.size();

    if (_fs_option != "") { // perform nested feature selection

        fs = true;
        fs_option = _fs_option;
        unsigned int selection_threshold = 0; // # subgraphs to select
        unsigned int correspondence_threshold = 0; // # correspondences allowed
        if ( (fs_option.size() > 4 && fs_option.substr(0,4) == "CORK") ||
                (fs_option.size() > 5 && fs_option.substr(0,5) == "FCORK") ) {
            if (verbose)
                std::cout<<"reading in class labels"<<std::endl;
            if (fs_option.size() > 5 && fs_option.substr(0,5) == "FCORK")
                delete_resolved = true;
            selection_threshold = strtoul(fs_option.substr((delete_resolved ? 5 : 4)).c_str(), NULL, 10);
            if (fs_option.size() > 4 && fs_option.find('C', 4) != std::string::npos && fs_option.find('C', 4) != fs_option.size())
                correspondence_threshold = strtoul(fs_option.substr(fs_option.find('C', 4)+1).c_str(), NULL, 10);
        } else if (fs_option.size() < 4 ||
                   ( (fs_option.size() == 4 && fs_option.substr(0,4) != "CORK") ||
                     (fs_option.size() == 5 && fs_option.substr(0,5) != "FCORK") ) ) {
            std::cerr << "error: do not know option \""<<fs_option<<"\" (gSpan::run_gSpan_CORK:"<< __LINE__<<")"<<std::endl;
            exit(1);
        }
        if (selection_threshold == 0) {
            selection_threshold = std::numeric_limits<unsigned int>::max();
            *infoStream << "warning: no upper limit for the maximum number of selected subgraphs given;"<<std::endl;
            *infoStream << "         using " << selection_threshold << " (gSpan::run_gSpan_CORK:" << __LINE__ << ")" << std::endl;
        }

        if (verbose)
            std::cout<<"starting with nested feature selection"<<std::endl;

        std::vector<unsigned int> class_counts = assign_class_labels(class_label_file, graphFilter);

        if (verbose)
            std::cout<<"labels assigned"<<std::endl;

        equ_classes = std::vector<unsigned int>(TRANS.size());
        possible_changes = std::vector<bool>(TRANS.size());
        correspondence_classes.resize(1); // 1st correspondence
        correspondence_classes[0].father = -1; // no original corresponding
        correspondence_classes[0].counts = class_counts;
        correspondences = oneAgainstRestCorresondences(class_counts);

        selectedSubgraphs.clear();

        *infoStream<<"\tcorrespondences\tcorrespondence_classes\tnum_unresolved\tsize_unresolved\tmax_unresolved\titeration\ttested\tfrequent\tminimal\twinners\tpruned\ttime_[s]"<<std::endl;

        time_now = time(NULL);

        // perform greedy forward selection: one feature == one pruned gSpan DFS code tree traversal
        for (unsigned int i = 0; i < selection_threshold; i++) {

            tested_subgraphs = 0;
            frequent_subgraphs = 0;
            minimal_subgraphs = 0;
            winner_subgraphs = 0;
            prunedFS_subgraphs = 0;

            current_best_subgraph = DFSCode(); // initialize with the empty subgraph

            if (verbose)
                std::cout<<"iteration "<<i<<std::endl;

            // run gSpan
            run_intern ();

            if (current_best_subgraph.size() == 0) { // no further best graph found
                *infoStream<<"note: terminating early: "<<i<<" instead of "<<selection_threshold<<" frequent subgraphs selected"<<std::endl;
                break;
            }

            // report selected subgraph
            ID = i;
            // if wanted: subgraphs can be reported online here
            if (report_result) // For gSpanSubs
                report(current_best_subgraph, current_best_support_counts);
            pattern_num ++;
            // determine the number of unresolved correspondences
            double numUnres = 0, corrsUnres = 0;
            unsigned int unresMax = 0;
            unsigned int numNotNull;
            for (std::vector<Correspondence>::iterator cor_iter = correspondence_classes.begin(); cor_iter != correspondence_classes.end(); cor_iter++) {
                unsigned int currentCorrSize = 0;
                numNotNull = 0;
                for (std::vector<unsigned int>::const_iterator cIt = cor_iter->counts.begin();
                        cIt != cor_iter->counts.end(); cIt++) {
                    currentCorrSize += *cIt;
                    if (*cIt != 0)
                        numNotNull++;
                }
                if (numNotNull > 1) { // at least two more class conflicts remain to be resolved
                    corrsUnres += currentCorrSize;
                    numUnres++;
                    if (unresMax < currentCorrSize)
                        unresMax = currentCorrSize;
                }
            }
            corrsUnres /= numUnres;
            *infoStream<<"\t"<<correspondences<<"\t"<<correspondence_classes.size()<<"\t"<<numUnres<<"\t"<<corrsUnres<<"\t"<<unresMax<<"\t"<< i<<"\t"<<tested_subgraphs<<"\t"<<frequent_subgraphs<<"\t"<<minimal_subgraphs<<"\t"<<winner_subgraphs<<"\t"<<prunedFS_subgraphs<<"\t"<<difftime(time(NULL), time_now)<<std::endl;

            // collect output graphs
            selectedSubgraphs.push_back(current_best_subgraph);
            if (supportMap != NULL) // and the traceback
                supportMap->push_back(translateSupportCounts(current_best_support_counts));

            time_now = time(NULL);

            // reset variables
            DFS_CODE.clear();

            // ensure the extension of the correspondence classes by the currently
            // selected subgraph
            getAskCORK(NULL, NULL, true); // call for extension only
            extendCORK(); // extend

            if (delete_resolved) {
                // Adapting minsup to the number of remaining graphs is not helpful for
                // small graph collections and ought to be used with care!
                minsup = static_cast<unsigned int>(floor(minsupRatio*TRANS.size()+.5));
                if (minsup == 0) // should not happen
                    minsup = 1;
            }

            if (correspondence_threshold > correspondences)
                break;
        }

        // collect final reporting output
        double numUnres = 0, corrsUnres = 0;
        unsigned int unresMax = 0, numNotNull;
        for (std::vector<Correspondence>::iterator cor_iter = correspondence_classes.begin(); cor_iter != correspondence_classes.end(); cor_iter++) {
            unsigned int currentCorrSize = 0;
            numNotNull = 0;
            for (std::vector<unsigned int>::const_iterator cIt = cor_iter->counts.begin();
                    cIt != cor_iter->counts.end(); cIt++) {
                currentCorrSize += *cIt;
                if (*cIt != 0)
                    numNotNull++;
            }
            if (numNotNull > 1) {
                corrsUnres += currentCorrSize;
                numUnres++;
                if (unresMax < currentCorrSize)
                    unresMax = currentCorrSize;
            }
        }
        corrsUnres /= numUnres;
        *infoStream<<"\t"<<correspondences<<"\t"<<correspondence_classes.size()<<"\t"<<numUnres<<"\t"<<corrsUnres<<"\t"<<unresMax<<"\tlast\t"<<tested_subgraphs<<"\t"<<frequent_subgraphs<<"\t"<<minimal_subgraphs<<"\t"<<winner_subgraphs<<"\t"<<prunedFS_subgraphs<<"\t"<<difftime(time(NULL), time_now)<<std::endl;

    } else { // no feature selection
        // For gSpanSubs
        if (subs_file != "") {
            std::cerr << "reading substructure constraints......" << std::endl;
            start_time = clock();
            read_SUBS_GRAPH(subs_file);
            std::cerr << SUBS_SET.size()
                      << " substructure constraints, time: " 
                      << (float)(clock() -  start_time) / CLOCKS_PER_SEC 
                      << " s" << std::endl;
            std::cerr << "starting with base substructure ......" << std::endl;
            BASE_SUBS_GRAPH.write(std::cerr);
        }

        if (verbose)
            *infoStream << "running gSpan without feature selection" << std::endl;

        run_intern ();
    }
    return selectedSubgraphs;
}

std::vector<DFSCode> gSpan::run_gSpan_CORK (std::istream &is, std::ostream &_os,
        unsigned int _minsup,
        unsigned int _maxpat_min, unsigned int _maxpat_max,
        bool _enc, unsigned int _where,
        bool _xml, bool _directed,
        std::set<unsigned int> const * graphFilter,
        std::string _fs_option,
        std::map<unsigned int,unsigned int> const * _instances2classLabels,
        std::vector<std::map<unsigned int, unsigned int> > * supportMap,
        std::string subs_file, // For gSpanSubs
        bool _report_result,    // For gSpanSubs
        bool _verbose) {
    if (_instances2classLabels != NULL)
        instances2classLabels = std::map<unsigned int,unsigned int>(*_instances2classLabels);
    return run_gSpan_CORK (is, _os, _minsup, _maxpat_min, _maxpat_max, _enc, _where, _xml, _directed, graphFilter, _fs_option, "", supportMap, subs_file, _report_result, _verbose); // For gSpanSubs
}


void gSpan::run_intern (void) {
    /* In case 1 node subgraphs should also be mined for, do this as
     * preprocessing step.
     */
    if (maxpat_min <= 1) {
        /* Do single node handling, as the normal gspan DFS code based processing
         * cannot find subgraphs of size |subg|==1.  Hence, we find frequent node
         * labels explicitly.
         */
        for (unsigned int id = 0; id < TRANS.size(); ++id) {
            for (unsigned int nid = 0 ; nid < TRANS[id].size() ; ++nid) {
                if (singleVertex[id][TRANS[id][nid].label] == 0) {
                    // number of graphs it appears in
                    singleVertexLabel[TRANS[id][nid].label] += 1;
                }

                singleVertex[id][TRANS[id][nid].label] += 1;
            }
        }
        /* All minimum support node labels are frequent 'subgraphs'.
         * singleVertexLabel[nodelabel] gives the number of graphs it appears
         * in.
         *
         * 1/1.5-class case: All nodelabels that do not appear at all have a
         *    gain of zero, hence we do not need to consider them.
         *
         * 2-class case: Not appearing nodelabels are counted negatively.
         */
        for (std::map<unsigned int, unsigned int>::iterator it =
                    singleVertexLabel.begin () ; it != singleVertexLabel.end () ; ++it) {
            if ((*it).second < minsup)
                continue;

            unsigned int frequent_label = (*it).first;

            /* Found a frequent node label, report it.
             */
            Graph g(directed);
            g.resize (1);
            g[0].label = frequent_label;

            /* [graph_id] = count for current substructure
             */
            std::vector<unsigned int> counts (TRANS.size ());
            for (std::map<unsigned int, std::map<unsigned int, unsigned int> >::iterator it2 =
                        singleVertex.begin () ; it2 != singleVertex.end () ; ++it2) {
                counts[(*it2).first] = (*it2).second[frequent_label];
            }

            std::map<unsigned int, unsigned int> gycounts;
            for (unsigned int n = 0 ; n < counts.size () ; ++n)
                gycounts[n] = counts[n];

            if (report_result) // For gSpanSubs
                report_single (g, gycounts);
            pattern_num ++;

        }
    } // END maxpatmin <= 1

    if (SUBS_SET.size() != 0) 
    // For gSpanSubs
    {
        if(SUBS_SET.size() == 1 && BASE_SUBS_GRAPH.size() == 2
            && BASE_SUBS_GRAPH[0].label == BASE_SUBS_GRAPH[1].label)
        // constraint is one single edge, relabel it, then use gSpan
        {
            std::cerr << "symmetric one-edge substructural constraint" << std::endl;
            std::cerr << "relabelling vertex and edge lables ......" << std::endl; 
            relabel_time = clock();
            swapLabelPair swapedge(-1, -1), swapv1(-1, -1);

            /* swap edge label */
            if (BASE_SUBS_GRAPH[0].edge[0].elabel != 0)
            {
                swapedge = std::make_pair(0, BASE_SUBS_GRAPH[0].edge[0].elabel); 
                BASE_SUBS_GRAPH[0].edge[0].elabel = 0;
                BASE_SUBS_GRAPH[1].edge[0].elabel = 0;
            }

            /* swap vertex label */

                if (BASE_SUBS_GRAPH[0].label != 0)
                {
                    swapv1 = std::make_pair (0, BASE_SUBS_GRAPH[0].label);
                    BASE_SUBS_GRAPH[0].label = 0;
                    BASE_SUBS_GRAPH[1].label = 0;
                }
            relabel(swapedge, swapv1);
            relabel_time = clock() - relabel_time;
            if (swapedge.first != -1)
                std::cerr << "swap edge label: " 
                          << swapedge.first << " <-> "
                          << swapedge.second << std::endl;
            if (swapv1.first != -1)
                std::cerr << "swap vertex label: " 
                          << swapv1.first << " <-> "
                          << swapv1.second << std::endl;
            std::cerr << "swapping label time: "
                      << (float)relabel_time / CLOCKS_PER_SEC
                      <<" s" << std::endl;

            std::cerr << "running frequent pattern mining ......" <<std::endl;
            mining_time = clock();
            gspan_filter_dup_time = 0;
            gspan_filter_sup_time = 0;
            gspan_edge_grow_time = 0;
            /* start normal gSpan process, but only add smalllest edge*/
            EdgeList edges;
            Projected root;
            start_time = clock();
            for (unsigned int id = 0; id < TRANS.size(); ++id) 
            {
                Graph &g = TRANS[id];
                for (unsigned int from = 0; from < g.size() ; ++from) 
                    if (g[from].label == 0)
                        for (Vertex::edge_iterator it = g[from].edge.begin(); it != g[from].edge.end(); it ++) 
                            if (it->elabel == 0 && g[it->to].label == 0)
                                root.push (id, &(*it), 0);
            }
            gspan_edge_grow_time = clock() - start_time;
            DFS_CODE.push (0, 1, 0, 0, BASE_SUBS_GRAPH[1].label);
            project (root);
            DFS_CODE.pop ();

            mining_time = clock() - mining_time;
            std::cerr << pattern_num << " frequent patterns are found" << std::endl;
            std::cerr << "total mining time:" 
                      << (float)mining_time / CLOCKS_PER_SEC
                      <<" s" << std::endl;
            std::cerr << "filter duplicate candidates time:"
                      << (float)gspan_filter_dup_time / CLOCKS_PER_SEC
                      <<" s" << std::endl;
            std::cerr << "filter support threshold time:"
                      << (float)gspan_filter_sup_time / CLOCKS_PER_SEC
                      <<" s" << std::endl;
            std::cerr << "get edge extention candidates time:"
                      << (float)gspan_edge_grow_time / CLOCKS_PER_SEC
                      <<" s" << std::endl;
        }
        else
        {
            get_BASE_SUBS_MINDFS();
            mining_time = clock() - mining_time;
            std::cerr << pattern_num << " frequent patterns are found" << std::endl;
            std::cerr << "total mining time:" 
                      << (float)mining_time / CLOCKS_PER_SEC
                      <<" s" << std::endl;
            std::cerr << "filter duplicate candidates time:"
                      << (float)filter_dup_time / CLOCKS_PER_SEC
                      <<" s" << std::endl;
            std::cerr << "filter substructure constraints time:"
                      << (float)filter_subs_time / CLOCKS_PER_SEC
                      <<" s" << std::endl;
            std::cerr << "filter support threshold time:"
                      << (float)filter_sup_time / CLOCKS_PER_SEC
                      <<" s" << std::endl;
            std::cerr << "get edge extention candidates time:"
                      << (float)edge_grow_time / CLOCKS_PER_SEC
                      <<" s" << std::endl;
        }
        return;
    }

    std::cerr << "running frequent pattern mining ......" <<std::endl;
    EdgeList edges;
    Projected_map3 root;
    gspan_mining_time = clock(); // For gSpanSubs  
    gspan_filter_sup_time = 0; // For gSpanSubs
    gspan_edge_grow_time = 0; // For gSpanSubs
    gspan_filter_dup_time = 0; // For gSpanSubs
    
    start_time = clock();
    for (unsigned int id = 0; id < TRANS.size(); ++id) {
        Graph &g = TRANS[id];

        for (unsigned int from = 0; from < g.size() ; ++from) {
            if (get_forward_root (g, g[from], edges)) {
                for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it)
                    root[g[from].label][(*it)->elabel][g[(*it)->to].label].push (id, *it, 0);
            }
        }
    }

    gspan_edge_grow_time = gspan_edge_grow_time + clock() - start_time;
    for (Projected_iterator3 fromlabel = root.begin() ;
            fromlabel != root.end() ; ++fromlabel) {
        for (Projected_iterator2 elabel = fromlabel->second.begin() ;
                elabel != fromlabel->second.end() ; ++elabel)	{
            for (Projected_iterator1 tolabel = elabel->second.begin();
                    tolabel != elabel->second.end(); ++tolabel)	{
                /* Build the initial two-node graph.  It will be grown
                 *  recursively within the project.
                 */
                DFS_CODE.push (0, 1, fromlabel->first, elabel->first, tolabel->first);
                project (tolabel->second);
                DFS_CODE.pop ();
            }
        }
    }

    /* For gSpanSubs*/
    gspan_mining_time = clock() - gspan_mining_time;
    std::cerr << pattern_num << " frequent patterns are found" << std::endl;
    std::cerr << "total mining time:" 
              << (float)gspan_mining_time / CLOCKS_PER_SEC
              <<" s" << std::endl;
    std::cerr << "filter duplicate candidates time:"
              << (float)gspan_filter_dup_time / CLOCKS_PER_SEC
              <<" s" << std::endl;
    std::cerr << "filter support threshold time:"
              << (float)gspan_filter_sup_time / CLOCKS_PER_SEC
              <<" s" << std::endl;
    std::cerr << "get edge extention candidates time:"
              << (float)gspan_edge_grow_time / CLOCKS_PER_SEC
              <<" s" << std::endl;
}

void gSpan::read_SUBS_GRAPH( std::string subs_file )
/* For gSpanSubs: read substructure constraintf from file
 * choose the one with most vertices and most edges as the base
 * for gSpan to grow, put it in the front of SUBS_SET
 */
{
    SUBS_SET.clear();
    std::ifstream subsFile(subs_file.c_str());

    if (!subsFile)
    {
        std::cerr<< "error: cannot open file containing restricted substructure \""<< subs_file <<"\" (gSpan::read_SUBS_GRAPH:"<< __LINE__<<")"<<std::endl;
        exit(1);
    }

    unsigned int max_edge_num = 0;
    unsigned int max_vertex_num = 0;
    while (true)
    {
        Graph g(directed);
        std::string subsName;
        g.read (subsFile, &subsName);
        if (g.empty())
            break;
        SUBS_SET.push_back (g);
        SUBS_IGRAPH_SET[SUBS_SET.size() - 1] = g.get_IGraph();
        if (g.vertex_size() > max_vertex_num ||
                (g.vertex_size() == max_vertex_num &&
                 g.edge_size() > max_edge_num))
        {
            BASE_SUBS_GRAPH_ID = SUBS_SET.size() - 1;
            max_vertex_num = g.vertex_size();
            max_edge_num = g.edge_size();
        }
    }
    subsFile.close();
    /* put base subs graph in the first place of SUBS_SET */
    BASE_SUBS_GRAPH = SUBS_SET[BASE_SUBS_GRAPH_ID];
    IGraph swap_temp = SUBS_IGRAPH_SET[BASE_SUBS_GRAPH_ID];
    SUBS_SET[BASE_SUBS_GRAPH_ID] = SUBS_SET[0];
    SUBS_IGRAPH_SET[BASE_SUBS_GRAPH_ID] = SUBS_IGRAPH_SET[0];
    SUBS_SET[0] = BASE_SUBS_GRAPH;
    SUBS_IGRAPH_SET[0] = swap_temp;
    BASE_SUBS_GRAPH_ID = 0;
}

void gSpan::get_BASE_SUBS_MINDFS(void)
/* For gSpanSubs: get minDFS for subs, and all supported graphs
 * similar to run_intern
 */
{
    EdgeList edges;
    Projected_map3 subs_root, dataset_root;
    bool result = false;

    // For gSpanSubs: generate all roots based on subs
    for (unsigned int from = 0; from < BASE_SUBS_GRAPH.size() ; ++from)
        if (get_forward_root (BASE_SUBS_GRAPH, BASE_SUBS_GRAPH[from], edges))
            for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it)
            {
                subs_root[BASE_SUBS_GRAPH[from].label][(*it)->elabel][BASE_SUBS_GRAPH[(*it)->to].label].push (BASE_SUBS_GRAPH_ID, *it, 0);
            }

    // For gSpanSubs: push graphs to correspondance root
    for (unsigned int id = 0; id < TRANS.size(); ++id)
    {
        Graph &g = TRANS[id];
        for (unsigned int from = 0; from < g.size() ; ++from)
            if (get_forward_root (g, g[from], edges))
                for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it)
                    if (subs_root[g[from].label][(*it)->elabel][g[(*it)->to].label].size() != 0)
                        dataset_root[g[from].label][(*it)->elabel][g[(*it)->to].label].push (id, *it, 0);
    }

    for (Projected_iterator3 fromlabel = dataset_root.begin() ;
            fromlabel != dataset_root.end() ; ++fromlabel)
        for (Projected_iterator2 elabel = fromlabel->second.begin() ;
                elabel != fromlabel->second.end() ; ++elabel)
            for (Projected_iterator1 tolabel = elabel->second.begin();
                    tolabel != elabel->second.end(); ++tolabel)
            {
                /* Build the initial two-node graph.  It will be grown
                *  recursively within the project.
                */
                DFS_CODE.push (0, 1, fromlabel->first, elabel->first, tolabel->first);
                result = project_SUBS_MINDFS (
                             subs_root[fromlabel -> first][elabel -> first][tolabel -> first], tolabel->second);
                if (result)
                    return;
                DFS_CODE.pop ();
            }

    std::cerr << "Program Exit: support of substructure's ancestor is smaller than minimum support threshold: (gSpan::get_SUBS_MINDFS:"<< __LINE__ <<")"<<std::endl;
    exit(0);
}

bool gSpan::project_SUBS_MINDFS (Projected &subs_projected, Projected &dataset_projected)
// For gSpanSubs: recursively compute min dfs for substructure
{
    bool result = false;

    unsigned int sup = support (dataset_projected);
    if (sup < minsup)
    {
        if (DFS_CODE.edgeCount() == BASE_SUBS_GRAPH.edge_size())
        {
            /* For gSpanSubs: a DFS code for BASE_SUBS_GRAPH is found
             * it may not be minimum DFS code
             * but the support does not change
             * since sup < minsup, exit
             */
            std::cerr << "Program Exit: support of all substructure constraints is smaller than minimum support threshold, no pattern will be generated: (gSpan::project_SUBS_MINDFS:"<< __LINE__ <<")"<<std::endl;
            exit(0);
        }
        else
            return false;
    }

    if (verbose)
        std::cout<<"frequent ("<<sup<<")"<<std::endl;

    if (DFS_CODE.is_min() == false)
    {
        return false;
    }

    if (verbose)
        std::cout<<"minimal"<<std::endl;

    if (DFS_CODE.edgeCount() == BASE_SUBS_GRAPH.edge_size())
    {
        // For gSpanSubs: min DFS code is found
        SUBS_MINDFS = DFS_CODE;
        SUBS_APPEARANCE = dataset_projected;
        std::cerr << sup << " graphs support base substructure" 
                  << ", time: " << (float)(clock() - start_time) / CLOCKS_PER_SEC
                  << " s" << std::endl;
        if (SUBS_SET.size() > 1)
        {
            std::cerr << "filtering graphs by other substructures ......" << std::endl;
            start_time = clock();
            SUBS_APPEARANCE = ask_all_subs(SUBS_APPEARANCE); 
            sup = support(SUBS_APPEARANCE);
            if (sup < minsup)
            {
                std::cerr << "Program Exit: support of all substructure constraints is smaller than minimum support threshold, no pattern will be generated: (gSpan::project_SUBS_MINDFS:"<< __LINE__ <<")"<<std::endl;
                exit(0);
            }
            std::cerr << sup << " graphs support all substructures found" 
                      << ", time: " << (float)(clock() - start_time) / CLOCKS_PER_SEC
                      << " s" << std::endl;
        }

        std::cerr << "running frequent pattern mining ......" <<std::endl;
        mining_time = clock();
        DFS_CODE.back().cover_subs_num = 1;
        filter_dup_time = 0;
        insert_pattern(DFS_CODE.get_IGraph());
        // minDFSTree.insert(SUBS_MINDFS);
        // the strategy of using minimum DFS code eliminate duplicates
        if(DFS_CODE.cover_all_subs(SUBS_IGRAPH_SET))
        {
            if (report_result) 
                report (SUBS_APPEARANCE, sup); // output the subs itself as a pattern
            pattern_num ++;
        }
        project_subs_get_root(SUBS_APPEARANCE);
        return true;
    }
    /* In case we have a node upper bound
    * We do not check for equality as the DFS exploration may
    * still only add edges but not nodes.
    */
    if (maxpat_max > maxpat_min && DFS_CODE.nodeCount () > maxpat_max)
        return false;

    /* We just find a minDFS code, now need to get its extensions,
    * and find the graphs that include this extended DFS code.
    */
    const RMPath &rmpath = DFS_CODE.buildRMPath ();
    int minlabel = DFS_CODE[0].fromlabel;
    int maxtoc = DFS_CODE[rmpath[0]].to;

    Projected_map3 new_subs_fwd_root;
    Projected_map2 new_subs_bck_root;
    Projected_map3 new_dataset_fwd_root;
    Projected_map2 new_dataset_bck_root;
    EdgeList edges;

    // Enumerate all possible one edge extensions for subs
    for (unsigned int n = 0; n < subs_projected.size(); ++n)
    {
        // for all appearances in the graph
        unsigned int id = subs_projected[n].id;
        PDFS *cur = &subs_projected[n];
        History history (BASE_SUBS_GRAPH, cur);

        // backward
        for (int i = (int)rmpath.size()-1; i >= 1; --i)
        {
            Edge *e = get_backward (BASE_SUBS_GRAPH, history[rmpath[i]], history[rmpath[0]], history);
            if (e)
            {
                new_subs_bck_root[DFS_CODE[rmpath[i]].from][e->elabel].push (id, e, cur);
            }
        }

        // forward
        if (get_forward_pure (BASE_SUBS_GRAPH, history[rmpath[0]], minlabel, history, edges))
            for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
            {
                new_subs_fwd_root[maxtoc][(*it)->elabel][BASE_SUBS_GRAPH[(*it)->to].label].push (id, *it, cur);
            }

        // backtracked forward
        for (int i = 0; i < (int)rmpath.size(); ++i)
            if (get_forward_rmpath (BASE_SUBS_GRAPH, history[rmpath[i]], minlabel, history, edges))
                for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it)
                {
                    new_subs_fwd_root[DFS_CODE[rmpath[i]].from][(*it)->elabel][BASE_SUBS_GRAPH[(*it)->to].label].push (id, *it, cur);
                }
    }

    // put all fitting graphs
    for (unsigned int n = 0; n < dataset_projected.size(); ++n)
    {
 
        unsigned int id = dataset_projected[n].id;
        PDFS *cur = &dataset_projected[n];
        History history (TRANS[id], cur);

        // backward
        for (int i = (int)rmpath.size()-1; i >= 1; --i)
        {
            Edge *e = get_backward (TRANS[id], history[rmpath[i]], history[rmpath[0]], history);
            if (e &&
                    new_subs_bck_root[DFS_CODE[rmpath[i]].from][e->elabel].size() != 0 )
                // For gSpanSubs: push only if BASE_SUBS_GRAPH already generates this pattern
                new_dataset_bck_root[DFS_CODE[rmpath[i]].from][e->elabel].push (id, e, cur);
        }

        // pure forward
        if (get_forward_pure (TRANS[id], history[rmpath[0]], minlabel, history, edges))
            for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
                if (new_subs_fwd_root[maxtoc][(*it)->elabel][TRANS[id][(*it)->to].label].size() != 0)
                    // For gSpanSubs: push only if BASE_SUBS_GRAPH already generates this pattern
                    new_dataset_fwd_root[maxtoc][(*it)->elabel][TRANS[id][(*it)->to].label].push (id, *it, cur);

        // backtracked forward
        for (int i = 0; i < (int)rmpath.size(); ++i)
            if (get_forward_rmpath (TRANS[id], history[rmpath[i]], minlabel, history, edges))
                for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it)
                    if (new_subs_fwd_root[DFS_CODE[rmpath[i]].from][(*it)->elabel][TRANS[id][(*it)->to].label].size() != 0)
                        // For gSpanSubs: push only if BAES_SUBS_GRAPH already generates this pattern
                        new_dataset_fwd_root[DFS_CODE[rmpath[i]].from][(*it)->elabel][TRANS[id][(*it)->to].label].push (id, *it, cur);
    }

    /* Test all extended substructures. */
    // backward
    for (Projected_iterator2 to = new_dataset_bck_root.begin(); to != new_dataset_bck_root.end(); ++to)
    {
        for (Projected_iterator1 elabel = to->second.begin(); elabel != to->second.end(); ++elabel)
        {
            DFS_CODE.push (maxtoc, to->first, -1, elabel->first, -1);
            result = project_SUBS_MINDFS (new_subs_bck_root[to -> first][elabel -> first], elabel->second);
            if (result)
                return true;
            DFS_CODE.pop();
        }
    }

    // forward
    for (Projected_riterator3 from = new_dataset_fwd_root.rbegin() ;
            from != new_dataset_fwd_root.rend() ; ++from)
        for (Projected_iterator2 elabel = from->second.begin() ;
                elabel != from->second.end() ; ++elabel)
            for (Projected_iterator1 tolabel = elabel->second.begin();
                    tolabel != elabel->second.end(); ++tolabel)
            {
                DFS_CODE.push (from->first, maxtoc+1, -1, elabel->first, tolabel->first);
                result = project_SUBS_MINDFS (new_subs_fwd_root[from -> first][elabel->first][tolabel->first], tolabel->second);
                if (result)
                    return true;
                DFS_CODE.pop ();
            }
    return false;
}

void gSpan::project_subs(Projected &projected)
/* For gSpanSubs: the BFS-based pattern grow function
 * recursively grow pattern
 */
{
    start_time = clock();
    unsigned int sup = support (projected);
    filter_sup_time = filter_sup_time + clock() - start_time;
    bool filter_subs_result = false;
    if (sup < minsup)
    {
        return;
    }

    if (verbose)
        std::cout<<"frequent ("<<sup<<")"<<std::endl;

//        if (minDFSTree.insert(DFS_CODE.get_min_DFS()) == false)
    if (insert_pattern(DFS_CODE.get_IGraph()) == false)
    {
//            minDFSTree.show(0); // lei debug
        return; // prune current pattern because of duplicate
    }
    // std::cerr << "after insert a dfscode into tree true" << std::endl;
//        minDFSTree.show(0); // lei debug
    if (verbose)
        std::cout<<"not duplicate"<<std::endl;

    /* In case we have a node upper bound
    * We do not check for equality as the DFS exploration may
    * still only add edges but not nodes.
    */
    if (maxpat_max > maxpat_min && DFS_CODE.nodeCount () > maxpat_max)
        return;
    start_time = clock();
    filter_subs_result = SUBS_SET.size() == 1 ||
        (SUBS_SET.size() > 1 && 
         DFS_CODE.cover_all_subs(SUBS_IGRAPH_SET));
    filter_subs_time = filter_subs_time + clock() - start_time;
    if (filter_subs_result)
    // only if current pattern cover all subs, report it
    // otherwise, grow without report
    {
        if (report_result) // For gSpanSubs
            report (projected, sup); // Output the frequent substructure
        pattern_num ++;
    }

    /* We just find a minDFS code, now need to get its extensions,
    * and find the graphs that include this extended DFS code.
    */
    Projected_map3 fwd_root;
    Projected_map3 bck_root;
    EdgeList edges;
    unsigned int maxtoc = DFS_CODE.nodeCount();

    // Enumerate all possible one edge extensions for subs
    for (unsigned int n = 0; n < projected.size(); ++n)
    {
        // for all appearances in the graph
        unsigned int id = projected[n].id;
        PDFS *cur = &projected[n];
        History history (TRANS[id], cur);
        history.buildSubsInfo();

        // backward edges
        start_time = clock();
        if (get_backward_subs(TRANS[id], history.back(), SUBS_MINDFS.nodeCount(), history, edges))
            for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
                bck_root[history.get_vertex_rank((*it) -> from)][history.get_vertex_rank((*it)->to)][(*it)->elabel].push (id, *it, cur);
                // use from id, to id, elabel as index

        //forward edges
        if (get_forward_subs(TRANS[id], history.lastFwdEdge, SUBS_MINDFS.nodeCount(), history, edges))
            for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
                fwd_root[history.get_vertex_rank((*it) -> from)][(*it)->elabel][TRANS[id][(*it)->to].label].push (id, *it, cur);
                // use from id, elabel, tolabel as index
        edge_grow_time = edge_grow_time + clock() - start_time;
    }
    
    /* recursively grow edges from backward root */
    for (Projected_iterator3 from = bck_root.begin();
            from != bck_root.end() ; ++from)
        for (Projected_iterator2 to = from->second.begin() ;
                to != from->second.end() ; ++to)
            for (Projected_iterator1 elabel = to->second.begin();
                    elabel != to->second.end(); ++elabel)
            {
                DFS_CODE.push (from->first, to->first, -1, elabel->first, -1);
                DFS_CODE.back().cover_subs_num = DFS_CODE[DFS_CODE.size()-2].cover_subs_num;
                project_subs(elabel -> second);
                DFS_CODE.pop();
            }

    /* recursively grow edges from forward root */
    for (Projected_iterator3 from = fwd_root.begin();
            from != fwd_root.end() ; ++from)
        for (Projected_iterator2 elabel = from->second.begin() ;
                elabel != from->second.end() ; ++elabel)
            for (Projected_iterator1 tolabel = elabel->second.begin();
                    tolabel != elabel->second.end(); ++tolabel)
            {
                DFS_CODE.push (from->first, maxtoc, -1, elabel->first, tolabel->first);
                DFS_CODE.back().cover_subs_num = DFS_CODE[DFS_CODE.size()-2].cover_subs_num;
                project_subs(tolabel -> second);
                DFS_CODE.pop();
            }
    return;
}

void gSpan::project_subs_get_root(Projected &projected)
/* For gSpanSubs: the BFS-based pttern grow function
 * set the backward edge and forward roots
 */
{
    edge_grow_time = 0;
    filter_subs_time = 0;
    filter_sup_time = 0;
    Projected_map3 fwd_root;
    Projected_map3 bck_root;
    EdgeList edges;
    unsigned int maxtoc = DFS_CODE.nodeCount();

    // Enumerate all possible one edge extensions for subs
    for (unsigned int n = 0; n < projected.size(); ++n)
    {
        // for all appearances in the graph
        unsigned int id = projected[n].id;
        PDFS *cur = &projected[n];
        History history (TRANS[id], cur);
        history.buildSubsInfo();

        start_time = clock();
        for (unsigned int i = 0; i < maxtoc; i ++)
        {
            // backward edges within subs
            if (get_backward_root_subs(TRANS[id], TRANS[id][history.get_vertex(i)], edges, history))
                for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it)
                    bck_root[i][history.get_vertex_rank((*it)->to)][(*it)->elabel].push (id, *it, cur);
                    // use from id, to id, e label as index

            // forward edges with one vertex in subs
            if(get_forward_root_subs(TRANS[id], TRANS[id][history.get_vertex(i)], edges, history))
                for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it) 
                    fwd_root[i][(*it)->elabel][TRANS[id][(*it) -> to].label].push (id, *it, cur);
                    // use from id, elabel, tolabel as index
        }
        edge_grow_time = edge_grow_time + clock() - start_time;
    }

    // recursively grow edges from backward root
    for (Projected_iterator3 from = bck_root.begin();
            from != bck_root.end(); from ++)
        for (Projected_iterator2 to = from->second.begin();
                to != from->second.end(); to ++)
            for (Projected_iterator1 elabel = to->second.begin();
                    elabel != to->second.end(); elabel ++)
            {
                DFS_CODE.push (from -> first, to -> first, -1, elabel -> first, -1);
                DFS_CODE.back().cover_subs_num = 1;
                // lei debug
                project_subs(elabel -> second);
                DFS_CODE.pop();
            }

    // recursively grow edges from forward root
    for (Projected_iterator3 from = fwd_root.begin();
            from != fwd_root.end(); from ++)
        for (Projected_iterator2 elabel = from->second.begin();
                elabel != from->second.end(); elabel ++)
            for (Projected_iterator1 tolabel = elabel->second.begin();
                    tolabel != elabel->second.end(); tolabel ++)
            {
                DFS_CODE.push (from -> first, maxtoc, -1, elabel -> first, tolabel->first);
                DFS_CODE.back().cover_subs_num = 1;
                project_subs(tolabel -> second);
                DFS_CODE.pop();
            }
}

bool gSpan::insert_pattern(IGraph graph)
// For gSpanSubs:  insert a discovered pattern to discoverd pattern pool
{
    start_time = clock();
    igraph_bool_t result;
    std::vector <IGraph> freq_patterns = FREQ_PATTERNS[igraph_vcount(&graph.graph)][igraph_ecount(&graph.graph)];
    for (unsigned int i = 0; i < freq_patterns.size(); i ++)
    {
        igraph_isomorphic_vf2(
                &(freq_patterns[i].graph), 
                &(graph.graph), 
                &(freq_patterns[i].vertex_labels), 
                &(graph.vertex_labels), 
                &(freq_patterns[i].edge_labels), 
                &(graph.edge_labels), 
                &result, NULL, NULL);
        if(result)
            return false;
    }
    FREQ_PATTERNS[igraph_vcount(&graph.graph)][igraph_ecount(&graph.graph)].push_back(graph);
    filter_dup_time = clock() - start_time + filter_dup_time;
    return true;
}

Projected gSpan::ask_all_subs(Projected projected)
// For gSpanSubs: filter seed
{
    unsigned int oid = 0xffffffff;
    IGraph ig;
    igraph_bool_t result = true;
    Projected filtered_projected;
    for (Projected::iterator cur = projected.begin(); cur != projected.end(); cur ++) {
        if (cur->id != oid) {
            ig = TRANS[cur -> id].get_IGraph();
            for(std::map<unsigned int, IGraph>::iterator it = SUBS_IGRAPH_SET.begin(); it != SUBS_IGRAPH_SET.end(); it ++)
            {
                igraph_subisomorphic_vf2(
                    &ig.graph,
                    &(it->second.graph),
                    &(ig.vertex_labels),
                    &(it->second.vertex_labels),
                    &(ig.edge_labels),
                    &(it->second.edge_labels),
                    &result, NULL, NULL);
                if (!result)
                    break;
            }
            oid = cur->id;
        }
        if (result)
            filtered_projected.push_back(*cur);
    }
    return filtered_projected;
}

void gSpan::relabel(swapLabelPair &swapedge, swapLabelPair &swapv1)
// For gSpanSubs: relabel strategy
{
    if (swapedge.first != -1) // include edge swapping
    {
        if (swapv1.first == -1) // no vertex swap
        {
            for (unsigned int gi = 0; gi < TRANS.size(); gi ++)
                for (unsigned int vi = 0; vi <TRANS[gi].size(); vi ++)
                    for (Vertex::edge_iterator ei = TRANS[gi][vi].edge.begin(); 
                         ei != TRANS[gi][vi].edge.end(); ei ++)
                        if (ei->elabel == swapedge.first)
                             ei->elabel = swapedge.second;
                        else if (ei->elabel == swapedge.second)
                            ei->elabel = swapedge.first;
        }
        else //if (swapv2.first == -1 && swapv1.first != -1)
        // one vertex swap
        {
            for (unsigned int gi = 0; gi < TRANS.size(); gi ++)
                for (unsigned int vi = 0; vi <TRANS[gi].size(); vi ++)
                {
                    if (TRANS[gi][vi].label == swapv1.first)
                        TRANS[gi][vi].label = swapv1.second;
                    else if (TRANS[gi][vi].label == swapv1.second)
                        TRANS[gi][vi].label = swapv1.first;
                    for (Vertex::edge_iterator ei = TRANS[gi][vi].edge.begin(); 
                         ei != TRANS[gi][vi].edge.end(); ei ++)
                        if (ei->elabel == swapedge.first)
                             ei->elabel = swapedge.second;
                        else if (ei->elabel == swapedge.second)
                            ei->elabel = swapedge.first;
                }
        }
    }
    else
    // no edge swapping
    {
        if (swapv1.first != -1)
        // one vertex swapping
            for (unsigned int gi = 0; gi < TRANS.size(); gi ++)
                for (unsigned int vi = 0; vi <TRANS[gi].size(); vi ++)
                    if (TRANS[gi][vi].label == swapv1.first)
                        TRANS[gi][vi].label = swapv1.second;
                    else if (TRANS[gi][vi].label == swapv1.second)
                        TRANS[gi][vi].label = swapv1.first;
    }
}
}
