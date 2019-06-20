/**
 * \file bentley_ottmann.hpp
 * \brief Some declarations use in bentley_ottmann.cpp
 * \author Pierre Cagne
 * \date 07/11/2011
 */

#ifndef H_BEN_OTT
#define H_BEN_OTT

#include <vector>
#include <map>
#include <utility>
#include "geometry.hpp"
#include "structures.hpp"

/**
 * \brief Map for events
 *
 * Map containing events associated to a vector of segment lines :
 * those which are crossing in this point.
 */
typedef std::map<Event, std::vector<Segment*> > map_event;

/**
 * \brief Element type of map_event.
 */
typedef std::pair<Event, std::vector<Segment*> > pair_event;

/**
 * \brief Vector of Segment*
 */
typedef std::vector<Segment*> vector_seg;

map_event bentley_ottmann(std::vector<Segment>&);
void find_neighboors(const Point&,BST<Segment,Point>::Type&,Segment*&,Segment*&);
Segment* find_leftmost(vector_seg&,const Point&);
Segment* find_rightmost(vector_seg&,const Point&);
Segment* find_left_neighboor(Segment*,BST<Segment,Point>::Type&);
Segment* find_right_neighboor(Segment*,BST<Segment,Point>::Type&);
std::pair<vector_seg, vector_seg> get_sets(const Point&,BST<Segment,Point>::Type&);
void compute_new_events(Segment*, Segment*,const Event&,PriorityQueue<Event,Segment*>&);


#endif //H_BEN_OTT
