/**
 * \file bentley_ottmann.cpp
 * \brief Implementation of Bentley-Ottmann's algorithm
 * \author Pierre Cagne
 * \date 07/11/2011
 *
 * This file presents an implementation of Bentley-Ottmann algorithm
 * as shown \ref algo "here".
 */

/**
 * \mainpage
 * 
 * This project implements Bentley-Ottmann's algorithm. For a set of
 * segment lines, it computes and returns all the intersections of
 * segment lines of this set. In order to do that, we use the
 * generalized algorithm given in <em>M.de Berg, 0.Cheong, M.van
 * Kreveld and M.Overmars, Computational Geometry: Algorithms and
 * Applications, Third Edition</em>.
 *
 * This algorithm is based on the principle of <b>sweeping
 * line</b>. We <em>sweep</em> a line from the left to the right
 * (in some papers, it can be found from the top to the bottom) in the 2D
 * plan. At every position of the line, the intersections on the left
 * of the line has been computed and those on the right are waiting to
 * be treated. The interesting positions for the sweeping line are the
 * endpoints of the segment lines of the set and the crossing
 * points. We call those points <em>events</em>.
 *
 * We treat each event in a natural way : <ul> <li> if it is a left
 * endpoint of a segment line, this segment line can now potentially
 * be involved in an intersection. </li> <li> if it is a right
 * endpoint of a segment line, we don't have to care about this
 * segment line anymore : all the intersections that involves this
 * segment line are on the left of th sweeping line, so already
 * computed. </li> <li> if it is a crossing point of two or more
 * segment, we compute this intersection and move on. </li> </ul> (Of
 * course, a single event can be of more than one type.)
 *
 * In order to realize that, we have to maintain two structures : one
 * which gives the presently involved segment (those which are
 * crossing the sweeping line) called the Y-structure and one giving
 * the next event at every time called the X-structure. Since we want
 * to compute all the intersections with a time complexity \f$ O((n+k)
 * \log n) \f$ with \f$ n,k \f$ respectively the number of segments
 * and the number of intersections, we have to compute all the actions
 * on the X,Y-structures with a time complexity \f$ O(\log n)
 * \f$. Naturally, we choose to work with some trees. The X-structure
 * is a priority queue determining the next event at every time. For
 * some algorithmic reasons (we want to be able to know if an object
 * is a member of the X-structure in a logarithmic way), we will use a
 * binary search tree with the lexographic order on the
 * (x,y)-coordinates of the events. For the Y-structure, a binary
 * search tree is a natural choice to order the segment, insert and
 * remove them.
 *
 * \anchor algo
 * \image html "../../pic/bo_event.png" "An example of event with the three types of treatment."
 * \image latex "../../pic/bo_event.png" "An example of event (sweeping from the top)."
 *
 * At start, we push in the queue all the events corresponding to the
 * endpoints of the segment lines of the set. Then, for each event :
 * <ul> <li> we denote \f$ \mathcal L \f$ the set of segment having
 * the current event as left endpoints. (We suppose that this set is
 * stored with the event to keep our time complexity.) </li> <li>
 * Finding all segment lines in the binary search tree that contains
 * the event point, we partition those in \f$ \mathcal R \f$, which is
 * the set of the segment lines having the event point as right
 * endpoint, and \f$ \mathcal I \f$, which is the set of those
 * strictly containing the event point. </li> <li> If \f$ \mathcal L
 * \cup \mathcal I \cup \mathcal R \f$ contains strictly more than one
 * element, we can report one intersection. </li> <li> We now delete
 * from the binary search tree the segment lines in \f$ \mathcal I
 * \cup \mathcal R \f$ (those from \f$ \mathcal R \f$ because we does
 * not need them anymore ; those from \f$ \mathcal I \f$ to exchange
 * the order on them when re-inserting) </li> <li> It remains to treat
 * the possibly new events introduced by the segment lines from \f$
 * \mathcal L \cup \mathcal I\f$ : <ul> <li> if \f$ \mathcal L \cup
 * \mathcal I = \emptyset \f$ (i.e. the event point was just a right
 * endpoint of some segment lines, but was no left endpoint nor
 * crossing point), then we find s and s' the segment immediately
 * above and below the event point and check for any intersection
 * between them ; if there is one greather than the current event and
 * not already in the events's queue, then we push it. </li> <li> Else
 * we find the lower segment line of \f$ \mathcal L \cup \mathcal I
 * \f$ s and the upper one s' ; searching for the below neighboor r of
 * s and the above one r' of s', we check and insert if necessary any
 * new intersection between r and r'.</li> </ul> </ul> (Of course, the
 * several researches for neighboors depend of the existence of those
 * neighboors. If they don't, we just have nothing to do relatively to
 * it.)
 *
 * The files bentley_ottmann.hpp and bentley_ottmann.cpp implement
 * this algorithm. structures.hpp deals with the X,Y-structures using
 * the STL and some customizations. All of that depends of a
 * representation of some geometric objects, what is done in
 * geometry.hpp and geometry.cpp.
 */


#include "structures.hpp"
#include "geometry.hpp"
#include "bentley_ottmann.hpp"
#include <iostream>
#include <utility>
#include <iterator>
#include <algorithm>

/**
 * \brief Computes all the intersections between segment lines of set.
 * \param set Set of segment lines (implemented as Segment instances)
 * \return A map_event containing all the event crossing points
 *
 * Remark : the only Segment instances used are those in set ; to compute the
 * algorithm, we use pointers on Segment-s in order no to copy all
 * those data.
 */
map_event bentley_ottmann(std::vector<Segment>& set) {
  map_event inter;
  //std::cout << "Map : declared\n";
  // priority queue : insert all endpoints
  PriorityQueue<Event, Segment*> queue;

  std::vector<Segment>::iterator it;
  for(it = set.begin(); it != set.end(); it++) {
    Event e_left(it->get_left()), e_right(it->get_right());
    queue.push(e_left, &(*it)), queue.push(e_right);
  }

  // binary search tree : empty
  Point p_sweep(rat(I(0)),rat(I(0)));
  compare<Segment,Point> comp(&p_sweep);
  BST<Segment,Point>::Type btree(comp);

  // treatment of events
  while(queue.size()) {
    // current event : p event on top of set
    Event p = queue.top();
    Point point = p.get_point();
    vector_seg l_set = queue.begin()->second;    
    queue.pop();
    
    p_sweep.assign(point.get_abscissa(), point.get_ordinate());

    //checking for intersection
    std::pair<vector_seg, vector_seg> ir_set = get_sets(point, btree);
    
    if(l_set.size() + ir_set.first.size() + ir_set.second.size() > 1) {
      map_event::iterator it_ev = (inter.insert(pair_event(p, vector_seg()))).first;
      
      //it_ev->second : union of l_set, ir_set.first and ir_set.second
      vector_seg::iterator it_seg = l_set.begin();
      for(; it_seg != l_set.end(); it_ev->second.push_back(*it_seg++));
      for(it_seg = ir_set.first.begin(); it_seg != ir_set.first.end();
	  it_ev->second.push_back(*it_seg++));
      for(it_seg = ir_set.second.begin(); it_seg != ir_set.second.end();
	  it_ev->second.push_back(*it_seg++));
      //std::cout << it_ev->second << "\n";
    }
    

    //delete ir_set.first and ir_set.second from btree
    vector_seg::iterator it_seg = ir_set.first.begin();
    for(; it_seg != ir_set.first.end(); btree.erase(*it_seg++));
    for(it_seg = ir_set.second.begin(); it_seg != ir_set.second.end();
	btree.erase(*it_seg++));
        
    //update sl
    p_sweep.assign(point.get_abscissa(), point.get_ordinate());
    
    //insert l_set and ir_set.first in btree
    for(it_seg = l_set.begin(); it_seg != l_set.end();
	btree.insert(*it_seg++));
    for(it_seg = ir_set.first.begin(); it_seg != ir_set.first.end();
	btree.insert(*it_seg++));
    
    //treat new events
    if(l_set.size() + ir_set.first.size() == 0) {
      //p is only the rightpoint of several segments
      Segment* s_a, * s_b;
      find_neighboors(point, btree, s_a, s_b); 
      
      compute_new_events(s_a, s_b, p, queue);
    }
    
    else {
      vector_seg v(l_set.size() + ir_set.first.size());
      set_union(l_set.begin(), l_set.end(),
		ir_set.first.begin(), ir_set.first.end(), v.begin());
      Segment* sl = find_leftmost(v, point), * sr = find_rightmost(v, point);
      Segment* s_b = find_left_neighboor(sl, btree), * s_a = find_right_neighboor(sr, btree);
      
      compute_new_events(sl, s_b, p, queue); compute_new_events(sr, s_a, p, queue);
    }
  }
      
  return inter;
}

/**
 * \brief Find the two neighboors of a point
 * \param p Point that neighboors are searched for.
 * \param btree Binary search tree used to order segments.
 * \param above,below Two pointers use two store the potential
 * neighboors.
 *
 * This function create a segment of length 0, so a point, anchored in
 * p. In the tree, we search for the upper bound of the segment,
 * getting this way the upper neighboor (if exists). The lower
 * neighboor is now (if exists) the one just before the upper
 * neighboor. That's what this function strictely does.
 */
void find_neighboors(const Point& p, BST<Segment,Point>::Type& btree,
		     Segment*& above, Segment*& below) {
  //create a segment of length 0 representing p :
  rat x = p.get_abscissa(), y = p.get_ordinate();
  Segment s(x, y, x, y);
  
  //search for upper neighboor
  BST<Segment,Point>::Type::iterator it = btree.upper_bound(&s);
  
  if(it == btree.end())
    above = NULL;
  else
    above = *it;
  if(it == btree.begin())
    below = NULL;
  else 
    below = *--it;
}

/**
 * \brief Finds the lower segment.
 * \param v Non empty vector which for the lower segment is wanted.
 * \param p Point setting the order on segment
 * \return Pointer on the lower segment in the set.
 *
 * This is just a classical \f$ O(|v|) \f$ search for minimum in a vector v.
 */
Segment* find_leftmost(vector_seg& v, const Point& p) {
  //asuming v isn't empty
  vector_seg::iterator it = v.begin();
  Segment* min = *v.begin();

  while(++it != v.end()) {
    if((*it)->less(*min,p))
      min = *it;
  }
  return min;
}

/**
 * \brief Finds the upper segment.
 * \param v Non empty vector which for the upper segment is wanted.
 * \param p Point setting the order on segment
 * \return Pointer on the upper segment in the set.
 *
 * This is just a classical \f$ O(|v|) \f$ search for maximum in a vector v.
 */
Segment* find_rightmost(vector_seg& v, const Point& p) {
  //assuming v isn't empty
  vector_seg::iterator it = v.begin();
  Segment* max = *v.begin();
  
  while(++it != v.end()) {
    if(max->less(**it,p))
      max = *it;
  }
  return max;
}

/**
 * \brief Finds the lower neighboor of a segment in a BST.
 * \param s Pointer on the segment whose lower neighboor is wanted.
 * \param btree BST ordering the segment at the given time.
 * \return Pointer on the lower neighboor of s (if exists, else NULL).
 */
Segment* find_left_neighboor(Segment* s, BST<Segment,Point>::Type& btree) {
  BST<Segment,Point>::Type::iterator it = btree.find(s);
  if(it == btree.begin())
    return NULL;
  else
    return *--it;
}

/**
 * \brief Finds the upper neighboor of a segment in a BST.
 * \param s Pointer on the segment whose upper neighboor is wanted.
 * \param btree BST ordering the segment at the given time.
 * \return Pointer on the upper neighboor of s (if exists, else NULL).
 */
Segment* find_right_neighboor(Segment* s, BST<Segment,Point>::Type& btree) {
  BST<Segment,Point>::Type::iterator it = btree.find(s);
  if(++it == btree.end()) 
    return NULL;
  else
    return *it;
}

/**
 * \brief All segment lines containing a given point.
 * \param p Point contained in the wanted segment lines.
 * \param btree BST ordering the segment at the given time.
 * \return A pair of vector of segment lines : the first is \f$
 * \mathcal I \f$, the second \f$ \mathcal R \f$ defined \ref algo
 * "here".
 *
 * We use here the following trick : we first create a vertical
 * segment line of length 0 anchored in p ; this segment line have a
 * slope of value infinity and so is greather than any segment line
 * going through p except maybe some other vertical segment line
 * containing p ; that kind of vertical segment line is limited to a
 * single element because we do not accept overlapping segment
 * (undefined comportement of the algorithm) ; so having the upper
 * bound of the segment line with length 0 assure us to point on the
 * greatest segment line containing p ; it just remains to back
 * browsing the BST structure because all the segment containing p are
 * adjacent in the tree.
 */
std::pair<vector_seg, vector_seg> get_sets(const Point& p, BST<Segment,Point>::Type& btree) {
  vector_seg i,r;
  if(btree.empty())
    return std::pair<vector_seg, vector_seg> (i,r);
  //create a segment of lentgh zero representing p :
  rat x = p.get_abscissa(), y = p.get_ordinate();
  Segment* s = new Segment(x, y, x, y);
  
  
  BST<Segment,Point>::Type::iterator it = btree.upper_bound(s);
  std::reverse_iterator<BST<Segment,Point>::Type::iterator> rit(it);

  Point q;
  while(rit != btree.rend() && (*rit)->high(p) == y) {
    if((*rit)->is_in(p))
      i.push_back(*rit);
    else if((*rit)->is_rend(p))
      r.push_back(*rit);
    rit++;
  }

  return std::pair<vector_seg, vector_seg> (i,r);
}

/**
 * \brief Compute new events.  
 * \param s0,s1 Pointers on segment lines potentially introducing new
 * events.
 * \param p Current event : only events on the right of this one are
 * relevant.
 * \param q Priority queue which into the new events are pushed.
 *
 * First are checked s0 and s1 for equality to NULL. Next, lazily, s0
 * and s1 are checked for intersection. Next, lazily again, is checked
 * the intersection to be relevant or not. Finally, the event are
 * pushed into the queue.
 */
void compute_new_events(Segment* s0, Segment* s1, const Event& p,
			PriorityQueue<Event,Segment*>& q) {
  Point i;

  if(s0 && s1 && s0->intersect(*s1, i)) {
    Event ev_i(i);
    if(p < ev_i && ! q.mem(ev_i))
      q.push(ev_i);
  }
}
