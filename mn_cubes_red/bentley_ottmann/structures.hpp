/** 
 * \file structures.hpp 
 * \brief Declaration of structures needed in Bentley-Ottmann's
 * algorithm
 * \author Pierre Cagne 
 * \date 07/10/2001
 *
 * This file defines two structures used in Bentley-Ottmann's
 * algorithm. First, a binary search tree to maintain an order on some
 * objects. Next, a priority queue with associated elements and the
 * possibility of a quick research of a element. Those two structures
 * require the STL and assure that all actions are made in a
 * logarithmic way.
 */

#ifndef H_BST
#define H_BST

#include <set>
#include <vector>
#include <map>

/**
 * \class compare
 * \brief Implementation of an evolving order
 *
 * Bentley-Ottmann algorithm requires to make the order change with
 * time. This implementation will provide it by using a pointer on an
 * evolving paramater.
 *
 * \sa Segment
 */

template <class T, class Comp> struct compare {
  Comp* comp; /*!< Evolving parameter*/
  
  /**
   * \fn compare(Comp* c)
   * \brief Constructor
   *
   * \param c Initializer
   */
  compare(Comp* c) : comp(c) { }; 

  /**
   * \fn bool operator() (T* a, T* b)
   * \brief Define a comparator
   *
   * \param a,b Instances of T* to compare
   * \return True if and only if a < b regarding to the 'less' method
   */
  bool operator() (T* a, T* b) {
    return a->less(*b,*comp); 
  }
};

/**
 * \struct BST
 * \brief Red-black tree with evolving order
 *
 * Type 'set' of STL is implemented by a red-black tree (wich is a
 * binary search tree with amortized time complexity in the worst
 * case).
 *
 * The structure of tree is not harmed in any way by the changement on
 * the order. This assertion is pretty clear and purely algorithmic
 * and do not depend on the implementation.
 */

template <class T, class Comp> struct BST {
  typedef std::set <T*, compare<T,Comp> > Type;
};


/** 
 * \struct PriorityQueue_t 
 * \brief Basic structure for a priority queue
 *
 * This is the ground of our priority queue's definition. Type 'map'
 * of the STL uses a red-black tree assuring a logarithmic time
 * complexity. We use that kind of tree (binary search tree) instead
 * of the traditional heap to be able to search for a element in a
 * logarithmic way (heaps make it linear).
 */

template <class T, class A> struct PriorityQueue_t {
  typedef std::map<T, std::vector<A> > Type;
};


/**
 * \class PriorityQueue 
 * \brief Implementation of custom priority queues
 *
 * Inheritance assures the tree structure. We just define some other
 * methods to simplify the usage of the instances and keep the
 * vocabulary of the usual priority queues.
 */

template <class T, class A> class PriorityQueue : public PriorityQueue_t<T,A>::Type {
public:
  /**
   * \fn void push(const T& value)
   * \brief Push a new value with associated vector empty.
   *
   * \param value New value to push in the queue.
   */
  void push(const T& value) {
    std::vector<A> v; //empty vector
    this->insert(std::pair<T, std::vector<A> > (value, v));
  }

  /**
   * \fn void push(const T& value, const A& assoc)
   * \brief Push a new value in the queue or add an associated element
   * to an existing queue's member.
   *
   * \param value New value to push in the queue or existing member of
   * the queue.  
   * \param assoc Associated element to add to vector of value
   */
  void push(const T& value, const A& assoc) {
    std::vector<A> v; v.push_back(assoc);
    std::pair<typename PriorityQueue_t<T,A>::Type::iterator, bool> ret =
      this->insert(std::pair<T, std::vector<A> > (value, v));
    if(false == ret.second)
      ret.first->second.push_back(assoc);
    return;
  }

  /**
   * \fn T top()
   * \brief Extract the minimum of the tree.
   *
   * \return Top member of type T
   */
  T top() {
    return this->begin()->first;
  }

  /**
   * \fn void pop()
   * \brief Erase the top of the queue.
   */
  void pop() {
    this->erase(this->begin());
  }

  /**
   * \fn bool mem(const T& value)
   * \brief Test if a value is a member of the queue.
   *
   * \param value Value to test
   * \return True if and only if value is a member of the queue (independently to
   * the associated vector).
   */
  bool mem(const T& value) {
    return (this->find(value) != this->end());
  }
};

#endif //H_BST
