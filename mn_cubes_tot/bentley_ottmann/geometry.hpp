/**
 * \file geometry.hpp
 * \brief Definition of geometric objects
 * \author Pierre Cagne
 * \date 07/11/2011
 *
 * Here are implemented classes for geometric objects. First, a simple
 * class for 2D points. Next, a class for 2D segments with possibility
 * of ordering. Last, a class for 'events', an ameloriated type of
 * point used in Bentley-Ottmann algorithm.
 *
 * \sa geometry.cpp 
 */

#ifndef H_GEO
#define H_GEO

#include <vector>
#include <iostream>
#include <gmpxx.h>
#include <boost/rational.hpp>
#include "integer_type.hpp"

/** 
 * \brief Rational type.
 *
 * This type can be a rational type with arbitrary precision (if
 * USE_BIGINT is true) using the GMP library. Otherwise, it just is a
 * rational type with numerator and denominator of type int.
 *
 * Another library can be use instead of GMP by changing only this two
 * typedef (and of course the included header(s) of the library). The
 * only requirement for a good behaviour are : I should have a
 * constructor from the paramater type int.
 */

#if USE_BIGINT

typedef mpz_class I;
typedef mpq_class rat;

#else 

typedef int I;
typedef boost::rational<int> rat;

#endif

/**
 * \brief 2D points with rational coordinates
 */

class Point {
private:
  rat x; //!< abscissa
  rat y; //!< ordinate

public:
 
  /**
   * \brief Default constructor
   */
  Point() : x(), y() {};

  /**
   * \brief Copy constructor
   * \param p Point to copy
   */
  Point(const Point& p) : x(p.x), y(p.y) {};
  
  /**
   * \brief Initializer constructor
   * \param xx,yy Abscissa and ordinate to affect to the point
   */ 
  Point(const rat& xx, const rat& yy) : x(xx), y(yy) {};
  rat get_abscissa() const;
  rat get_ordinate() const;
  void assign(const rat&,const rat&);
};


/**
 * \brief Kind of slope
 *
 * Defines a kind of slope for the slope structure. undef is only used
 * in default construction.
 */
enum slope_t {undef, rational, infty};

/**
 * \brief Slope of a segment line.
 */
struct Slope {
  slope_t type; //!< Kind of slope
  rat value; //!< if type = rational, value of that rational slope
  void make(const rat&,const rat&,const rat&,const rat&);
  bool operator<(const Slope&) const;
};

/**
 * \brief 2D segment lines
 *
 * Using class Point for the endpoints of the segment line, and struct
 * Slope to indicate the direction of the segment line.
 *
 * For a couple of point, one is said to be the left one and the other
 * the right one regarding to the lexographic order on the
 * (x,y)-coordinates.
 */
class Segment {
private:
  Point left; //!< Left endpoint
  Point right; //!< Right endpoint
  Slope slope; //!< Slope

public:

  /**
   * \brief Default constructor
   */
  Segment() : left(), right() { slope.type = undef; };

  /** 
   * \brief Copy constructor
   * \param seg Segment to copy
   */
  Segment(const Segment& seg) : left(seg.left), right(seg.right) {
    slope.type = seg.slope.type;
    if(slope.type == rational)
      slope.value = seg.slope.value;
  };

  /** 
   * \brief Assign constructor
   * \param a,b Points to be left and right ones of the segment line
   */
  Segment(const Point& a, const Point& b) : left(a), right(b) {
    slope.make(a.get_abscissa(), a.get_ordinate(), b.get_abscissa(), b.get_ordinate());
  };

  /**
   * \brief Assign constructor
   * \param xa,ya Coordinates of the left endpoint
   * \param xb,yb Coordinates of the right endpoint
   */
  Segment(const rat& xa, const rat& ya, const rat& xb, const rat& yb)
    : left(xa, ya), right (xb, yb) {
    slope.make(xa,ya,xb,yb);
  };
  Point get_left() const;
  Point get_right() const;
  rat high(const Point&) const;
  bool less(const Segment&, const Point&) const;
  bool intersect(const Segment&, Point&) const;
  bool is_rend(const Point&) const;
  bool is_lend(const Point&) const;
  bool is_in(const Point&) const;
};

// Pretty-printing
std::ostream& operator<<(std::ostream&, const Segment&);
std::ostream& operator<<(std::ostream&, const Segment*&);
std::ostream& operator<<(std::ostream&, const Point&);
//std::ostream& operator<<(std::ostream&, const mpz_class&);
//std::ostream& operator<<(std::ostream&, const mpq_class&);

/**
 * \brief Event for Bentley-Ottmann's algorithm
 *
 * Events are a certain kind of points. They are ordered (lex order)
 * and are remarkable segment points, that is left endpoints, right
 * endpoints and crossing points. So is one event refering to a set of
 * segment lines having it as remarkable point.
 */
class Event {
private:
  Point point; //!< Point
public:
  //! \brief Default constructor
  Event() : point() {};

  //! \brief Copy constructor \param e Event to copy
  Event(const Event& e) : point(e.point) {};

  //! \brief Assign constructor \param p Point to assign to the event
  Event(const Point& p) : point(p) {};
  Point get_point() const;
  bool operator<(const Event&) const;
};

#endif //H_GEO
