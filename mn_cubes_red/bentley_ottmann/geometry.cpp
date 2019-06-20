/**
 * \file geometry.cpp
 * \brief Implementation of the methods and functions declared in
 * geometry.hpp.
 * \author Pierre Cagne
 * \date 07/11/2011
 */

#include "geometry.hpp"

/**
 * \brief Reads private member x.
 * \return Abcissa of the point
 */
rat Point::get_abscissa() const {
  return x;
}

/**
 * \brief Reads private member y.
 * \return Ordinate of the point
 */
rat Point::get_ordinate() const {
  return y;
}

/** 
 * \brief Assigns coordinates to the point.
 * \param x,y Coordinates to assign
 */
void Point::assign(const rat& x, const rat& y) {
  this->x = x, this->y = y;
  return ;
}

/**
 * \brief Initializes a slope.
 * \param x0,y0 Left endpoint defining the slope
 * \param x1,y1 Right endpoint defining the slope
 */
void Slope::make(const rat& x0, const rat& y0, const rat& x1, const rat&y1) {
  if(x0 == x1)
    type = infty;
  else {
    type = rational;
    value = (y0 - y1) / (x0 - x1) ;
  }

  return;
}

/**
 * \brief Defines an order on the slopes.
 * \param s Slope to compare with
 * \return True if and only is current slope less than s
 *
 * Slope with type 'infty' never are less than any other slope. If
 * rational, the slope is less than an 'infty' one and follows the
 * order on rational numbers for the general case.
 */
bool Slope::operator<(const Slope& s) const {
  if(type == infty)
    return false;
  else 
    return (s.type == infty || value < s.value);
}

/**
 * \brief Reads private memeber left
 * \return Left endpoint
 */
Point Segment::get_left() const {
  return left;
}

/**
 * \brief Reads private member right
 * \return Right endpoint
 */
Point Segment::get_right() const {
  return right;
}

/**
 * \brief Search for intersections
 * \param seg Segment to check if it is crossing the current one.
 * \param inter If there is an intersection, used to report it.
 * \return True if and only if there is an intersection.
 *
 * The searching is done by computing the determinant of the classic
 * system induced by the segment lines.
 */
bool Segment::intersect(const Segment& seg, Point& inter) const {
  // names
  Point s_left = seg.get_left(); Point s_right = seg.get_right();
  
  rat xa = left.get_abscissa(), ya = left.get_ordinate(),
    xb = right.get_abscissa(), yb = right.get_ordinate(),
    xc = s_left.get_abscissa(), yc = s_left.get_ordinate(),
    xd = s_right.get_abscissa(), yd = s_right.get_ordinate();

  static rat zero(I(0)), one(I(1));

  // determinant
  rat delta = (xb - xa)*(yc - yd) - (yb - ya)*(xc - xd);

  // Delta = 0 => no instersection
  if(delta == zero)
    return false;

  // computing r,s : Cramer's solutions
  rat r = (xc - xa)*(yc - yd) - (yc - ya)*(xc - xd); r /= delta;
  rat s = (yc - ya)*(xb - xa) - (xc - xa)*(yb - ya); s /= delta;
  
  // r, s in [0,1] else no intersection
  if(r<zero || r>one || s<zero || s>one)
    return false;

  rat xi = xa + r*(xb - xa), yi = ya + r*(yb - ya);
  inter.assign(xi, yi);
  
  return true;
}

/**
 * \brief Computes the high of a segment line.
 * \param p_sweep Sweeping point to use in the comparison.
 * \return High of the segment line
 *
 * We define the high of a segment line regarding to a point as :
 * <UL> <LI> not defined if the segment line does not cross the
 * vertical line with abscissa the abscissa of the point. </LI> <LI>
 * if the segment line is vertical <UL> <LI> the y-coordinate of the
 * point if the segment line contains the point. </LI> <LI> the
 * y-coordinate of the left endpoint if the point is laying below the
 * segment line. </LI> <LI> the y-coordinate of the right endpoint if
 * the point is laying above the segment line. </LI> </UL> <LI> the
 * y-coordinate of the crossing point of the segment and the vertical
 * line with abscissa the abscissa of the point else. </LI> </UL>
 */
rat Segment::high(const Point& p_sweep) const {
  static rat py, ly, ry;
  if(slope.type == infty) {
    py = p_sweep.get_ordinate(), ly = left.get_ordinate(),
      ry = right.get_ordinate();
    if(py < ly)
      return ly;
    else if(py > ry)
      return ry;
    else
      return py;
  }
  else {
    static rat xa, xb, ya, yb;
    xa = left.get_abscissa(), ya = left.get_ordinate(),
      xb = right.get_abscissa(), yb = right.get_ordinate();
    
    return (yb - ya)/(xb - xa)*(p_sweep.get_abscissa() - xa) + ya;
  }
}

/**
 * \brief Order of the segment relatively to a point.
 * \param s Segment to compare with.
 * \param p_sweep Sweeping point.
 * \return True if and only if the current segment is less than s
 *
 * A segment is less than an other if its high is less than the
 * other's high or if its high is equal to the other's high and its
 * slope is less than the other's slope.
 */
bool Segment::less(const Segment& s, const Point& p_sweep) const {
   return (
	   this->high(p_sweep) < s.high(p_sweep)
	   || (
	       this->high(p_sweep) == s.high(p_sweep)
	       && this->slope < s.slope
	       )
	   );
}

/**
 * \brief Checks if a point is the right endpoint of the segment line.
 * \param p Point to test.
 * \return True if and only if p is the right endpoint.
 */
bool Segment::is_rend(const Point& p) const {
  return (p.get_abscissa() == right.get_abscissa()
	  && p.get_ordinate() == right.get_ordinate());
}

/**
 * \brief Checks if a point is the left endpoint of the segment line.
 * \param p Point to test.
 * \return True if and only if p is the left endpoint.
 */
bool Segment::is_lend(const Point& p) const {
  return (p.get_abscissa() == left.get_abscissa()
	  && p.get_ordinate() == left.get_ordinate());
}

/**
 * \brief Checks if a point is in the interior of the segment line.
 * \param p Point to test.
 * \return True if and only if p is in the interior.
 */
bool Segment::is_in(const Point& p) const {
  //assuming p is on the segment
  return ((! is_rend(p)) && (! is_lend(p)));
}

// /**
//  * \brief Pretty-printing for mpz_class
//  * \param begin Ostream
//  * \param z Integer to print
//  * \return Ostream begin concatenated with the printing of z 
//  */
// std::ostream& operator<<(std::ostream& begin, const mpz_class& z) {
//   begin << z.get_str();
//   return begin;
// }

// /**
//  * \brief Pretty-printing for mpq_class
//  * \param begin Ostream
//  * \param f Rational to print
//  * \return Ostream begin concatenated with the printing of f 
//  */
// std::ostream& operator<<(std::ostream& begin, const mpq_class& f) {
//   begin << f.get_num() << "/" << f.get_den();
//   return begin;
// }

/**
 * Pointer version of the segment's \ref pp_seg "pretty-printing" 
 */
std::ostream& operator<<(std::ostream& begin, const Segment*& s) {
  begin << *s;
  return begin;
}

/**
 * \brief Pretty-printing for Point
 * \param begin Ostream
 * \param p Point to print
 * \return Ostream begin concatenated with the pretty-printing of p 
 */
std::ostream& operator<<(std::ostream& begin, const Point& p) {
  begin << "(" << p.get_abscissa()
	<< "," << p.get_ordinate()
	<< ")";
  return begin;
}

/**
 * \brief Pretty-printing for Event
 * \param begin Ostream
 * \param e Event to print
 * \return Ostream begin concatenated with the pretty-printing of e 
 */
std::ostream& operator<<(std::ostream& begin, const Event& e) {
  begin << e.get_point() ;
  return begin;
}

/**
 * \brief Pretty-printing for Segment \anchor pp_seg
 * \param begin Ostream
 * \param s Segment to print
 * \return Ostream begin concatenated with the pretty-printing of s 
 */
std::ostream& operator<<(std::ostream& begin, const Segment& s) {
  begin << "[(" << s.get_left().get_abscissa()
	<< "," << s.get_left().get_ordinate()
	<< "),(" << s.get_right().get_abscissa()
	<< "," << s.get_right().get_ordinate()
	<< ")]";
  return begin;
}

/**
 * \brief Order on events.
 * \param e Event to compare with.
 * \return True if and only if the current event is less than e
 *
 * The order is induced by the lexicographic order on the
 * (x,y)-coordinate of points.
 */
bool Event::operator<(const Event& e) const {
  return (
	  (point.get_abscissa() < e.point.get_abscissa())
	  || (
	      (point.get_abscissa() == e.point.get_abscissa())
	      && (point.get_ordinate() < e.point.get_ordinate())
	      )
	  );
}

/**
 * \brief Reads private member point.
 * \return Point of the event.
 */
Point Event::get_point() const {
  return point;
}
