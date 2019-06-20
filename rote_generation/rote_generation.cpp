/*! 
 * \file rote_generation.cpp
 * \brief Generation of the Rote sequences representing the m,n-cubes
 * of the discrete plans through (0,0,0).
 * \author Pierre Cagne
 * \date 08/04/2011
 */

/*! \mainpage
 *
 * <h2>Algorithm</h2>
 * Considering the graph presented in [DJVV08], one can fast generate
 * all the Rote sequences for m,n-cubes of plans through (0,0,0). Here
 * are presented the way to do it.
 *
 * First, each faces correpond to one and only one such m,n-cube. So,
 * to reach our goal, one have to visite each face at least once (and
 * for more efficiency only once).
 *
 * Secondly, knowing the Rote sequence of a given face, it is easy to
 * compute the Rote sequence of a neighboor face (i.e. a face sharing
 * an edge with the given one). It is a usefull property. Indeed,
 * randomly scanning the graph's faces, it is really not efficient to
 * track the corresponding m,n-cube <em>per se</em>. So, having this
 * property, the problem is to travel over all the faces, only once,
 * and through the edges : it is the search of a spanning
 * tree. Knowing that the lower most left face corresponds to the flat
 * m,n-cube (i.e. of the Rote sequence the matrix full of 0s), one can
 * perform a depth-first search (or a breadth-first one) to solve the
 * problem. That is indeed what it is done here.
 *
 * Last but not least, let's explain the computing of a new Rote
 * sequence going from a known face to a neighboor other through an
 * edge. Let \f$ f \f$ be the known face (i.e. the Rote sequence if
 * the m,n-cube representing by this face is known) and \f$ M \f$ the
 * matrix representation of the known Rote sequence. Let \f$ f' \f$
 * and \f$ M' \f$ be the homologous unknown ones. Finally, let \f$ e
 * \f$ be the sharing edge (our graph is such that only one edge can
 * be shared by two given bounded faces). Let also the equation of the
 * standing line of \f$ e \f$ be : \f$ ix+jy = k \f$. <br/> Then, by
 * construction of the graph, for a plan \f$ \mathcal P : z = -(\alpha
 * x + \beta y) \f$ whose dual point is in \f$ f \f$ (let's say \f$ f
 * \subset \{(x,y)\in \mathbb R^2 / ix+jy < k\} \f$), and a plan \f$
 * \mathcal P' : z = -(\alpha' x + \beta' y) \f$ whose dual point is
 * in \f$ f' \f$ (so is \f$ f' \subset \{(x,y)\in \mathbb R^2 / ix+jy
 * > k \} \f$), we have 
 *
 * \f[ \begin{aligned} \forall (p,q) \notin i\mathbb Z \times j\mathbb
 * Z,\ &-\lfloor \alpha p + \beta q \rfloor = -\lfloor \alpha' p +
 * \beta' q \rfloor\\ \forall r \geq 1,\ &-(\lfloor \alpha ri + \beta
 * rj \rfloor +1) = -rk = -\lfloor \alpha' ri + \beta' rj \rfloor
 * \end{aligned} \f]
 * 
 * So, translating this modulo 2, we have \anchor flip_algo
 *
 * \f[ \begin{aligned} \forall (p,q) \notin i\mathbb Z \times j\mathbb
 * Z,\ & M_{p,q} = M'_{p,q}\\ \forall r \geq 1,\ & M_{ri,rj} =
 * \overline{M'_{ri,rj}} \end{aligned} \f]
 *
 * where \f$ \overline{0} = 1 \f$ and \f$ \overline{1} = 0 \f$.
 *
 * And so, computing a Rote sequence from another going through an
 * edge with equation \f$ix+jy=k\f$ is done by flipping the values in
 * position \f$(ri,rj), r \geq 1, ri < m, rj < n\f$ in the known Rote
 * sequence.
 *
 * <h2>Implementation</h2>
 * 
 * In order to implement this algorithm, one have to maintain the
 * structure of the graph while inserting the segment lines. We will
 * here use <tt>CGAL</tt> (Computational Geometry Algorithms Library)
 * and its module <tt>Arrangement_2</tt>.<br/> <tt>CGAL</tt> requires
 * a kernel to do some geometric statements : we will work in the 2D
 * plan with rational (based on <tt>GMP</tt>) cartesian
 * coordinates. That the meaning of the three first
 * <tt>typedef</tt>. We then have to deal with the arrangements : the
 * initial structure of arrangement in <tt>CGAL</tt> is not helpfull
 * in the sense that the vertices, edgesand faces can not contains
 * data (other that the graph's structure). Hopefully, <tt>CGAL</tt>
 * offers a way to extend them called <tt>DCEL extended</tt>
 * (Double-Connected Edge List extended). Data in the vertices are
 * irrelevant for us. Data in the edges will be i and j from the
 * equation (see above). Data in the faces will be a Rote sequence (a
 * double entry array of char), and a boolean used in the depth-first
 * search.<br/> We will construct the graph by incremental insertions
 * (that avoids a preliminary stock of the segment lines), then set
 * the data on the edges and finally perform the depth-first search.
 *
 * The \ref dfs "depth-first search" is computed on
 * <tt>Face_handle</tt> (that is kind of a pointer on a face) and does
 * as follows : <ul><li> Create an empty list of face handlers to
 * return at the end.</li> <li> Create an empty stack of face
 * handlers, push the root face into it and mark th root face as
 * discovered.</li> <li> While the stack isn't empty : <ul> <li> Pop
 * the head, mark it as discovered, and push it into the list.</li>
 * <li> For every neighboor face of the head : if not discovered and
 * if bounded, mark as discovered, push it into the stack and \ref
 * flip "flip" its Rote sequence from the head's Rote sequence (as
 * shown \ref flip_algo "here").</li></ul></li></ul> The \ref dfs
 * "dfs" function will return a list cointaining the face handlers in
 * a depth-first order.
 *
 * Let's go : rote_generation.cpp
 */

#include <iostream>
#include <cstdlib>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <boost/math/common_factor_rt.hpp>
#include <stack>

using namespace std;

int m, n;

/*! Data on vertices for DCEL extended
 * 
 * DCEL extended needs a data structure for vertices, even a empty
 * one.
 */
struct VertexData {};

/*! Data on edges for DCEL extended
 *
 * The associated edge lies on the line with equation : ix+jy=k for
 * some k. (Specify k is irrelevant.)
 */
struct EdgeData {
  int i, j;
  EdgeData() : i(0), j(0) {};
  EdgeData(const int& i_, const int& j_) : i(i_), j(j_) {};
};

/*! Data on faces for DCEL extended */
struct FaceData {
  bool discovered; //!< Mark for depth-first search. \sa dfs
  char** matrix; //!< Two dimensional Rote sequence of the m,n-cubes
		 //!representing by the face.
  /*! Default constructor.
   *
   * 'discovered' is false, 'matrix' is allocated with only zeros.
   */
  FaceData() : discovered(false) {
    typedef char* char_ptr;
    matrix = new char_ptr [m];
    matrix[0] = new char[m*n];
    for(int i = 1; i < m; ++i)
      matrix[i] = matrix[i-1] + n;
    for(int i = 0; i < m; i++)
      for(int j = 0; j < n; ++j)
  	matrix[i][j] = 0;
  }
  void delete_face() {
    delete [] matrix[0];
    delete [] matrix;
  }
};

/*! Pretty-print for matrices.
 *
 * \pre matrix has to be of size \f$ m \times n \f$.
 * \return Updated stream.
 */
ostream& operator<< (ostream& stream, char** matrix) {
  // assuming the matrix has size m x n
  for(int i = 0; i < m; ++i) {
    for(int j = 0; j < n; ++j)
      stream << "[" << (int) matrix[i][j] << "]";
    stream << endl;
  }
  return stream;
}

/*! CGAL's arrangements require an exact integer type. */
typedef CGAL::Gmpz Integer;

/*! All intersections and endpoints are known to be rational. */
typedef CGAL::Gmpq rat;

/*! Compute a geometric kernel for CGAL's facilities. */
typedef CGAL::Cartesian<rat> Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel> Traits_2;
typedef Traits_2::Point_2 Point_2;
typedef Traits_2::X_monotone_curve_2 Segment_2;

/*! DCEL extended : graph with data on his components (vertices, edges, faces). */
typedef CGAL::Arr_extended_dcel<Traits_2,
				VertexData,
				EdgeData,
				FaceData> Dcel;

/*! Defining the arrangement type. */
typedef CGAL::Arrangement_2<Traits_2, Dcel> Arrangement_2;

/*! Handle a face. Facility for handling dual graphs.
 * 
 * \anchor Face_ptr 
 */
typedef Arrangement_2::Face_handle Face_ptr;


/*! Flip on m,n-cubes \anchor flip
 * 
 * \par i,j Parametrizing the flip.  
 * \par matrix_src Source matrix (Rote sequence of a m,n-cube) ; start
 * for the flip.
 * \par matrix_dest Initially full with 0s. Will contain the flip from
 * matrix_src.
 *
 * Going from a face to another through an edge with (standing line
 * of) equation \f$ ix+jy=k (i\wedge j\wedge k=1)\f$ will change from
 * 0 to 1 or from 1 to 0 every position \f$ (ri,rj) (r \geq 1)\f$ of
 * the Rote sequence representing the m,n-cube of the source face.
 */
void flip (int const& i, int const& j, char** matrix_src, char** matrix_dest) {
  // assuming i, j and k are in the correct interval, and matrices have
  // size m x n
  for(int r = 0; r < m; ++r)
    for(int s = 0; s < n; ++s)
      matrix_dest[r][s] = matrix_src[r][s];
  for (int r = int(1); r*i < int(m) && r*j < int(n); ++r)
    matrix_dest[r*i][r*j] = 1 - matrix_dest[r*i][r*j];
}

/*! Depth-first search on the dual graph \anchor dfs
 * 
 * \par graph Primal graph
 * \par f Handle for the root face
 * \return List of \ref Face_ptr "Face_ptr" performing a depth-first search (and so a
 * spanning tree).
 *
 * This will compute a depth-first search on the dual graph of
 * graph. This way we visit all the faces (so all the Rote sequence)
 * of the graph only once. Handling a face with some Rote sequence, we
 * can update all its still unvisited neighboor faces. This is the
 * core of the computation, everything else is either purely technical
 * or cosmetic.
 */
std::list <Face_ptr> dfs(Arrangement_2 const& graph, Face_ptr f) {
  std::list <Face_ptr> depth_first_result;
  std::stack <Face_ptr> vertices;

  vertices.push(f);
  f->data().discovered = true;
  while ( ! vertices.empty() ) {
    Face_ptr head = vertices.top();
    vertices.pop();
    
    depth_first_result.push_back(head);
    //    cout << head->data().matrix << endl;
    
    Arrangement_2::Ccb_halfedge_circulator first, curr, flip_op;
    first = curr = head->outer_ccb();
    do {
      Face_ptr son = curr->twin()->face();
      if( ! son->is_unbounded() && ! son->data().discovered) {
	vertices.push(son);
	son->data().discovered = true;
	flip(curr->data().i, curr->data().j, head->data().matrix, son->data().matrix);
      }
    } while(++curr != first);
  }
  
  return depth_first_result;
}

/*! Main function
 *
 * \function int main(int argc, char** argv) 
 * \param argc, argv Classic paramters.
 * \pre argc = 3
 * \pre argv contains 2 strings, each containing an <tt>int</tt>.
 * \pre Violating one of the conditions will exit with a error
 * message.
 *
 * Might be implement as a function to call in a larger program.
 * 
 * \todo Try to set data on edges while inserting segment lines. This
 * would avoid to cover them all after.
 */
int main(int argc, char** argv) {
  if(argc != 3) cout << "usage : ./generator [m] [n]\nExit 1.\n", exit(1);

  string mm(*++argv); string nn(*++argv);
  (void) sscanf(nn.c_str(),"%d",&n);
  (void) sscanf(mm.c_str(),"%d",&m);
  
  cout << " -- " << m << "," << n << "-cubes generation --" << endl;

  Arrangement_2 arr;

  //general case
  for(int j(1); j < n; j++)
    for(int i(1); i < m; i++)
      for(int k(1); k < i+j; k++) {
  	if(boost::math::gcd(j,boost::math::gcd(i,k)) == 1) {
  	  rat x0, y0, x1, y1;
  	  if(k > j) {
  	    y0 = rat(1);
  	    x0 = rat (k - j, i);
  	  }
  	  else if (k < j)
  	    x0 = rat(0), y0 = rat(k, j);
	  else
	    x0 = rat(0), y0 = rat(1);
  	  if(k < i) {
  	    y1 = rat(0);
  	    x1 = rat(k , i);
  	  }
  	  else if (k > i)
  	    x1 = rat(1), y1 = rat(k - i, j);
	  else 
	    x1 = rat(1), y1 = rat(0);
	  
	  insert(arr, Segment_2 (Point_2 (x0,y0), Point_2 (x1,y1)));
  	}
      }


  //horizontal segments
  for(int j(1); j < n; j++)
    for(int k(1); k < j; k++)
      if(boost::math::gcd(j,k) == 1)
  	insert(arr, Segment_2 (Point_2 (rat (0), rat (k, j)),
			       Point_2 (rat(1), rat(k, j))));

  insert(arr, Segment_2 (Point_2 (rat(0), rat(0)),
			 Point_2 (rat(1), rat(0))));
  
  insert(arr, Segment_2 (Point_2 (rat(0), rat(1)),
			 Point_2 (rat(1), rat(1)))); 


  //vectical segments
  for(int i(1); i < m; i++)
    for(int k(1); k < i; k++)
      if(boost::math::gcd(i,k) == 1)
  	insert(arr, Segment_2 (Point_2 (rat(k, i), rat(0)),
			       Point_2 (rat(k, i), rat(1))));

  insert(arr, Segment_2 (Point_2 (rat(0), rat(0)),
			 Point_2 (rat(0), rat(1))));

  insert(arr, Segment_2 (Point_2 (rat(1), rat(0)),
			 Point_2 (rat(1), rat(1))));

  Arrangement_2::Face_handle start; //flat m,n-cube
  Arrangement_2::Edge_iterator eit;
  

  for(eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
    Point_2 a = eit->source()->point(), b = eit->target()->point();

    if((a.x() == 0 && a.y() == 0) || (b.x() == 0 && b.y() == 0))
      if( ! eit->face()->is_unbounded())
	start = eit->face();
      else
	start = eit->twin()->face();

    if((a.x() - b.x()) != 0) { 
      rat alpha = (a.y() - b.y()) / (a.x() - b.x());
      rat beta = a.y() - alpha*a.x();

      Integer p = alpha.numerator(), q = alpha.denominator(),
  	p_ = beta.numerator(), q_ = beta.denominator();
      Integer gcd_ = boost::math::gcd(q*q_, boost::math::gcd(-p*q_, p_*q));
      int i = (int) (-p*q_/gcd_).to_double(), j = (int) (q*q_/gcd_).to_double();
      eit->set_data(EdgeData(i,j));
      eit->twin()->set_data(EdgeData(i,j));
    }
    else {
      rat k = a.x();
      Integer p = k.numerator(), q = k.denominator();
      Integer gcd_ = boost::math::gcd(p,q);
      int i = (int) (q/gcd_).to_double();
      eit->set_data(EdgeData(i, 0));
      eit->twin()->set_data(EdgeData(i, 0));
    }
  }
 
  std::list <Face_ptr> rote_sequences = dfs(arr, start);
  for(list<Face_ptr>::iterator it = rote_sequences.begin();
      it != rote_sequences.end(); it++) {
    cout << (*it)->data().matrix << endl;
    (*it)->data().delete_face();
  }
  cout << "Number of m,n-cubes of plans through (0,0,0) : " << rote_sequences.size() << endl;

  return 0;
}
