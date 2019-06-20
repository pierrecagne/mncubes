#include <iostream>
#include <vector>
#include "bentley_ottmann/bentley_ottmann.hpp"
using namespace std;

int gcd(int a, int b) {
  if(a < b) return gcd (b, a);
  else if(b == 0) return a;
  else return gcd(b, a % b);
}

int main(int argc, char** argv) {
  if(argc != 3) cout << "usage : ./mncubes [m] [n]\nExit 1.\n", exit(1);

  string mm(*++argv); string nn(*++argv);
  int n,m;
  (void) sscanf(nn.c_str(),"%d",&n);
  (void) sscanf(mm.c_str(),"%d",&m);

  cout << m << " , " << n << endl;

  vector<Segment> vec_seg;

  //general case
  for(int j(1); j < n; j++)
    for(int i(1); i < m; i++)
      for(int k(1); k < i+j; k++) {
	if(gcd(j,gcd(i,k)) == 1) {
	  rat x0, y0, x1, y1;
	  if(k > j) {
	    y0 = rat(I(1));
	    x0 = rat (I(k - j), I(i));
	  }
	  else 
	    x0 = rat(I(0)), y0 = rat(I(k), I(j));
	  if(k < i) {
	    y1 = rat(I(0));
	    x1 = rat(I(k) , I(i));
	  }
	  else
	    x1 = rat(I(1)), y1 = rat(I(k - i), I(j));

#if USE_BIGINT
	  x0.canonicalize(), x1.canonicalize(), 
	    y0.canonicalize(), y1.canonicalize();
#endif

	  Segment seg (x0, y0, x1, y1);
	  vec_seg.push_back (seg);
	}
      }
  
  //horizontal segments
  for(int j(1); j < n; j++)
    for(int k(1); k < j; k++)
      if(gcd(j,k) == 1)
	vec_seg.push_back (Segment (rat(I(0)), rat(I(k), I(j)),rat(I(1)), rat(I(k),I(j))));
  vec_seg.push_back (Segment (rat(I(0)), rat(I(0)), rat(I(1)), rat(I(0))));
  vec_seg.push_back (Segment (rat(I(0)), rat(I(1)), rat(I(1)), rat(I(1))));

  //vectical segments
  for(int i(1); i < m; i++)
    for(int k(1); k < i; k++)
      if(gcd(i,k) == 1)
	vec_seg.push_back (Segment (rat(I(k), I(i)), rat(I(0)), rat(I(k), I(i)), rat(I(1))));
  vec_seg.push_back (Segment (rat(I(0)), rat(I(0)), rat(I(0)), rat(I(1))));
  vec_seg.push_back (Segment (rat(I(1)), rat(I(0)), rat(I(1)), rat(I(1))));

  map_event vertices = bentley_ottmann(vec_seg);

  int sum_mul_vint = 0, nb_vint = 0;
  map_event::iterator it = vertices.begin();
  for(; it != vertices.end(); it++) {
    Point p = it->first.get_point();
    rat px = p.get_abscissa(), py = p.get_ordinate();
    if( px != rat(I(0)) && px != rat(I(1)) && py != rat(I(0)) && py != rat(I(1)) ) //p is in the inside
      sum_mul_vint += it->second.size(), nb_vint += 1;
  }

  int nb_lint = vec_seg.size() - 4;
  cout << "vec_seg : " << vec_seg.size() << "\n";
  cout << "nb_lint : " << nb_lint << "\n";
  cout << "total_inter : " << vertices.size() << "\n";
  cout << "nb_vint : " << nb_vint << "\n";
  cout << "sum_mul_vint : " << sum_mul_vint << "\n";
  cout << "M_{m,n}^0 = " << 1 + nb_lint + sum_mul_vint - nb_vint << "\n";
  
  return 0;
}
