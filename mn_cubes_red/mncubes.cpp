#include <iostream>
#include <vector>
#include "bentley_ottmann/bentley_ottmann.hpp"
using namespace std;

int gcd(int a, int b) {
  if(a < b) return gcd (b, a);
  else if(b == 0) return a;
  else return gcd(b, a % b);
}

int totient(int n) {
  int res = 0;
  for(int i = 0; i < n; ++i)
    if(gcd(i,n) == 1)
      ++res;

  return res;
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
  int gen = 0;
  for(int j(1); j < n; j++)
    for(int i(1); i < m; i++)
      for(int k(1); 2*k < 2*i+j; k++) {
	if(gcd(j,gcd(i,k)) == 1) {
	  rat x0, y0, x1, y1;
	  if(2*k > j) {
	    y0 = rat(1, 2);
	    x0 = rat (2*k - j, 2*i);
	  }
	  else { 
	    x0 = rat(0);
	    y0 = rat(k, j);
	  }
	  if(k < i) {
	    y1 = rat(0);
	    x1 = rat(k , i);
	  }
	  else {
	    x1 = rat(1);
	    y1 = rat(k - i, j);
	  }

#if USE_BIGINT
	  x0.canonicalize(), x1.canonicalize(), 
	    y0.canonicalize(), y1.canonicalize();
#endif

	  Segment seg (x0, y0, x1, y1); gen++;
	  vec_seg.push_back (seg);
	}
      }
  
  //horizontal segments
  int hor = 0;
  for(int j = 2; j < n; j++)
    vec_seg.push_back (Segment (rat(0), rat(1, j), rat(1), rat(1, j))), hor++;
  vec_seg.push_back (Segment (rat(0), rat(0), rat(1), rat(0)));
  //vec_seg.push_back (Segment (rat(0), rat(1, 2), rat(1), rat(1, 2)));

  //vectical segments
  int vert = 0;
  for(int i(1); i < m; i++)
    for(int k(1); k < i; k++)
      if(gcd(i,k) == 1)
	vec_seg.push_back (Segment (rat(k, i), rat(0), rat(k, i), rat(1, 2))), vert++;
  vec_seg.push_back (Segment (rat(0), rat(0), rat(0), rat(1, 2)));
  vec_seg.push_back (Segment (rat(1), rat(0), rat(1), rat(1, 2)));

  map_event vertices = bentley_ottmann(vec_seg);

  int sum_mul_vint = 0, nb_vint = 0;
  map_event::iterator it = vertices.begin();
  for(; it != vertices.end(); it++) {
    Point p = it->first.get_point();
    rat px = p.get_abscissa(), py = p.get_ordinate();
    if(py.numerator() == 1) {
      if( px != rat(0) && px != rat(1) && py != rat(0) && py != rat(1)) { //p is in the inside
	int phi_d = totient(py.denominator());
	sum_mul_vint += phi_d*it->second.size(), nb_vint += phi_d*1;
      }
    }
  }

  int nb_lint = 0;
  for(int i = 1; i < m; i++)
    for(int j = 1; j < n; j++)
      for(int k = 1; k < i+j; k++)
	if(gcd(i, gcd(j, k)) == 1)
	  nb_lint++;
  for(int i = 1; i < m; i++)
    for(int k = 1; k < i; k++)
      if(gcd(i, k) == 1)
	nb_lint++;
  for(int j = 1; j < n; j++)
    for(int k = 1; k < j; k++)
      if(gcd(j, k) == 1)
	nb_lint++;

  cout << "Nb gen : " << gen << ", Nb hor : " << hor << ", Nb vert : " << vert << endl;
  cout << "vec_seg : " << vec_seg.size() << "\n";
  cout << "nb_lint : " << nb_lint << "\n";
  cout << "total_inter : " << vertices.size() << "\n";
  cout << "nb_vint : " << nb_vint << "\n";
  cout << "sum_mul_vint : " << sum_mul_vint << "\n";
  cout << "M_{m,n}^0 = " << 1 + nb_lint + sum_mul_vint - nb_vint << "\n";
  
  return 0;
}
