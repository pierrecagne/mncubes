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
  if(argc != 4) cout << "usage : ./mncubes [m] [n] [mu]\nExit 1.\n", exit(1);

  string mm(*++argv); string nn(*++argv); string mumu(*++argv);
  int n, m, num_mu, den_mu;
  (void) sscanf(nn.c_str(),"%d",&n);
  (void) sscanf(mm.c_str(),"%d",&m);
  (void) sscanf(mumu.c_str(),"%d/%d", &num_mu, &den_mu);
  rat mu(num_mu, den_mu);

  int p = num_mu, q = den_mu;

  cout << m << " , " << n << endl;

  vector<Segment> vec_seg;

  //general case
  int count_gen = 0;
  for(int j(1); j < n; j++) {
    int qj = q*j;
    for(int i(1); i < m; i++) {
      int qi = q*i;
      for(int k(1); q*k - p < qi + qj; k++) {
	int qk = q*k;
	if(gcd (i,gcd (j, qk-p)) <= q) {
	  // cout << i << "," << j << "," << k << endl;
	  rat x0, y0, x1, y1;
	  if(qk - p > qj) {
	    y0 = rat(1);
	    x0 = rat (qk-qj-p, qi);
	  }
	  else {
	    x0 = rat(0);
	    y0 = rat(qk-p, qj);
	  }
	  if(qk - p < qi) {
	    y1 = rat(0);
	    x1 = rat(qk-p , qi);
	  }
	  else {
	    x1 = rat(1); 
	    y1 = rat(qk-qi-p, qj);
	  }

#if USE_BIGINT
	  x0.canonicalize(), x1.canonicalize(), 
	    y0.canonicalize(), y1.canonicalize();
#endif

	  Segment seg (x0, y0, x1, y1);
	  count_gen++;
	  vec_seg.push_back (seg);
	}
      }
    }
  }
  
  //horizontal segments
  int count_horiz = 0;
  for(int j(1); j < n; j++) {
    int qj = q*j; 
    for(int k(1); q*k - p < qj; k++) {
      int qk = q*k;
      if(gcd (j, qk-p) <= q)
	count_horiz++, vec_seg.push_back (Segment (rat(0), rat(qk-p, qj), rat(1), rat(qk - p, qj)));
    }
  }
  vec_seg.push_back (Segment (rat(0), rat(0), rat(1), rat(0)));
  vec_seg.push_back (Segment (rat(0), rat(1), rat(1), rat(1)));

  //vectical segments
  int count_vert = 0;
  for(int i(1); i < m; i++) {
    int qi = q*i; 
    for(int k(1); q*k - p < qi; k++) {
      int qk = q*k;
      if(gcd (i, qk-p) <= q)
	count_vert++, vec_seg.push_back (Segment (rat(qk - p, qi), rat(0), rat(qk - p, qi), rat(1)));
    }
  }
  vec_seg.push_back (Segment (rat(0), rat(0), rat(0), rat(1)));
  vec_seg.push_back (Segment (rat(1), rat(0), rat(1), rat(1)));

  map_event vertices = bentley_ottmann(vec_seg);

  int sum_mul_vint = 0, nb_vint = 0;
  map_event::iterator it = vertices.begin();
  for(; it != vertices.end(); it++) {
    Point p = it->first.get_point();
    rat px = p.get_abscissa(), py = p.get_ordinate();
    if( px != rat(0) && px != rat(1) && py != rat(0) && py != rat(1) ) //p is in the inside
      sum_mul_vint += it->second.size(), nb_vint += 1;
  }

  int nb_lint = vec_seg.size() - 4;
  cout << "Nb general : " << count_gen << ", Nb horizontal : " << count_horiz
       << ", Nb vertical : " << count_vert << endl;
  cout << "vec_seg : " << vec_seg.size() << "\n";
  cout << "nb_lint : " << nb_lint << "\n";
  cout << "total_inter : " << vertices.size() << "\n";
  cout << "nb_vint : " << nb_vint << "\n";
  cout << "sum_mul_vint : " << sum_mul_vint << "\n";
  cout << "M_{" << m << "," << n << "}^(" << mu << ") = " << 1 + nb_lint + sum_mul_vint - nb_vint << "\n";
  
  return 0;
}
