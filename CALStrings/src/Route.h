#ifndef ROUTE_H_
#define ROUTE_H_

#include <vector>
#include <string>

using namespace std;

struct nameAndCost { string name; unsigned int cost; };
bool compareByCost(const nameAndCost &a, const nameAndCost &b)
{
    return a.cost < b.cost;
}
class Route {
	string name;
	vector<string> stops;


public:
	Route(string name);
	void addStop(const string stop);
	bool operator==(const Route & r);
	void debug();
	int levDistance(const string source, const string target);
	nameAndCost approximateStringMatching(const string source);
	vector<string> getStops();
};

Route::Route(const string name): name(name) {}

void Route::addStop(const string stop) {stops.push_back(stop);}

bool Route::operator==(const Route & r) {return name == r.name;}

vector<string> Route::getStops(){return stops;}

int Route::levDistance(const string source, const string target)
	{
	  // Step 1
	  const int n = source.length();
	  const int m = target.length();
	  if (n == 0) {
	    return m;
	  }
	  if (m == 0) {
	    return n;
	  }
	  // Good form to declare a TYPEDEF
	  typedef vector<vector<int>> Tmatrix;
	  Tmatrix matrix(n+1);
	  // Size the vectors in the 2.nd dimension. Unfortunately C++ doesn't
	  // allow for allocation on declaration of 2.nd dimension of vec of vec
	  for (int i = 0; i <= n; i++) {
	    matrix[i].resize(m+1);
	  }
	  // Step 2
	  for (int i = 0; i <= n; i++) {
	    matrix[i][0]=i;
	  }
	  for (int j = 0; j <= m; j++) {
	    matrix[0][j]=j;
	  }
	  // Step 3
	  for (int i = 1; i <= n; i++) {
	    const char s_i = source[i-1];
	    // Step 4
	    for (int j = 1; j <= m; j++) {
	      const char t_j = target[j-1];
	      // Step 5
	      int cost;
	      if (s_i == t_j) {
	        cost = 0;
	      }
	      else {
	        cost = 1;
	      }
	      // Step 6
	      const int above = matrix[i-1][j];
	      const int left = matrix[i][j-1];
	      const int diag = matrix[i-1][j-1];
	      int cell = min( above + 1, min(left + 1, diag + cost));
	      // Step 6A: Cover transposition, in addition to deletion,
	      // insertion and substitution. This step is taken from:
	      // Berghel, Hal ; Roach, David : "An Extension of Ukkonen's
	      // Enhanced Dynamic Programming ASM Algorithm"
	      // (http://www.acm.org/~hlb/publications/asm/asm.html)
	      if (i>2 && j>2) {
	        int trans=matrix[i-2][j-2]+1;
	        if (source[i-2]!=t_j) trans++;
	        if (s_i!=target[j-2]) trans++;
	        if (cell>trans) cell=trans;
	      }
	      matrix[i][j]=cell;
	    }
	  }
	  // Step 7
	  return matrix[n][m];
	}

 	nameAndCost Route::approximateStringMatching(const string source)
	{
 	nameAndCost nc;
	const int maxCost = 20;
	int currCost=0;
	int cost = 100;
	//vector<string> stations = this->getAllStations();
	string finalTarget;
	string target;
	vector<string> stops = getStops();
	for(unsigned int i = 0; i < stops.size(); i++){
		target = stops[i];
		currCost = levDistance(source, target);
		if(currCost < maxCost){
			if(currCost < cost){
				finalTarget = target;
				cost = currCost;
			}
		}

	}
	nc.name = finalTarget;
	nc.cost = cost;
	return nc;
}

#endif /* ROUTE_H_ */
