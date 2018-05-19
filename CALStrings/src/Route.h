#ifndef ROUTE_H_
#define ROUTE_H_

#include <vector>
#include <string>

using namespace std;


class Stop {
	string name;
	int cost = 0;

	int queueIndex = 0; 				// Required by MutablePriorityQueue

public:
	Stop(string name);
	bool operator<(Stop & stop) const; 	// Required by MutablePriorityQueue

	string getName() const;
};

Stop::Stop(string name): name(name) {}

bool Stop::operator<(Stop & stop) const {return (this->cost < stop.cost);}

string Stop::getName() const {return this->name;}


class Route {
	string name;
	vector<Stop> stops;

	int levDistance(const string source, const string target);
	int StringMatching(const string text, const string pattern);
public:
	Route(string name);

	void addStop(const string stop);

	string getName() const;
	vector<Stop> getStops() const;

	bool hasStop(const string source);
	//nameAndCost approximateStringMatching(const string source);
};

Route::Route(const string name): name(name) {}

void Route::addStop(const string stop) {stops.push_back(Stop(stop));}

string Route::getName() const {return name;}

vector<Stop> Route::getStops() const {return stops;}

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

int Route::StringMatching(const string text, const string pattern) {
	int T = text.length();
	int P = pattern.length();
	int matches = 0;

	// Compute Prefix Function
	int prefix[P];
	prefix[0] = 0;

	int k = 0;
	for (int q = 1; q < P; q++) {
		while (k > 0 && pattern[k] != pattern[q])
			k = prefix[k-1];
		if (pattern[k] == pattern[q])
			k++;
		prefix[q] = k;
	}

	// Search pattern in text
	int q = 0;
	for (int i = 0; i < T; i++) {
		while (q > 0 && pattern[q] != text[i])
			q = prefix[q-1];
		if (pattern[q] == text[i])
			q++;
		if (q == P) { // full match
			matches++;
			q = prefix[q-1];
		}
	}

	return matches;
}

bool Route::hasStop(const string source) {
	for (auto stop : getStops())
		if (StringMatching(source, stop.getName()) > 0)		return true;

	return false;
}

/*
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
*/
#endif /* ROUTE_H_ */
