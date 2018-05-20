#ifndef ROUTE_H_
#define ROUTE_H_

#include <vector>
#include <string>

#define MAX_ACCEPTED_COST	5

using namespace std;


class Stop {
	string name;
	int cost = 0;
public:
	int queueIndex = 0; 						// Required by MutablePriorityQueue

	Stop(string name);
	Stop(string name, int cost);
	Stop(const Stop &s);

	bool operator <(const Stop & stop) const; 	// Required by MutablePriorityQueue

	string getName() const;
	int getCost() const;

	void setCost(const int cost);
};

Stop::Stop(string name): name(name) {}

Stop::Stop(string name, int cost): name(name), cost(cost) {}
Stop::Stop(const Stop &s): name(s.getName()), cost(s.getCost()) {}

bool Stop::operator <(const Stop & stop) const {return (cost < stop.cost);}

string Stop::getName() const {return this->name;}

int Stop::getCost() const {return this->cost;}

void Stop::setCost(const int cost) {this->cost = cost;}

class Route {
	string name;
	vector<Stop> stops;


public:
	Route(string name);

	void addStop(const string stop);

	string getName() const;
	vector<Stop> getStops() const;

	bool hasStop(const string source);
	vector<Stop> possibleStops(const string source);
	int AproximateStringMatching(const string source, const string target);
	int StringMatching(const string text, const string pattern);
};

Route::Route(const string name): name(name) {}

void Route::addStop(const string stop) {stops.push_back(Stop(stop));}

string Route::getName() const {return name;}

vector<Stop> Route::getStops() const {return stops;}

int Route::AproximateStringMatching(const string text, const string pattern)
{
	// Step 1
	const int n = text.length();
	const int m = pattern.length();
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
		const char s_i = text[i-1];
		// Step 4
		for (int j = 1; j <= m; j++) {
			const char t_j = pattern[j-1];
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
				if (text[i-2]!=t_j) trans++;
				if (s_i!=pattern[j-2]) trans++;
				if (cell>trans) cell=trans;
			}
			matrix[i][j]=cell;
		}
	}
	// Step 7
	return matrix[n][m];
	/*
	int T = text.length();
	int P = pattern.length();

	if (T == 0)		return P;
	if (P == 0)		return T;
	if (T + P == 2)	return (pattern[0] == text[0]) ? 0 : 1;

	int dist[T], oldDist, newDist;

	for (int j = 0; j < T; j++)		dist[j] = j;

	for (int i = 1; i < P; i++) {
		oldDist = dist[0];
		dist[0] = i;
		for (int j = 1; j < T; j++) {
			newDist = (pattern[i] == text[j]) ? oldDist : 1 + min(oldDist, min(dist[j], dist[j-1]));
			oldDist = dist[j];
			dist[j] = newDist;
		}
	}

	return dist[T-1];
	*/
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

	return StringMatching(source, getName()) > 0;
}

vector<Stop> Route::possibleStops(const string source) {
	vector<Stop> v;
	int cost;

	for (auto stop : getStops()) {
		cost = AproximateStringMatching(source, stop.getName());

		if (cost <= MAX_ACCEPTED_COST)	{
			stop.setCost(cost);
			v.push_back(stop);
		}
	}

	cost = AproximateStringMatching(source, getName());
	if (cost <= MAX_ACCEPTED_COST)		v.push_back(Stop(getName(), cost));

	return v;
}

#endif /* ROUTE_H_ */
