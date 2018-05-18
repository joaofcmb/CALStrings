#ifndef ROUTE_H_
#define ROUTE_H_

#include <vector>
#include <string>

using namespace std;

class Route {
	string name;
	vector<string> stops;

public:
	Route(string name);
	void addStop(const string stop);

	bool operator==(const Route & r);
	void debug();
};

Route::Route(const string name): name(name) {}

void Route::addStop(const string stop) {stops.push_back(stop);}

bool Route::operator==(const Route & r) {return name == r.name;}

#endif /* ROUTE_H_ */
