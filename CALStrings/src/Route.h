#ifndef ROUTE_H_
#define ROUTE_H_

#include <vector>
#include <string>

using namespace std;

class Route {
	string name;
	string stops;

public:
	Route(string name);
	void addStop(string stop);

	bool operator==(const Route & r);
	void debug();
};

Route::Route(const string name): name(name) {stops = "";}

void Route::addStop(const string stop) {stops += "|" + stop + "|";}

bool Route::operator==(const Route & r) {return name == r.name;}

#endif /* ROUTE_H_ */
