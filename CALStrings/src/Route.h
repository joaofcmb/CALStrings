#ifndef ROUTE_H_
#define ROUTE_H_

#include <vector>
#include <string>

class Route {
	string name;
	string stops;

public:
	Route(string name);
	void addStop(string stop);
};

Route::Route(const string name): name(name) {stops = "";}

void Route::addStop(const string stop) {stops += "|" + stop + "|";}


#endif /* ROUTE_H_ */
