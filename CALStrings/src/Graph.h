/*
 * Graph.h
 */
#ifndef GRAPH_H_
#define GRAPH_H_

#include <vector>
#include <algorithm>
#include <queue>
#include <list>
#include <limits>
#include <cmath>
#include <cfloat>
#include <pthread.h>
#include "MutablePriorityQueue.h"
#include "Route.h"


using namespace std;

template <class T> class Edge;
template <class T> class Graph;
template <class T> class Vertex;

#define INF numeric_limits<double>::max()

/************************* Vertex  **************************/

template <class T>
class Vertex {
	T info;                // contents
	vector<Edge<T> *> adj;  // outgoing edges
	bool visited = false;          // auxiliary field

	// Improved algorithm uses dijkstra thrice, concurrently
	double dist = 0;
	double distT[3] = {0, 0, 0}; // metro dist uses default dist
	Vertex<T> *path = NULL;
	Vertex<T> *pathT[3] = {NULL, NULL, NULL}; // metro path uses default path

	int queueIndex = 0; 		// required by MutablePriorityQueue
	bool processing = false;

	string type = "";

public:
	Vertex(T in);
	Vertex(T in, string type);
	bool operator<(Vertex<T> & vertex) const; // // required by MutablePriorityQueue
	T getInfo() const;
	double getDist() const;
	Vertex *getPath() const;
	friend class Graph<T>;
	friend class MutablePriorityQueue<Vertex<T>>;
};


template <class T>
Vertex<T>::Vertex(T in): info(in) {}

template <class T>
Vertex<T>::Vertex(T in, string type): info(in), type(type) {}
/*
 * Auxiliary function to add an outgoing edge to a vertex (this),
 * with a given destination vertex (d) and edge weight (w).
 */

template <class T>
bool Vertex<T>::operator<(Vertex<T> & vertex) const {
	return this->dist < vertex.dist;
}

template <class T>
T Vertex<T>::getInfo() const {
	return this->info;
}

template <class T>
double Vertex<T>::getDist() const {
	return this->dist;
}

template <class T>
Vertex<T> *Vertex<T>::getPath() const {
	return this->path;
}

/********************** Edge  ****************************/

template <class T>
class Edge {
	Vertex<T> * dest;      	// destination vertex
	double dist;
	double avg_speed;
	double price;
	double weight; 			// edge weight, default = time
	boolean transport_transfer = false;

public:
	Edge(Vertex<T> *d, double dist, double vel, double price, bool transfer);
	Edge(Vertex<T> *d, double dist, double vel, double price);
	void setWeight(double w);
	friend class Graph<T>;
	friend class Vertex<T>;
};

template <class T>
Edge<T>::Edge(Vertex<T> *d, double dist, double vel, double price, bool transfer): dest(d), dist(dist), avg_speed(vel), price(price), weight(dist/vel), transport_transfer(transfer) {}
template <class T>
Edge<T>::Edge(Vertex<T> *d, double dist, double vel, double price): dest(d), dist(dist), avg_speed(vel), price(price), weight(dist/vel) {}

template <class T>
void Edge<T>::setWeight(double w) {
	weight = w;
}

/*************************** Graph  **************************/

template <class T>
class Graph {
	vector<Vertex<T> *> vertexSet;    // vertex set

public:
	Vertex<T> *findVertex(const T &in) const;
	bool addVertex(const T &in);
	bool addVertex(const T &in, string type);
	bool addEdge(const T &sourc, const T &dest, double dist, double vel, double price);
	bool addEdge(const T &sourc, const T &dest, double dist, double vel, double price, bool transfer);
	bool addUEdge(const T &sourc, const T &dest, double dist, double vel, double price);
	bool addStationEdge(const T &sourc, const T &dest, double dist, double vel, double price);
	int getNumVertex() const;
	vector<Vertex<T> *> getVertexSet() const;

	// Fp05 - single source
	void dijkstraShortestPath(const T &s);
	void dijkstraShortestPath(const T &s, const T &dest);

	void dijkstraShortestPathOld(const T &s);
	void unweightedShortestPath(const T &s);
	void bellmanFordShortestPath(const T &s);

	vector<T> getPath(const T &origin, const T &dest) const;
	double getDist(const T &origin, const T &dest) const;
	void setWeights(const string pref);

	// Fp05 - all pairs
	void floydWarshallShortestPath();
	vector<T> getfloydWarshallPath(const T &origin, const T &dest) const;

	// Graph Display
	void displayGraph(GraphViewer *gv);
	void displayPath(GraphViewer *gv, const T &origin, const T &dest);

	//String Matching
	int levDistance(const string source, const string target);
	string approximateStringMatching(const string source);
};

template <class T>
int Graph<T>::getNumVertex() const {
	return vertexSet.size();
}

template <class T>
vector<Vertex<T> *> Graph<T>::getVertexSet() const {
	return vertexSet;
}

/*
 * Auxiliary function to find a vertex with a given content.
 */
template <class T>
Vertex<T> * Graph<T>::findVertex(const T &in) const {
	for (auto v : vertexSet)
		if (v->info == in)
			return v;
	return NULL;
}

/*
 *  Adds a vertex with a given content or info (in) to a graph (this).
 *  Returns true if successful, and false if a vertex with that content already exists.
 */
template <class T>
bool Graph<T>::addVertex(const T &in) {
	if ( findVertex(in) != NULL)
		return false;
	vertexSet.push_back(new Vertex<T>(in));
	return true;
}

template <class T>
bool Graph<T>::addVertex(const T &in, string type) {
	if ( findVertex(in) != NULL)
		return false;
	vertexSet.push_back(new Vertex<T>(in, type));
	return true;
}

/*
 * Adds an edge to a graph (this), given the contents of the source and
 * destination vertices and the edge weight (w).
 * Returns true if successful, and false if the source or destination vertex does not exist.
 */
template <class T>
bool Graph<T>::addEdge(const T &sourc, const T &dest, double dist, double vel, double price) {
	auto v1 = findVertex(sourc);
	auto v2 = findVertex(dest);
	if (v1 == NULL || v2 == NULL)
		return false;

	v1->adj.push_back(new Edge<T>(v2, dist, vel, price));
	return true;
}

template <class T>
bool Graph<T>::addEdge(const T &sourc, const T &dest, double dist, double vel, double price, bool transfer) {
	auto v1 = findVertex(sourc);
	auto v2 = findVertex(dest);
	if (v1 == NULL || v2 == NULL)
		return false;

	v1->adj.push_back(new Edge<T>(v2, dist, vel, price, transfer));
	return true;
}


/*
 * Adds two edges to a graph (this), given the contents of the source and
 * destination vertices and the edge weight (w), going both ways (source -> destination and destination -> source).
 * Returns true if successful, and false if the source or destination vertex does not exist.
 */
template <class T>
bool Graph<T>::addUEdge(const T &sourc, const T &dest, double dist, double vel, double price) {
	return (addEdge(sourc, dest, dist, vel, price) && addEdge(dest, sourc, dist, vel, price));
}

template <class T>
bool Graph<T>::addStationEdge(const T &sourc, const T &dest, double dist, double vel, double price) {
	return (addEdge(sourc, dest, dist, vel, price, true) && addEdge(dest, sourc, dist, vel, price, true));
}

/**************** Single Source Shortest Path algorithms ************/

template<class T>
void Graph<T>::dijkstraShortestPath(const T &origin) {
	dijkstraShortestPath(origin, "");
}

template<class T>
void Graph<T>::dijkstraShortestPath(const T &origin, const T &dest) {
	// Setup Graph (Reset paths and distances to new point of origin) and Queue
	MutablePriorityQueue<Vertex<T>> priorityQueue = MutablePriorityQueue<Vertex<T>>();

	for (auto vertex : vertexSet) {
		vertex->dist = INF;
		vertex->path = NULL;
		priorityQueue.insert(vertex);
	}

	// Setup origin (whose distance is 0)
	Vertex<T>* originV = findVertex(origin);
	if (originV == NULL)	return;

	originV->dist = 0;
	priorityQueue.decreaseKey(originV);


	// Iterate each vertex starting from the closest to the origin to the farthest
	while(!priorityQueue.empty()) {
		Vertex<T> *v = priorityQueue.extractMin();
		if (v->info == dest)	return; // Already got shortest path to dest, so no need to continue processing

		for (auto edge : v->adj) {
			if (edge->dest->dist > v->dist + edge->weight) {
				edge->dest->dist = v->dist + edge->weight;
				edge->dest->path = v;
				priorityQueue.decreaseKey(edge->dest);
			}
		}
	}
}

template<class T>
vector<T> Graph<T>::getPath(const T &origin, const T &dest) const{
	vector<T> res;

	Vertex<T>* v = findVertex(dest);
	while(v != NULL) {
		res.insert(res.begin(), v->info);
		v = v->path;
	}

	return res;
}

template<class T>
double Graph<T>::getDist(const T &origin, const T &dest) const{
	return findVertex(dest)->dist;
}

template <class T>
void Graph<T>::setWeights(const string pref) {
	for (auto vertex : vertexSet)  {
		for (auto edge : vertex->adj) {
			if 		(pref == "price")		edge->setWeight(edge->price);
			else if (pref == "time")		edge->setWeight(edge->dist / edge->avg_speed);
			else if (pref == "distance")	edge->setWeight(edge->dist);
			else if (pref == "transfer")	edge->setWeight(edge->transport_transfer ? 1 : 0);
			else							edge->setWeight(0.60 * edge->dist / edge->avg_speed + 0.30 * edge->price + 0.10 * (edge->transport_transfer ? 1 : 0)); // default function

		}
	}
}



template<class T>
void Graph<T>::unweightedShortestPath(const T &orig) {
	// TODO
}

template<class T>
void Graph<T>::bellmanFordShortestPath(const T &orig) {
	// TODO
}


/**************** All Pairs Shortest Path  ***************/

template<class T>
void Graph<T>::floydWarshallShortestPath() {
	// TODO
}

template<class T>
vector<T> Graph<T>::getfloydWarshallPath(const T &orig, const T &dest) const{
	vector<T> res;
	// TODO
	return res;
}

/*************** Graph Display **************/

template <class T>
void Graph<T>::displayGraph(GraphViewer *gv) {
	//gv->defineVertexColor("blue");
	gv->defineEdgeColor("black");

	int nodeC = 0, edgeC = 0;
	vector<string> vertexId;

	// Add all vertices to GraphViewer
	for (auto vertex : vertexSet) {
		vertex->processing = false;

		gv->addNode(nodeC);
		gv->setVertexLabel(nodeC, vertex->info);

		vertexId.push_back(vertex->info);

		if(vertex->type == "M")
			gv->setVertexColor(nodeC, "blue");
		else if(vertex->type == "T")
			gv->setVertexColor(nodeC, "green");
		else if(vertex->type == "B")
			gv->setVertexColor(nodeC, "red");
		else
			gv->setVertexColor(nodeC, "gray");

		nodeC++;
	}

	// Add edges to GraphViewer
	nodeC = 0;
	for (auto vertex : vertexSet) {
		for (auto edge : vertex->adj) {
			if (edge->dest->processing)	continue; // if already been processed, ignore since node is already added

			int dest = find(vertexId.begin(), vertexId.end(), edge->dest->info) - vertexId.begin();

			gv->addEdge(edgeC, nodeC, dest, EdgeType::UNDIRECTED);
			gv->setEdgeWeight(edgeC, edge->weight);


			if(vertex->type == "M" && edge->dest->type == "M")
				gv->setEdgeColor(edgeC, "blue");
			else if(vertex->type == "T" && edge->dest->type == "T")
				gv->setEdgeColor(edgeC, "green");
			else if(vertex->type == "B" && edge->dest->type == "B")
				gv->setEdgeColor(edgeC, "red");
			else
				gv->setEdgeColor(edgeC, "gray");
			edgeC++;
		}
		vertex->processing = true;
		nodeC++;
	}

	gv->rearrange();
}

template <class T>
void Graph<T>::displayPath(GraphViewer *gv, const T &origin, const T &dest) {
	vector<T> path = this->getPath(origin, dest);

	for (auto string : path) {
		Vertex<T> *vertex = findVertex(string);
		gv->setVertexColor(find(vertexSet.begin(), vertexSet.end(), vertex) - vertexSet.begin(), "ORANGE");
	}

	gv->rearrange();
}


template<class T>
class TransportGrid: public Graph<T> {
	vector<Route *> routes;
	double metroDist, busDist, trainDist;

	double getSpeed(string type);

	void *metroWrapper(void *queue);
	void *busWrapper(void *queue);
	void *trainWrapper(void *queue);
	double improvedAlgorithmThread(const Vertex<T> *origin, const Vertex<T> *dest, MutablePriorityQueue<Vertex<T>> *priorityQueue, string type, int iT);

	Route* findRoute(const string source);
public:
	void setRoutes(vector<Route *> r);
	vector<Route *> getRoutes();

	MutablePriorityQueue<Stop> ApproximateStringMatching(const string source);
	bool ExactStringMatching(const string source);

	vector <Route *> getRoutesWith(const string source);
	string ShowInfo(const string source, int type);

	bool addStation(const T &in, const string mode);
	bool addConnection(const T &origin, const T &dest, double dist, double price, string type);
	bool addUConnection(const T &origin, const T &dest, double dist, double price, string type);

	double dijkstraAlgorithm(const T &origin, const T &dest);
	double improvedAlgorithm(const T &origin, const T &dest);
};

template <class T>
Route* TransportGrid<T>::findRoute(const string source) {
	for (auto route : routes) {
		if (route->getName() == source)		return route;
	}
	return NULL;
}

template <class T>
void TransportGrid<T>::setRoutes(vector<Route *> r) {routes = r;}

template <class T>
vector<Route *> TransportGrid<T>::getRoutes(){return routes;}

/*
template<class T>
MutablePriorityQueue<Stop> TransportGrid<T>::approximateStringMatching(const string source)
{
	vector<Route *> routes = this->getRoutes();
	vector<nameAndCost> ncvec;
	for(unsigned int i = 0; i < routes.size(); i++){
		ncvec.push_back(routes[i]->approximateStringMatching(source));
	}
	sort(ncvec.begin(), ncvec.end(), compareByCost);
	// Select top 3
	vector <string> top;
	for(int i = 0; i < 4; i++){
		top.push_back(ncvec[i].name);
	}
	return top;
}
*/
template<class T>
MutablePriorityQueue<Stop> TransportGrid<T>::ApproximateStringMatching(const string source) {
	MutablePriorityQueue<Stop> q;

	for (auto route : getRoutes()) {
		vector<Stop> v = route->possibleStops(source); // possibleStops() also tests the source against the route name

		for (auto stop : v) {
			q.insert(new Stop(stop));
		}
	}

	return q;
}



template <class T>
bool TransportGrid<T>::ExactStringMatching(const string source) {
	for (auto route : getRoutes())
		if (route->hasStop(source))	return true; // hasStop() also tests the source against the route name

	return false;
}

template <class T>
vector<Route *> TransportGrid<T>::getRoutesWith(const string source){
	vector<Route *> v;

	for (auto route : getRoutes())
		if (route->hasStop(source))	v.push_back(route);

	return v;
}

template <class T>
string TransportGrid<T>::ShowInfo(const string source, int type) {
	// Type: 0 - Origin Choice, 1 - Destination Choice

	string input;
	Route* r = findRoute(source);

	// ROUTE
	if (r != NULL) {
		cout << "## " << source << " (Route) ##\n\n";
		cout << "> This route passes through the following stops:" << endl;

		vector<Stop> stops = r->getStops();
		int i = 1;
		for (auto stop : stops) {
			cout << i << ". " << stop.getName() << "\n";
			i++;
		}
		cout << "\n> Type the stop number to get more info:" << endl;

		cin.clear();
		cin >> input;

		unsigned int val = stoi(input, nullptr);

		if (val > 0 && val <= stops.size())		return ShowInfo(stops[val - 1].getName(), type);
		else									return "";
	}

	// STOP
	else {
		cout << "## " << source << " (Stop) ##\n\n";
		cout << "> This stop passes through the following routes:\n";

		vector<Route *> routes = getRoutesWith(source);
		int i = 1;
		for (auto route : routes) {
			cout << i << "." << route->getName() << "\n";
			i++;
		}
		cout << "\n> If you wish to choose this stop press '0'.\n";
		cout << "> Otherwise type the route number to get more info:" << endl;

		cin.clear();
		cin >> input;

		unsigned int val = stoi(input, nullptr);

		if (val == 0)								return source;
		else if (val > 0 && val <= routes.size())	return ShowInfo(routes[val - 1]->getName(), type);
		else										return "";
	}
}

template <class T>
double TransportGrid<T>::getSpeed(string type) {
	if 		(type == "T")	return 100;
	else if	(type == "M")	return 40;
	else 					return 30;
}

/*
 * Adds a station supporting various types of transport specified in a string (mode)
 * (M) - Metro
 * (B) - Bus
 * (T) - Train
 *
 * Graph type T must support string conversion (or be a string) and it is assumed that no specific element was added
 * Ex: If "stationX(M)" doesn't exist, it assumes "stationX(B)" and "stationX(T)" don't exist either
 */
template <class T>
bool TransportGrid<T>::addStation(const T &in, const string mode) {
	const bool hasMetro = (mode.find('M', 0) != mode.npos);
	const bool hasBus = (mode.find('B', 0) != mode.npos);
	const bool hasTrain = (mode.find('T', 0) != mode.npos);

	if (!Graph<T>::addVertex(in, "central")) 	return false;


	if (hasMetro) {
		if (!Graph<T>::addVertex(in + "(M)", "M"))	return false;
		Graph<T>::addStationEdge(in + "(M)", in, 50, 1, 0);
	}
	if (hasBus) {
		if (!Graph<T>::addVertex(in + "(B)", "B"))	return false;
		Graph<T>::addStationEdge(in + "(B)", in, 10, 1, 0);
	}
	if (hasTrain) {
		if (!Graph<T>::addVertex(in + "(T)", "T"))	return false;
		Graph<T>::addStationEdge(in + "(T)", in, 60, 1, 0);
	}

	return true;
}

template <class T>
bool TransportGrid<T>::addConnection(const T &origin, const T &dest, double dist, double price, string type) {
	return Graph<T>::addEdge(origin + "(" + type + ")", dest + "(" + type + ")", dist, getSpeed(type), price);
}

template <class T>
bool TransportGrid<T>::addUConnection(const T &origin, const T &dest, double dist, double price, string type) {
	return Graph<T>::addUEdge(origin + "(" + type + ")", dest + "(" + type + ")", dist, getSpeed(type), price);
}


template <class T>
double TransportGrid<T>::dijkstraAlgorithm(const T &origin, const T &dest) {
	Graph<T>::dijkstraShortestPath(origin, dest);
	return Graph<T>::getDist(origin, dest);
}

template <class T>
void *TransportGrid<T>::metroWrapper(void *queue) {
	Vertex<T> *origin = ((MutablePriorityQueue<Vertex<T>>) queue)->extractMin();
	Vertex<T> *dest = ((MutablePriorityQueue<Vertex<T>>) queue)->extractMin();

	((MutablePriorityQueue<Vertex<T>>) queue)->insert(origin);
	dest->dist = INF;
	((MutablePriorityQueue<Vertex<T>>) queue)->insert(dest);

	double metroDist = improvedAlgorithmThread(origin, dest, queue, "M");
	return NULL;
}

template <class T>
void *TransportGrid<T>::busWrapper(void *queue) {
	Vertex<T> *origin = ((MutablePriorityQueue<Vertex<T>>) queue)->extractMin();
	Vertex<T> *dest = ((MutablePriorityQueue<Vertex<T>>) queue)->extractMin();

	((MutablePriorityQueue<Vertex<T>>) queue)->insert(origin);
	dest->dist = INF;
	((MutablePriorityQueue<Vertex<T>>) queue)->insert(dest);

	double busDist = improvedAlgorithmThread(origin, dest, queue, "B");
	return NULL;
}

template <class T>
void *TransportGrid<T>::trainWrapper(void *queue) {
	Vertex<T> *origin = ((MutablePriorityQueue<Vertex<T>>) queue)->extractMin();
	Vertex<T> *dest = ((MutablePriorityQueue<Vertex<T>>) queue)->extractMin();

	((MutablePriorityQueue<Vertex<T>>) queue)->insert(origin);
	dest->dist = INF;
	((MutablePriorityQueue<Vertex<T>>) queue)->insert(dest);

	double trainDist = improvedAlgorithmThread(origin, dest, queue, "T");
	return NULL;
}



template <class T>
double TransportGrid<T>::improvedAlgorithm(const T &origin, const T &dest) {
	// Graph Pre-processing
	MutablePriorityQueue<Vertex<T>> priorityQueue = MutablePriorityQueue<Vertex<T>>();

	priorityQueue.insert(findVertex(dest));
	for (auto vertex : Graph<T>::vertexSet) {
		if 		(vertex == findVertex(origin))	vertex->dist = 0;
		else if (vertex == findVertex(dest))	vertex->dist = 1;
		else									vertex->dist = INF;

		vertex->path = NULL;
		priorityQueue.insert(vertex);
	}

	// call each thread to process the 3 greedy paths
	pthread_t metroT, busT, trainT;
	pthread_create(&metroT, NULL, metroWrapper, (void *) priorityQueue);
	pthread_create(&busT, NULL, busWrapper, (void *) priorityQueue);
	pthread_create(&trainT, NULL, trainWrapper, (void *) priorityQueue);

	pthread_join(metroT, NULL);
	pthread_join(busT, NULL);
	pthread_join(trainT, NULL);

	// get shortest dist and path
	double minDist = min(metroDist, min(busDist, trainDist));

	int iT;
	if 		(minDist == metroDist)	iT = 0;
	else if (minDist == busDist)	iT = 1;
	else if (minDist == trainDist)	iT = 2;

	Vertex<T>* v = findVertex(dest);
	while(v != NULL) {
		v->path = v->path[iT];
		v = v->path;
	}

	return minDist;
}

template <class T>
double TransportGrid<T>::improvedAlgorithmThread(const Vertex<T> *origin, const Vertex<T> *dest, MutablePriorityQueue<Vertex<T>> *priorityQueue, string type, int iT) {
	// check if origin and dest have connection to requested transportation type
	Vertex<T> *originType = NULL, *destType = NULL;
	double currDist = 0;

	for (auto edge : origin->adj) {
		if (edge->dest->type == type)	{
			originType = edge->dest;
			currDist += edge->weight;
			edge->dest->pathT[iT] = origin;
		}
	}
	for (auto edge : dest->adj) {
		if (edge->dest->type == type)	{
			destType = edge->dest;
			currDist += edge->weight;
			dest->pathT[iT] = edge->dest;
		}
	}

	// if both stations have the requested transportation type, use dijkstra on that specific transport
	while (originType == NULL || destType == NULL) {
		// otherwise look for the nearest stations that support that transportation type

		// find nearest from origin
		MutablePriorityQueue<Vertex<T>> secQueue = MutablePriorityQueue<Vertex<T>>();

		for (auto vertex : Graph<T>::vertexSet) {
			vertex->distT[iT] = INF;
			vertex->pathT[iT] = NULL;
			secQueue.insert(vertex);
		}

		// Setup origin (whose distance is 0)
		origin->distT[iT] = 0;
		secQueue.decreaseKey(origin);


		// Iterate each vertex starting from the closest to the origin to the farthest
		while(!secQueue.empty()) {
			Vertex<T> *v = secQueue.extractMin();
			if (v->type == type) {
				originType = v;
				currDist += v->distT[iT];
			}
			if (v->info == dest)	return currDist + v ->distT[iT]; // Already got shortest path to dest, so no need to continue processing

			for (auto edge : v->adj) {
				if (edge->dest->distT[iT] > v->distT[iT] + edge->weight) {
					edge->dest->distT[iT] = v->distT[iT] + edge->weight;
					edge->dest->pathT[iT] = v;
					secQueue.decreaseKey(edge->dest);
				}
			}
		}

		while(!secQueue.empty()) {
			secQueue.extractMin();
		}

		// find nearest from dest
		for (auto vertex : Graph<T>::vertexSet) {
			vertex->distT[iT] = INF;
			secQueue.insert(vertex);
		}

		// Setup origin (whose distance is 0)
		dest->distT[iT] = 0;
		secQueue.decreaseKey(dest);


		// Iterate each vertex starting from the closest to the origin to the farthest
		while(!secQueue.empty()) {
			Vertex<T> *v = secQueue.extractMin();
			if (v->type == type) {
				destType = v;
				currDist += v->distT[iT];
			}
			if (v->info == dest)	return currDist + v ->distT[iT]; // Already got shortest path to dest, so no need to continue processing

			for (auto edge : v->adj) {
				if (edge->dest->distT[iT] > v->distT[iT] + edge->weight) {
					edge->dest->distT[iT] = v->distT[iT] + edge->weight;
					v->pathT[iT] = edge->dest;
					secQueue.decreaseKey(edge->dest);
				}
			}
		}
	}

	return currDist + singleTransportationDijkstra(originType, destType, priorityQueue, type, iT);
}

template <class T>
double singleTransportationDijkstra(const Vertex<T> *origin, const Vertex<T> *dest, MutablePriorityQueue<Vertex<T>> *priorityQueue, string type, int iT) {
	// Setup origin (whose distance is 0)
	origin->distT[iT] = 0;
	priorityQueue->decreaseKey(origin);

	// Iterate each vertex starting from the closest to the origin to the farthest
	while(!priorityQueue->empty()) {
		Vertex<T> *v = priorityQueue->extractMin();
		if (v->type != type)	continue;
		if (v->info == dest)	return v ->distT[iT]; // Already got shortest path to dest, so no need to continue processing

		for (auto edge : v->adj) {
			if (edge->dest->distT[iT] > v->distT[iT] + edge->weight) {
				edge->dest->distT[iT] = v->distT[iT] + edge->weight;
				edge->dest->pathT[iT] = v;
				priorityQueue->decreaseKey(edge->dest);
			}
		}
	}
}

#endif /* GRAPH_H_ */
