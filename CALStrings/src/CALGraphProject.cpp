#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <chrono>
#include <Windows.h>

#include "GraphViewer\GraphViewer.h"
#include "CALGraphProject.h"
#include "Graph.h"
#include "Route.h"

using namespace std;

/*
 * 3 Steps of Program
 *
 * 1 - Construct the graph
 *
 * 2 - Apply an algorithm to it to find the min distance between 2 points
 *
 * 3 - Display it using the GUI (Both the Graph and an highlighted path)
 */

int main() {

	char input;
	while(1){
	// This is where the algorithm would be chosen (and then put in a function pointer)
	cout << "-- What do you want to do? --" << endl;
	cout << "-- 1. Porto Preset Test --" << endl;
	cout << "-- 2. Procedurally Generated Graph Test --" << endl;

	cin >> input;
	// Choose what to do (Test Preset, Procedurally Generated Graph Test, etc..)
	if(input=='1')
		presetTest("Porto");
	else if(input =='2')
		genTest();
	else
		system("CLS");
		cout << "Choose '1' or '2' !!! " << endl;

	}


}

void presetTest(string presetId) {
	TransportGrid<string>* preset = createGraphPreset("Presets\\" + presetId + ".txt");
	GraphViewer *gv = new GraphViewer(1200, 800, true);
	gv->createWindow(1200, 800);

	preset->displayGraph(gv);

	while (1) {
		char input;
		string origin, dest;
		double dist;

		cout << "-- Choose your starting location: --" << endl;
		cin >> origin;
		cout << "-- Choose your destination: --" << endl;
		cin >> dest;
		cout << "-- Choose your preference (1 - 4 | press any other key to exit) --" << endl;
		cout << "1. Lowest price\n" << "2. Shortest time\n" << "3. Shortest distance\n"
				<< "4. Number of transport transfers\n"  <<"5. Default (considers various parameters)" << endl;

		cin >> input;

		switch(input) {
		case '1':
			preset->setWeights("price");
			dist = preset->dijkstraAlgorithm(origin, dest);
			cout << "- Price: " << dist <<"€" << endl;
			break;
		case '2':
			preset->setWeights("time");
			dist = preset->dijkstraAlgorithm(origin, dest)/60;
			cout << "- Total Time: " << dist << " minutes" << endl;
			break;
		case '3':
			preset->setWeights("distance");
			dist = preset->dijkstraAlgorithm(origin, dest);
			cout << "- Distance: " << dist << " meters" <<endl;
			break;
		case '4':
			preset->setWeights("transfer");
			dist = preset->dijkstraAlgorithm(origin, dest);
			cout << "- Number of vehicle changes: " << dist << endl;
			break;
		case '5':
			preset->setWeights("default");
			dist = preset->dijkstraAlgorithm(origin, dest);
			break;
		default:
			return;
		}


		preset->displayGraph(gv); // to update weights
		preset->displayPath(gv, origin, dest);

	}
	return;
}

void generateRandomTransportGraph(int& n, TransportGrid<string>* g) {
	string type, type2;
	int randType, randDist, randPrice, aux;


		for(int i = 0; i < n; i++){
		/*	randType = rand()%(7-1 + 1) + 1;
			switch(randType){
			case 1: type = "BTM"; break;
			case 2: type = "BM"; break;
			case 3: type = "TM"; break;
			case 4: type = "M"; break;
			}*/
			g->addStation(to_string(i), "M");
		}

		for(int i = 0; i < n; i++){

			randDist = rand()%(10000-300 + 1) + 300;
			randPrice = rand()%(70-1 + 1) + 1;
			g->addUConnection(to_string(i), to_string(i+1), randDist, randPrice, "M");
			/*
			 * if(g->findVertex(to_string(i)+"(B)")!=NULL)
				type = "B";
			else if(g->findVertex(to_string(i)+"(M)")!=NULL)
				type = "M";
			else
				type = "T";
			aux = i+1;
			while(aux < n){
				if(g->findVertex(to_string(aux)+type) != NULL){
					g->addUConnection(to_string(i)+type, to_string(aux)+type, randDist, randPrice, type);
				}
				aux++;
			}*/
		}
		g->addUConnection(to_string(n-1), to_string(0), randDist, randPrice, "M");
}

void genTest() {

	TransportGrid<string> *g = new TransportGrid<string>();
	/*int test;
	int test2=5;

	GraphViewer *gv = new GraphViewer(1200, 800, true);
	gv->createWindow(1200, 800);
	generateRandomTransportGraph(test2, g);
	g->displayGraph(gv);
	while(1){
		cin >> test;
	}*/

	for (int n = 100; n <= 1000; n += 100) {
		cout << "generating graph" << n << " x " << n << " ..." << endl;
		generateRandomTransportGraph(n, g);
		cout << "processing graph" << n << " x " << n << " ..." << endl;
		auto start = std::chrono::high_resolution_clock::now();

		// roda o algoritmo a partir de todos os pontos
		for(int i=0; i < n; i++){
			g->dijkstraShortestPath(to_string(i));
		}

		auto finish = std::chrono::high_resolution_clock::now();
		auto elapsed = chrono::duration_cast<chrono::nanoseconds>(finish - start).count();
		cout << "processing grid " << n << " x " << n << " average time (nano-seconds)=" << (elapsed / (n*n)) << endl;
	}
}

TransportGrid<string> *createGraphPreset(string filename) {
	TransportGrid<string> *g = new TransportGrid<string>();
	ifstream infile;
	infile.open(filename, ios::in);

	vector<Route *> routes;

	double dist, price;
	string command, arg1, arg2, type;
	while (infile >> command >> arg1 >> arg2) {
		if (command == "STATION")
			g->addStation(arg1, arg2);
		else if (command == "U_EDGE") {
			infile >> dist >> price >> type;
			g->addUConnection(arg1, arg2, dist, price, type);

			if (routes.size() > 0) {
				auto r = routes.back();
				r->addStop(arg2);

				if (r == routes[0])	r->addStop(arg1);
			}
		}
		else if (command == "D_EDGE") {
			infile >> dist >> price>> type;
			g->addConnection(arg1, arg2, dist, price, type);

			if (routes.size() > 0) {
				auto r = routes.back();
				r->addStop(arg2);

				if (r == routes[0])	r->addStop(arg1);
			}
		}
		else if (command == "ROUTE") {
			arg1 += " " + arg2;
			getline(infile, arg2);
			routes.push_back(new Route(arg1 + arg2));
		}
	}

	g->setRoutes(routes);

	infile.close();
	return g;
}
