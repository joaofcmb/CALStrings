#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <chrono>
#include <Windows.h>
#include <vector>
#include <set>

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
	while(1) {
		// Choose what to do (Test Preset, Procedurally Generated Graph Test, etc..)
		cout << "-- What do you want to do? --\n";
		cout << "-- 1. Porto Preset Test --\n";
		cout << "-- 2. Procedurally Generated Graph Test --\n";
		cout << "-- 3. Exit Application --"<< endl;

		switch(getchar()) {
		case '1':
			presetTest("Porto");
			break;
		case '2':
			genTest();
			break;
		case '3':
			return 0;
		default:
			cout << "Number out of range\n\n";
			continue;
		}
	}
}

string AproximateStringMenu(string input, TransportGrid<string> *preset) {
	MutablePriorityQueue<Stop> q = preset->ApproximateStringMatching(input);

	if (q.empty()) {
		cout << "No match for " << input << " found.\n\n";
		return "";
	}

	Stop* firstStop = q.extractMin();

	if (firstStop->getCost() == 0)		cout << "Found match for " << firstStop->getName() << ".\n";
	else								cout << "Closest match was " << firstStop->getName() << ".\n";

	if (!q.empty())						cout << "\nHowever, the following entries were also considered:" << endl;

	// Take rest of unique stops and put them in a vector for easy access
	vector<string> candidates;
	set<string> dupCheck;

	while (!q.empty()) {
		string s = q.extractMin()->getName();

		if (dupCheck.insert(s).second && s != firstStop->getName())
			candidates.push_back(s);
	}

	int i = 1;
	for (string name : candidates)		cout << i++ << ". " << name << endl;

	cout << "\nPress '0' to get more info on the initial match." << endl;
	if (candidates.size() > 0)			cout << "To get more info on the other entries instead, enter the corresponding number." << endl;

	cin.clear();
	cin >> input;

	unsigned int val = stoi(input, nullptr);

	if (val == 0)									input = firstStop->getName();
	else if (val > 0 && val <= candidates.size())	input = candidates[val - 1];
	else											return "";

	input = preset->ShowInfo(input, 1);

	return input;
}


void presetTest(string presetId) {
	TransportGrid<string>* preset = createGraphPreset("Presets\\" + presetId + ".txt");
	GraphViewer *gv = new GraphViewer(1200, 800, true);
	gv->createWindow(1200, 800);

	preset->displayGraph(gv);

	while (1) {
		char c;
		string input;
		string origin, dest;
		double dist;


		cout << "1. Approximate String Matching\n";
		cout << "2. Exact String Matching\n";
		cout << "Press any other character to exit" << endl;

		cin.clear();
		cin >> c;

		if (c < '1' || c > '2')		break;


		cout << "-- Search for Stop or Route (Origin): --" << endl;

		cin.ignore();
		cin.clear();
		getline(cin, input);

		if (c == '1') {
			origin = AproximateStringMenu(input, preset);

			if (origin == "")		continue;
		}
		else {
			if (!preset->ExactStringMatching(input)) {
				cout << "No match for " << input << " found.\n\n";
				continue;
			}
			input = preset->ShowInfo(input, 1);

			if (input == "") continue;

			origin = input;
		}

		cout << "-- Search for Stop or Route (Destination): --" << endl;

		cin.ignore();
		cin.clear();
		getline(cin, input);

		if (c == '1') {
			dest = AproximateStringMenu(input, preset);

			if (dest == "")		continue;
		}
		else {
			if (!preset->ExactStringMatching(input)) {
				cout << "No match for " << input << " found.\n\n";
				continue;
			}
			input = preset->ShowInfo(input, 1);

			if (input == "") continue;

			dest = input;
		}

		cout << "\nPath from "<< origin << " to " << dest << "." << endl;


		cout << "\n-- Choose your preference (1 - 4 | press any other key to exit) --\n";
		cout << "1. Lowest price\n" << "2. Shortest time\n" << "3. Shortest distance\n";
		cout << "4. Number of transport transfers\n" << "5. Default (considers various parameters)" << endl;

		cin.ignore();
		cin.clear();
		cin >> c;
		switch(c) {
		case '1':
			preset->setWeights("price");
			dist = preset->dijkstraAlgorithm(origin, dest);
			cout << "- Price: " << dist << " euros" << endl;
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
			continue;
		}

		preset->displayGraph(gv); // to update weights
		preset->displayPath(gv, origin, dest);

		cout << "\nPath calculated.\n" << endl;
	}
	return;
}

void generateRandomTransportGraph(int& n, TransportGrid<string>* g) {
	string type, type2;
	int randDist, randPrice;


	for(int i = 0; i < n; i++){
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

string checkString(vector<string> top4){

	string name;
	cout << "Did we get it right or did you mean: "
			<< endl << "{"
			<< endl << "	1." << top4[1]
										<< endl << "	2." << top4[2]
																	<< endl << "	3." << top4[3] << "?"
																	<< endl << "}"
																	<< endl
																	<< endl << "4. Yes we got it right :3" << endl;

	int input2;
	cin >> input2;

	switch(input2){
	case 1:
		name = top4[1];
		return name;
	case 2:
		name = top4[2];
		return name;
	case 3:
		name = top4[3];
		return name;
	case 4:
		break;
	}
	return name;
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

				if (r->getStops().size() == 0)	r->addStop(arg1);
				r->addStop(arg2);
			}
		}
		else if (command == "D_EDGE") {
			infile >> dist >> price>> type;
			g->addConnection(arg1, arg2, dist, price, type);

			if (routes.size() > 0) {
				auto r = routes.back();

				if (r->getStops().size() == 0)	r->addStop(arg1);
				r->addStop(arg2);
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
