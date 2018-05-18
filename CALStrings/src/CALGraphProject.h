#ifndef CALGRAPHPROJECT_H_
#define CALGRAPHPROJECT_H_

#include "Graph.h"

void presetTest(string presetId);
void genTest();

/*
 * File Syntax:
 *
 * Commands: STATION, U_EDGE, D_EDGE
 *
 * Ex:
 * "STATION S.Bento BTM" 	 - Creates a station called S.Bento with all the transport types (Bus, Train and Metro)
 * "STATION Trindade BM" 	 - Creates a station called Trindade supporting Bus and Metro
 * "U_EDGE S.Bento Trindade  2 M" - Creates a metro connection between S.Bento and Trindade both ways (Undirected) with a weight of 2
 * "D_EDGE S.Bento Trindade  5 M" - Creates a metro connection from S.Bento to Trindade (Directed) with a weight of 5
 */
TransportGrid<string> *createGraphPreset(string filename);

#endif /* CALGRAPHPROJECT_H_ */
