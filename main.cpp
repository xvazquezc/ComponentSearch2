#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Graph.hpp"
#include "io.hpp"
#include "algorithms.hpp"
#include "splicing.hpp"


using namespace std;


int main (int argc, char * argv[]) {

	// --- Default parameters ---
	string nodesFilename = "nodes.csv";
	string edgesFilename = "edges.csv";
	string basename = "";
	int nodeFilter = 0;
	int edgeFilter = 0;
	int seed = time(0);
	bool optimizeComponents = false;


	// --- Arguments parsing ---
	for (int idx=1 ; idx<argc ; idx++) {
		if (argv[idx][0] == '-') {
			switch (argv[idx][1]) {
			case 'n':
				nodesFilename = string(argv[++idx]);
				break;
			case 'e':
				edgesFilename = string(argv[++idx]);
				break;
			case 'b':
				basename = string(argv[++idx]);
				break;
			case 'N':
				nodeFilter = atoi(argv[++idx]);
				break;
			case 'E':
				edgeFilter = atoi(argv[++idx]);
				break;
			case 'o':
				optimizeComponents = true;
				break;
			case 'r':
				seed = atoi(argv[++idx]);
				break;

			default:
				break;
			}
		}
	}
	// Parsed but unused, as upstream: filterEdges gets -N and splice never
	// optimizes. Kept so that results don't change.
	(void) edgeFilter;
	(void) optimizeComponents;

	// Init random generator
	srand(seed);


	// --- Loading ---
	cout << endl << "--- Loading ---" << endl;

	{
		ifstream nodesFile (nodesFilename.c_str());
		if (!nodesFile) {
			cerr << "The file " << nodesFilename << " seems to be missing." << endl;
			return 1;
		}

		ifstream edgesFile (edgesFilename.c_str());
		if (!edgesFile) {
			cerr << "The file " << edgesFilename << " seems to be missing." << endl;
			return 1;
		}
	}

	ReadGraph graph = load_graph(nodesFilename, edgesFilename);
	int readsNb = graph.size();
	cout << readsNb << " nodes loaded " << endl;


	// --- Algorithms ---
	cout << endl << "--- Run algorithms ---" << endl;
	cout << "-> BFS..." << endl;
	vector<int> bfs = BFS (graph);

	cout << "-> Contraction..." << endl;
	MetaGraph metagraph = contract (graph, bfs);
	vector<int>().swap(bfs);
	cout << metagraph.nodes.size() << " nodes in the contracted graph" << endl;
	cout << metagraph.getEdgesNb() << " edges in the contracted graph" << endl;

	cout << "-> Absorb fingers..." << endl;
	absorbFingers (metagraph);
	cout << metagraph.nodes.size() << " nodes in the contracted graph" << endl;
	cout << metagraph.getEdgesNb() << " edges in the contracted graph" << endl;

	cout << endl << "--- Filtering ---" << endl;
	cout << "-> Filtering nodes..." << endl;
	filterNodes (metagraph, nodeFilter);
	cout << metagraph.nodes.size() << " nodes in the filtered graph" << endl;

	cout << "-> Filtering edges..." << endl;
	filterEdges (metagraph, nodeFilter, graph);
	cout << metagraph.getEdgesNb() << " edges in the filtered graph" << endl;

	// The read links are not needed any more.
	graph = ReadGraph();


	// --- Splicing ---
	cout << endl << "--- Splicing ---" << endl;
	splice (metagraph, false);


	// --- Output ---
	cout << endl << "--- Saving results ---" << endl;
	cout << "-> Saving components" << endl;
	stringstream ss;
	ss << basename << ".components.csv";
	ofstream outStream (ss.str());
	ifstream inStream (nodesFilename);
	save_componants (metagraph, readsNb, inStream, outStream);

	cout << "-> Saving meta graph" << endl;
	stringstream nodes;
	stringstream edges;
	nodes << basename << ".metaNodes.csv";
	edges << basename << ".metaEdges.csv";
	ofstream nodesStream (nodes.str());
	ofstream edgesStream (edges.str());
	save_metagraph (metagraph, nodesStream, edgesStream);

	cout << endl << "--- Program ended ---" << endl << endl;

	return 0;
}
