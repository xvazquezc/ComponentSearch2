#include <fstream>
#include <string>

#include "Graph.hpp"


#ifndef IO_HPP
#define IO_HPP

/* One read per non-empty line of the nodes file after its header; the edges
   file is read twice (degrees, then neighbors) so that the graph is built in
   place, without a temporary edge list. */
ReadGraph load_graph (const std::string & nodesFilename, const std::string & edgesFilename);
void save_componants (MetaGraph & metagraph, int readsNb, std::ifstream & nodes, std::ofstream & componentsStream);
void save_metagraph (const MetaGraph & graph, std::ofstream & nodes, std::ofstream & edges);

#endif
