#include <vector>

#include "Graph.hpp"


#ifndef ALGOS_HPP
#define ALGOS_HPP

std::vector<int> BFS (const ReadGraph & graph);
MetaGraph contract (const ReadGraph & graph, const std::vector<int> & dfs);
void absorbFingers (MetaGraph & graph);
void filterNodes (MetaGraph & graph, int nodeThreshold);
void filterEdges (MetaGraph & graph, int edgeThreshold, const ReadGraph & reads);

#endif
