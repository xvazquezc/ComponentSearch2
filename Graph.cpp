#include "Graph.hpp"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <utility>


MetaNode & MetaGraph::getNodeFromIdx (int idx) {
	int position = positionOf(idx);
	if (position < 0) {
		std::cerr << "Meta node " << idx << " is missing from the graph" << std::endl;
		std::exit(1);
	}
	return nodes[position];
}

int MetaGraph::positionOf (int idx) const {
	if (idx < 0 || (size_t) idx >= positions.size())
		return -1;
	return positions[idx];
}

void MetaGraph::indexNodes () {
	int maxIdx = -1;
	for (const MetaNode & mn : nodes)
		maxIdx = std::max(maxIdx, mn.idx);

	positions.assign(maxIdx + 1, -1);
	for (size_t position = 0 ; position < nodes.size() ; position++)
		positions[nodes[position].idx] = position;
}

void MetaGraph::removeMarked (const std::vector<char> & removed) {
	size_t kept = 0;
	for (size_t position = 0 ; position < nodes.size() ; position++) {
		if (!removed[position]) {
			if (kept != position)
				nodes[kept] = std::move(nodes[position]);
			kept++;
		}
	}
	nodes.resize(kept);
	indexNodes();
}

long MetaGraph::getEdgesNb () const {
	long edgesNb = 0;
	for (const MetaNode & mn : nodes)
		edgesNb += mn.neighbors.size();
	return edgesNb / 2;
}
