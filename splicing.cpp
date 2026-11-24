#include "splicing.hpp"

#include <algorithm>
#include <vector>


using namespace std;


/* Hub meta nodes (more than two neighbors) are cut from their neighbors of
   at most two neighbors, in node order with the current neighbor counts. */
void splice (MetaGraph & graph, bool optimize) {
	vector<int> neiToRemove;
	for (MetaNode & mn : graph.nodes) {
		neiToRemove.clear();

		// Looks for hub nodes
		if (mn.neighbors.size() > 2) {
			for (int neiIdx : mn.neighbors) {
				MetaNode & nei = graph.getNodeFromIdx(neiIdx);

				// If the neighbor is not from the same component
				if (nei.neighbors.size() <= 2) {
					// Remove the neighbor to main node link
					nei.neighbors.erase(
						remove(nei.neighbors.begin(), nei.neighbors.end(), mn.idx),
						nei.neighbors.end()
					);

					neiToRemove.push_back(neiIdx);

					// optimize the nodes in the components that are not hubs
					if (optimize) {
						// TODO
					}
				}
			}
		}

		// Remove the main node to neighbor links (neighbors are distinct)
		for (int idx : neiToRemove) {
			mn.neighbors.erase(
				remove(mn.neighbors.begin(), mn.neighbors.end(), idx),
				mn.neighbors.end()
			);
		}
	}
}
