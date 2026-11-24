#include "algorithms.hpp"

#include <algorithm>
#include <cstdlib>


using namespace std;


namespace {

/* The reads not yet annotated by the BFS, with the k-th smallest found in
   O(log n) (a Fenwick tree), as advancing k steps in a std::set would. */
class RemainingReads {
public:
	explicit RemainingReads (int size) : tree(size + 1), size(size), topStep(1) {
		for (int i = 1 ; i <= size ; i++)
			tree[i] = i & -i;
		while (topStep * 2 <= size)
			topStep *= 2;
	}

	void remove (int read) {
		for (int i = read + 1 ; i <= size ; i += i & -i)
			tree[i]--;
	}

	int kth (size_t k) const {
		int position = 0;
		for (int step = topStep ; step > 0 ; step >>= 1) {
			if (position + step <= size && (size_t) tree[position + step] <= k) {
				position += step;
				k -= tree[position];
			}
		}
		return position;
	}

private:
	vector<int> tree;
	int size;
	int topStep;
};

}  // namespace


/* Levels of a breadth-first search from a random read of each connected
   component. The start read and its neighbors share a level; each component
   ends with an unused level number. */
vector<int> BFS (const ReadGraph & graph) {
	int readsNb = graph.size();
	vector<int> dfs (readsNb, -1);
	RemainingReads toAnnotate (readsNb);
	size_t toAnnotateNb = readsNb;

	vector<int> nextIdxs;
	vector<int> neighbors;
	int currentIdx = 1;
	while (toAnnotateNb > 0) {
		// Random init of the first node
		size_t offset = rand() % toAnnotateNb;
		int startIdx = toAnnotate.kth(offset);
		toAnnotate.remove(startIdx);
		toAnnotateNb--;

		// Init
		nextIdxs.assign(1, startIdx);
		dfs[startIdx] = currentIdx;

		// One level at a time; the order within a level changes nothing.
		while (nextIdxs.size() > 0) {
			neighbors.clear();

			for (int nodeIdx : nextIdxs) {
				for (const int * nei = graph.neighborsBegin(nodeIdx) ; nei != graph.neighborsEnd(nodeIdx) ; nei++) {
					if (dfs[*nei] == -1) {
						dfs[*nei] = currentIdx;
						neighbors.push_back(*nei);
						toAnnotate.remove(*nei);
						toAnnotateNb--;
					}
				}
			}

			nextIdxs.swap(neighbors);
			currentIdx++;
		}
	}

	return dfs;
}


/* One meta node per group of connected reads of the same BFS level, in the
   order of their smallest read; meta node links sorted by idx. */
MetaGraph contract (const ReadGraph & graph, const vector<int> & dfs) {
	int readsNb = graph.size();
	MetaGraph contracted;
	vector<int> affectations (readsNb, -1);

	// Meta nodes creation. The reads of a meta node are the same whatever
	// the order of the visit, and nothing depends on the order of subNodes.
	int metaIdx = 0;
	vector<int> currentSet;
	for (int idx=0 ; idx<readsNb ; idx++) {
		if (affectations[idx] == -1) {
			int annotation = dfs[idx];
			MetaNode mn;
			mn.idx = metaIdx++;

			affectations[idx] = mn.idx;
			currentSet.assign(1, idx);
			while (currentSet.size() > 0) {
				int nodeIdx = currentSet.back();
				currentSet.pop_back();
				mn.subNodes.push_back(nodeIdx);

				// Adding the neighbors
				for (const int * nei = graph.neighborsBegin(nodeIdx) ; nei != graph.neighborsEnd(nodeIdx) ; nei++) {
					if (affectations[*nei] == -1 && dfs[*nei] == annotation) {
						affectations[*nei] = mn.idx;
						currentSet.push_back(*nei);
					}
				}
			}

			contracted.nodes.push_back(std::move(mn));
		}
	}


	// MetaLinks creation
	vector<int> linkedFrom (metaIdx, -1);
	for (MetaNode & mn : contracted.nodes) {
		for (int read : mn.subNodes) {
			for (const int * nei = graph.neighborsBegin(read) ; nei != graph.neighborsEnd(read) ; nei++) {
				int metaNei = affectations[*nei];
				if (metaNei != mn.idx && linkedFrom[metaNei] != mn.idx) {
					linkedFrom[metaNei] = mn.idx;
					mn.neighbors.push_back(metaNei);
				}
			}
		}
		sort(mn.neighbors.begin(), mn.neighbors.end());
	}

	contracted.indexNodes();
	return contracted;
}


/* A meta node with a single neighbor that has more than two neighbors is
   merged into it. Neighbor counts are the current ones, in node order. */
void absorbFingers (MetaGraph & graph) {
	vector<char> toRemove (graph.nodes.size(), 0);

	// Check for removable meta nodes
	for (size_t position = 0 ; position < graph.nodes.size() ; position++) {
		MetaNode & mn = graph.nodes[position];
		if (mn.neighbors.size() == 1) {
			MetaNode & metaNei = graph.getNodeFromIdx(mn.neighbors.front());

			if (metaNei.neighbors.size() > 2) {
				metaNei.subNodes.insert(metaNei.subNodes.end(), mn.subNodes.begin(), mn.subNodes.end());
				metaNei.neighbors.erase(
					remove(metaNei.neighbors.begin(), metaNei.neighbors.end(), mn.idx),
					metaNei.neighbors.end()
				);

				toRemove[position] = 1;
			}
		}
	}

	// Remove metanodes
	graph.removeMarked(toRemove);
}


/* Meta nodes of at most nodeThreshold reads are removed, with their links. */
void filterNodes (MetaGraph & graph, int nodeThreshold) {
	vector<char> toRemove (graph.nodes.size(), 0);

	// Detect removable
	for (size_t position = 0 ; position < graph.nodes.size() ; position++) {
		MetaNode & mn = graph.nodes[position];

		// size_t against int, as upstream compares them
		if (mn.subNodes.size() <= (size_t) nodeThreshold) {
			// Remove links
			for (int neiIdx : mn.neighbors) {
				MetaNode & metaNei = graph.getNodeFromIdx(neiIdx);
				metaNei.neighbors.erase(
					remove(metaNei.neighbors.begin(), metaNei.neighbors.end(), mn.idx),
					metaNei.neighbors.end()
				);
			}

			toRemove[position] = 1;
		}
	}

	// Remove meta nodes filtered
	graph.removeMarked(toRemove);
}


/* The weight of an edge between two meta nodes is the number of read links
   between their reads. Edges of weight at most edgeThreshold are removed,
   from both meta nodes (the weights are the same from both sides); 0 keeps
   every edge.

   Upstream (ComponentSearch2 f705924) was given -N instead of -E, counted
   the links per outside read instead of per meta node, and removed the
   neighbor whose meta node number equalled such a read number, on one side
   only. On production graphs it removed none of the weight-1 edges. */
void filterEdges (MetaGraph & graph, int edgeThreshold, const ReadGraph & reads) {
	if (edgeThreshold <= 0)
		return;

	// Position of the meta node of each read (-1: filtered out)
	vector<int> readPosition (reads.size(), -1);
	for (size_t position = 0 ; position < graph.nodes.size() ; position++)
		for (int read : graph.nodes[position].subNodes)
			readPosition[read] = position;

	vector<long> weights (graph.nodes.size(), 0);
	for (size_t position = 0 ; position < graph.nodes.size() ; position++) {
		MetaNode & mn = graph.nodes[position];
		if (mn.neighbors.empty())
			continue;

		// Compute weights
		for (int read : mn.subNodes)
			for (const int * nei = reads.neighborsBegin(read) ; nei != reads.neighborsEnd(read) ; nei++) {
				int neiPosition = readPosition[*nei];
				if (neiPosition >= 0 && neiPosition != (int) position)
					weights[neiPosition]++;
			}

		// Deletion of edges under the threshold
		size_t kept = 0;
		for (int neiIdx : mn.neighbors) {
			int neiPosition = graph.positionOf(neiIdx);
			if (neiPosition >= 0 && weights[neiPosition] > edgeThreshold)
				mn.neighbors[kept++] = neiIdx;
		}
		mn.neighbors.resize(kept);

		for (int read : mn.subNodes)
			for (const int * nei = reads.neighborsBegin(read) ; nei != reads.neighborsEnd(read) ; nei++)
				if (readPosition[*nei] >= 0)
					weights[readPosition[*nei]] = 0;
	}
}
