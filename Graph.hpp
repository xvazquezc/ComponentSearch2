#include <cstdint>
#include <vector>


#ifndef GRAPH_HPP
#define GRAPH_HPP

/* Read graph in compressed sparse row form: the neighbors of read i are
   targets[offsets[i]] .. targets[offsets[i + 1] - 1], in the order of the
   edges file (an edge i;j adds j to i, then i to j; a self-loop twice).
   It costs 8 bytes per edge and nothing else per edge. */
class ReadGraph {
public:
	std::vector<uint64_t> offsets;  // number of reads + 1
	std::vector<int> targets;

	int size () const {
		return offsets.empty() ? 0 : (int) (offsets.size() - 1);
	}
	const int * neighborsBegin (int read) const {
		return targets.data() + offsets[read];
	}
	const int * neighborsEnd (int read) const {
		return targets.data() + offsets[read + 1];
	}
};


/* A group of reads (subNodes, by read index) and its links to other meta
   nodes (neighbors, by meta node idx). */
class MetaNode {
public:
	int idx;
	std::vector<int> neighbors;
	std::vector<int> subNodes;
};


class MetaGraph {
public:
	std::vector<MetaNode> nodes;  // by increasing idx

	// The node with this idx; exits if there is none. Nodes marked for
	// removal stay reachable until removeMarked.
	MetaNode & getNodeFromIdx (int idx);
	// Position of the node with this idx in nodes, or -1.
	int positionOf (int idx) const;
	// Rebuilds the idx -> position index after nodes changed.
	void indexNodes ();
	// Removes the nodes marked (by position), keeping the others in order.
	void removeMarked (const std::vector<char> & removed);
	long getEdgesNb () const;

private:
	std::vector<int> positions;
};

#endif
