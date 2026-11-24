#include "io.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>


using namespace std;


namespace {

/* Lines of a file, as std::getline gives them (split on '\n', a last line
   without '\n' included), read through a large buffer. */
class LineReader {
public:
	explicit LineReader (const string & filename)
		: file(fopen(filename.c_str(), "rb")), buffer(1 << 24), begin(0), end(0), atEnd(false), carryUsed(false) {
		if (!file) {
			cerr << "The file " << filename << " can't be read." << endl;
			exit(1);
		}
	}

	~LineReader () {
		fclose(file);
	}

	// Next line, without its '\n'; false after the last one. The line stays
	// valid until the next call.
	bool next (const char * & line, size_t & length) {
		if (carryUsed) {
			carry.clear();
			carryUsed = false;
		}

		while (true) {
			char * start = buffer.data() + begin;
			char * newline = (char *) memchr(start, '\n', end - begin);
			if (newline) {
				begin = newline + 1 - buffer.data();
				if (carry.empty()) {
					line = start;
					length = newline - start;
				} else {
					carry.append(start, newline);
					line = carry.data();
					length = carry.size();
					carryUsed = true;
				}
				return true;
			}

			carry.append(start, end - begin);
			begin = end = 0;
			if (atEnd) {
				if (carry.empty())
					return false;
				line = carry.data();
				length = carry.size();
				carryUsed = true;
				return true;
			}

			end = fread(buffer.data(), 1, buffer.size(), file);
			if (end < buffer.size())
				atEnd = true;
		}
	}

	void rewind () {
		std::rewind(file);
		begin = end = 0;
		atEnd = false;
		carry.clear();
		carryUsed = false;
	}

private:
	FILE * file;
	vector<char> buffer;
	size_t begin, end;
	bool atEnd;
	string carry;
	bool carryUsed;
};


// Split a string into multiple strings using the delim char as a separator
vector<string> split (const string & line, char delim) {
	vector<string> strs;

	string token;
	stringstream ss;
	ss << line;
	while (ss) {
		getline(ss, token, delim);

		if (ss)
			strs.push_back(token);
	}

	return strs;
}


// Field number `field` of a ';'-separated line, read as atoi would read it.
// False if the line has fewer fields.
bool fieldValue (const char * line, size_t length, int field, long & value) {
	const char * p = line;
	const char * end = line + length;
	for (int f = 0 ; f < field ; f++) {
		p = (const char *) memchr(p, ';', end - p);
		if (!p)
			return false;
		p++;
	}

	while (p < end && isspace((unsigned char) *p))
		p++;
	bool negative = false;
	if (p < end && (*p == '+' || *p == '-')) {
		negative = (*p == '-');
		p++;
	}
	long v = 0;
	while (p < end && *p >= '0' && *p <= '9' && v <= 1000000000000L) {
		v = v * 10 + (*p - '0');
		p++;
	}
	value = negative ? -v : v;
	return true;
}


// Both ends of an edge line, checked against the number of reads.
void edgeEnds (const char * line, size_t length, int srcIdx, int tgtIdx, int readsNb, int & idx1, int & idx2) {
	long v1, v2;
	if (!fieldValue(line, length, srcIdx, v1) || !fieldValue(line, length, tgtIdx, v2)) {
		cerr << "Malformed line in the edges csv file: " << string(line, length) << endl;
		exit(1);
	}
	if (v1 < 0 || v1 >= readsNb || v2 < 0 || v2 >= readsNb) {
		cerr << "Edge between unknown nodes in the edges csv file: " << string(line, length) << endl;
		exit(1);
	}
	idx1 = v1;
	idx2 = v2;
}

}  // namespace


ReadGraph load_graph (const string & nodesFilename, const string & edgesFilename) {
	ReadGraph g;
	const char * line;
	size_t length;

	// First line verticies, then one node per non-empty line
	int readsNb = 0;
	{
		LineReader verticies(nodesFilename);
		if (verticies.next(line, length)) {
			while (verticies.next(line, length))
				if (length > 0)
					readsNb++;
		}
	}
	g.offsets.assign(readsNb + 1, 0);


	// first line of edges
	LineReader edges(edgesFilename);
	string header;
	if (edges.next(line, length))
		header.assign(line, length);
	vector<string> tokens = split(header, ';');
	// Find the idx of the targets and sources in the header
	int srcIdx = -1, tgtIdx = -1;
	for (int idx=0 ; idx<(int) tokens.size() ; idx++) {
		// To lower case
		transform(tokens[idx].begin(), tokens[idx].end(), tokens[idx].begin(), ::tolower);

		if (tokens[idx] == "source") {
			srcIdx = idx;
		} else if (tokens[idx] == "target") {
			tgtIdx = idx;
		}
	}

	if (srcIdx == -1 || tgtIdx == -1) {
		cerr << "Source or Target label absent from the edges csv file" << endl;
		return g;
	}

	// Degrees
	int idx1, idx2;
	while (edges.next(line, length)) {
		if (length > 0) {
			edgeEnds(line, length, srcIdx, tgtIdx, readsNb, idx1, idx2);
			g.offsets[idx1 + 1]++;
			g.offsets[idx2 + 1]++;
		}
	}
	for (int i = 0 ; i < readsNb ; i++)
		g.offsets[i + 1] += g.offsets[i];

	// Neighbors, in the order of the file
	g.targets.resize(g.offsets[readsNb]);
	vector<uint64_t> next (g.offsets.begin(), g.offsets.end() - 1);
	edges.rewind();
	edges.next(line, length);
	while (edges.next(line, length)) {
		if (length > 0) {
			edgeEnds(line, length, srcIdx, tgtIdx, readsNb, idx1, idx2);
			if (next[idx1] == g.offsets[idx1 + 1] || next[idx2] == g.offsets[idx2 + 1]) {
				cerr << "The edges csv file changed while it was read" << endl;
				exit(1);
			}
			g.targets[next[idx1]++] = idx2;
			g.targets[next[idx2]++] = idx1;
		}
	}

	return g;
}


void save_componants (MetaGraph & metagraph, int readsNb, ifstream & nodes, ofstream & componentsStream) {
	// Components are grown from the meta node of smallest idx not yet
	// assigned, following the neighbor lists (which are not always symmetric
	// at this point). The set reached does not depend on the order of the
	// visit; nodes are sorted by idx, so the smallest idx comes first.
	int compIdx = 0;
	vector<int> components (readsNb, -1);
	vector<int> metaNodes (readsNb, -1);
	vector<char> toCompute (metagraph.nodes.size(), 1);
	vector<int> currentComponent;
	for (size_t start = 0 ; start < metagraph.nodes.size() ; start++) {
		if (!toCompute[start])
			continue;
		toCompute[start] = 0;
		currentComponent.push_back(start);

		while (currentComponent.size() > 0) {
			MetaNode & mn = metagraph.nodes[currentComponent.back()];
			currentComponent.pop_back();

			// Perform assignations
			for (int read : mn.subNodes) {
				components[read] = compIdx;
				metaNodes[read] = mn.idx;
			}

			// Add neighbors
			for (int neiIdx : mn.neighbors) {
				int position = metagraph.positionOf(neiIdx);
				if (position >= 0 && toCompute[position]) {
					toCompute[position] = 0;
					currentComponent.push_back(position);
				}
			}
		}

		compIdx++;
	}

	// Save components
	string line;
	getline(nodes, line);
	componentsStream << line << ";metanode;component" << '\n';

	for (int i=0 ; i<readsNb ; i++) {
		getline(nodes, line);
		componentsStream << line << ";" << metaNodes[i] << ";" << components[i] << '\n';
	}
}


void save_metagraph (const MetaGraph & graph, ofstream & nodes, ofstream & edges) {
	nodes << "node;size" << '\n';
	edges << "from;to" << '\n';

	for (const MetaNode & mn : graph.nodes) {
		nodes << mn.idx << ";" << mn.subNodes.size() << '\n';

		for (int neiIdx : mn.neighbors) {
			if (neiIdx > mn.idx)
				edges << mn.idx << ";" << neiIdx << '\n';
		}
	}
}
