#pragma once
#include <vector>
#include <math.h>
#include "AADNode.h"

/* 
 * Abstract class that uses the
 * Node class to record each node
 * in the computation graph on
 * a tape.
 */
class Tape 
{
public:

	// Record a new node on the tape and return its index.
	size_t tapeIdx() {
		size_t idx = tape.size();
		tape.push_back(Node(0.0, 0.0, idx, idx));
		return tape.size();
	}

	// Record a new node on the tape and return its index.
	size_t tapeIdx(size_t idx1, double der1)	{
		size_t idx = tape.size();
		tape.push_back(Node(der1, 0.0, idx1, idx));
		return idx;
	}
	
	// Record a new node on the tape and return its index.
	size_t tapeIdx(size_t idx1, double der1, size_t idx2, double der2) {
		size_t idx = tape.size();
		tape.push_back(Node(der1, der2, idx1, idx2));
		return idx;
	}

	// Class member variables, the tape
	std::vector<Node> tape;

};

Tape* GlobalTape = new Tape;


