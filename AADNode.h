#pragma once
#include <vector>
#include <math.h>

/* Abstract class implementation
 * of a node in a computation tree 
 * used for recording the Tape in 
 * Automatic Adjoint Differentiation.
 *
 * @Input: double der1: partial derivative to the firt argument
 *         double der2: partial derivative to the second argument
 *         size_t idx1: index of the first argument on the tape
 *         size_t idx2: index of the second argument on the tape
 */
class Node
{
	
public:

	// Constructor
	Node(double der1, double der2, size_t idx1, size_t idx2) : der{ der1, der2 }, idx{ idx1, idx2 } {}

	// Public class member variables
	double der[2];
	size_t idx[2];
};

