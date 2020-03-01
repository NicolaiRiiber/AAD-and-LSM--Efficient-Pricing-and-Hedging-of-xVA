#pragma once
#include <vector>
#include <math.h>

extern Tape *GlobalTape;

/* 
 * AADNumber.h
 * Custom number class responsible
 * for the actual Adjoint Automatic
 * Differentiation backprobagation.
 * The Number class contains the derivative
 * library built with operator overloading 
 * and a calculateAdjoints method that performs 
 * the backprobagation based on the operations
 * that are recorded on the tape.
 */
class Number
{
	
public:

	// Default constructor that does nothing.
	Number() {}

	// Constructs with an index.
	Number(double val) : tape(GlobalTape), value(val) 
	{
		idx = tape->tapeIdx();
	}

	// Constructs with an index and a value.
	Number(double val, size_t index) : tape(GlobalTape), value(val), idx(index) {}

	Number& operator+=(const Number& rhs) {
		idx = tape->tapeIdx(idx, 1.0, rhs.idx, 1.0);
		value += rhs.value;
		return *this;
	}

	Number operator+(const Number& rhs) {
		return Number(*this) += rhs;
	}
	
	Number& operator-=(const Number& rhs) {
		idx = tape->tapeIdx(idx, 1.0, rhs.idx, -1.0);
		value -= rhs.value;
		return *this;
	}

	Number operator-(const Number& rhs) {
		return Number(*this) -= rhs;
	}

	Number& operator*=(const Number& rhs) {
		idx = tape->tapeIdx(idx, rhs.value, rhs.idx, value);
		value *= rhs.value;
		return *this;
	}

	Number operator*(const Number& rhs) {
		return Number(*this) *= rhs;
	}

	Number& operator/=(const Number& rhs) {
		idx = tape->tapeIdx(idx, 1.0/rhs.value, rhs.idx, -value/(rhs.value*rhs.value));
		value /= rhs.value;
		return *this;
	}

	 Number operator/(const Number& rhs) {
		return Number(*this) /= rhs;
	}

	 std::vector<double> calculateAdjoints() {

		 // Initialize the backprobagation.
		 size_t N = tape->tape.size();            // Find N.
		 std::vector<double> adjoints(N, 0.0);    // Initialize all adjoints to 0.
		 adjoints[N - 1] = 1.0;                   // Seed the backprobagation algorithm.

		 // Backprobagation (recursive loop).
		 for (int i = N - 1; i > 0; i--)
		 {
			 // Probagate first argument.
			 adjoints[tape->tape[i].idx[0]] += adjoints[i] * tape->tape[i].der[0];

			 // Probagate second argument.
			 adjoints[tape->tape[i].idx[1]] += adjoints[i] * tape->tape[i].der[1];
		 }

		 // Return the adjoints.
		 return adjoints;
	 }

	 // Public class member variables
	 Tape* tape;
	 size_t idx;
	 double value;
	
}; 

inline Number operator+(Number& lhs, double rhs) {
	return Number(lhs.value + rhs, lhs.tape->tapeIdx(lhs.idx, 1.0));
}

inline Number operator+(double lhs, Number& rhs) {
	return Number(rhs.value + lhs, rhs.tape->tapeIdx(rhs.idx, 1.0));
}

inline Number operator-(Number& var) {
	return Number(-var.value, var.tape->tapeIdx(var.idx, -1.0));
}

inline Number operator-(Number& lhs, double rhs) {
	return Number(lhs.value - rhs, lhs.tape->tapeIdx(lhs.idx, 1.0));
}

inline Number operator-(double lhs, Number& rhs) {
	return Number(lhs - rhs.value, rhs.tape->tapeIdx(rhs.idx, -1.0));
}

inline Number operator*(Number& lhs, double rhs) {
	return Number(lhs.value * rhs, lhs.tape->tapeIdx(lhs.idx, rhs));
}

inline Number operator*(double lhs, Number& rhs) {
	return Number(rhs.value * lhs, rhs.tape->tapeIdx(rhs.idx, lhs));
}

inline Number operator/(Number& lhs, double rhs) {
	return lhs * (1.0 / rhs);
}

inline Number operator/(double lhs, Number& rhs) {
	return Number(lhs / rhs.value, rhs.tape->tapeIdx(rhs.idx, -lhs / (rhs.value*rhs.value)));
}

inline bool operator>(const Number& lhs, const Number& rhs) {
	return lhs.value > rhs.value ? true : false;
}

inline bool operator>=(const Number& lhs, const Number& rhs) {
	return lhs.value >= rhs.value ? true : false;
}

inline bool operator<(const Number& lhs, const  Number& rhs) {
	return lhs.value < rhs.value ? true : false;
}

inline bool operator<=(const Number& lhs, const  Number& rhs) {
	return lhs.value <= rhs.value ? true : false;
}

inline bool operator>(const Number& lhs, const double& rhs) {
	return lhs.value > rhs ? true : false;
}

inline bool operator>=(const Number& lhs, const double& rhs) {
	return lhs.value >= rhs ? true : false;
}

inline bool operator<(const Number& lhs, const double& rhs) {
	return lhs.value < rhs ? true : false;
}

inline bool operator<=(const Number& lhs, const double& rhs) {
	return lhs.value <= rhs ? true : false;
}

inline bool operator>(const double& lhs, const  Number& rhs) {
	return lhs > rhs.value ? true : false;
}

inline bool operator>=(const double& lhs, const Number& rhs) {
	return lhs >= rhs.value ? true : false;
}

inline bool operator<(const double& lhs, const Number& rhs) {
	return lhs < rhs.value ? true : false;
}

inline bool operator==(const Number& lhs, const Number& rhs) {
	return abs(lhs.value - rhs.value) < 0.000001 ? true : false;
}

inline bool operator==(const Number& lhs, const double& rhs) {
	return abs(lhs.value - rhs) < 0.000001 ? true : false;
}

inline bool operator==(const double& lhs, const Number& rhs) {
	return (rhs == lhs);
}

inline bool operator!=(const Number& lhs, const Number& rhs) {
	return abs(lhs.value - rhs.value) < 0.000001 ? false : true;
}

inline bool operator!=(const Number& lhs, const double& rhs) {
	return abs(lhs.value - rhs) < 0.000001 ? false : true;
}

inline bool operator!=(const double& lhs, const Number& rhs) {
	return (rhs != lhs);
}

inline Number exp(Number self) {
	return Number(exp(self.value),
		(self.tape)->tapeIdx(self.idx, exp(self.value)));
}

inline Number sqrt(Number self) {
	return Number(sqrt(self.value),
		(self.tape)->tapeIdx(self.idx, 1 / (2 * sqrt(self.value))));
}

inline Number abs(Number self) {
	return Number(abs(self.value),
		(self.tape)->tapeIdx(self.idx, self.value > 0 ? 1 : -1));
}

inline Number log(Number self) {
	return Number(log(self.value),
		(self.tape)->tapeIdx(self.idx, 1 / (self.value)));
}

double value(double& k) {
	return k;
}

double value(Number& k) {
	return k.value;
}

template<typename T>
double value(T& const k) {
	return value(k);
}

template<typename T>
double as_double(T& k) {
	return value(k);
}

double as_double(double k) {
	return k;
}

double as_double(Number k) {
	return k.value;
}

