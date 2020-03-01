#pragma once
#include <memory>
#include <vector>

template <class T, class U>
T modulo(T x, U y) {
	T r = x%y;
	return r < 0 ? r + y : r;
}

class mrg32k3a
{
	// Random generator object with the mrg32k3a methodology. 

public:

	mrg32k3a(long long s11, long long s12, long long s13, long long s21, long long s22, long long s23) {
		setseeds(s11, s12, s13, s21, s22, s23);
	};

	virtual std::unique_ptr<mrg32k3a> clone() const {
		std::unique_ptr<mrg32k3a> rng = std::unique_ptr<mrg32k3a>(new mrg32k3a(12345, 12345, 12345, 12346, 12346, 12346));
		return rng;
	}; //Makes a unique_ptr clone of the object, such that the threads don't "have their fingers on the same one"
	   //Thanks to raii semantics the use of unique_ptr's is easy and "cleans-up after itself"

	void setseeds(long long s11, long long s12, long long s13, long long s21, long long s22, long long s23) {
		x = { s11, s12, s13 };
		y = { s21, s22, s23 };
	};

	double next() 
	{
		//Outputs a single unif[0,1] from the mrg32k3a-stream of pseudo-random
		x.push_back(modulo(aOneTwo * x[1] + aOneThree * x[0], m1));
		y.push_back(modulo(aTwoOne * y[2] + aTwoThree * y[0], m2));
		x.erase(x.begin());
		y.erase(y.begin());

		return (double)(modulo(x[2] - y[2], m1)) / (double)m1;
	};

private:

	std::vector<long long> x, y; //state Number
	long long aOneTwo = 1403580;
	long long aOneThree = -810728;
	long long aTwoOne = 527612;
	long long aTwoThree = -1370589;
	long long m1 = 4294967087;
	long long m2 = 4294944443;

};
