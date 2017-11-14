#include <iostream>
#include <nfl/poly_p.hpp>
#include "tools.h"

template <size_t degree, size_t modulus, class T>
bool run() {
	using poly_p = nfl::poly_p_from_modulus<T, degree, modulus>;
	using poly_t = nfl::poly_from_modulus<T, degree, modulus>;

	bool ret = true;
	poly_t x{1,2,3};
	poly_t y{1};
	poly_t z = x * y; // I'd expect z to be {1,2,3,0...}
	// but it's {1,0,0...}
	
	ret &= (x==z);
	assert(x==z);
	
	std::cout << "x=" << x << std::endl;
	std::cout << "z=" << z << std::endl;
	
	// Seeing that x={1,2,3} and z = {1,0,0}
	// I'd expect that (x==z)=FALSE but (x==z)=TRUE !!!
	
	nfl::mul(z,x,y); // This is going to have the same strange behaviour

	ret &= (x==z);
	assert(x==z);
	
	std::cout << "x=" << x << std::endl;
	std::cout << "z=" << z << std::endl;
	
	// The equality test is completely bugged
	x={1,2,3,4,5,6,7,8};
	y={0};
	assert(x==y); // returns TRUE
	
	return ret;
}

int main(int argc, char const* argv[]) {
	return not run<CONFIG>();
}

