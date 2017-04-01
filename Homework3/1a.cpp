#include "Declarations.h"
#include "Definitions.h"
#include "Portfolio.h"

int main(){

	double r = 0.05, sigma = 0.3, S0 = 100, T = 0.5, N = 200, P0 = 1000.0, dt = T/N;

	// Initialize Portfolio struct from Portfolio.h
	Portfolio p;

	// The initial portfolio funds is 1000
	p.funds = P0;

	// Implement the strategy
	strategy(p, S0, T, T, r, 'a');
	
	printf("%8.4f", p.funds);
	
	Exit();
}