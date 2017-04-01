#include "Declarations.h"
#include "Definitions.h"
#include "Portfolio.h"

int main(){

	double r = 0.05, sigma = 0.3, S0 = 100, T = 0.5, N = 200, P0 = 1000.0,
	 B, pbar, p2bar, elapsed_time, t_star, stdhatP, error, epsilon=0.1,
	done, S, antiS, StdNormal, money;
	int n;
	double mu = r - 0.5*sigma*sigma;
	double dt = T/N;
	
	// Seed the RNG
	MTUniform(1);

	// Print column headings for output to execution window.
	printf ("         n       Port        +/-        t       t*\n");

	// Initialize counter variables
	done = n = pbar = p2bar = money = 0;

	// Start timer
	Time();

	// Initialize portfolio variables
	Portfolio p, antiP;

	while (!done){

		// Initialize stock price
		S = antiS = S0;

		// Initialize Brownian step
		B = 0;

		// Reset Portfolio
		clearAll(p, P0);

		// Reset Antithetic Portfolio
		clearAll(antiP, P0);

		// Generate a standard normal r.v.
		StdNormal = PsiInv(MTUniform(0));

		// Brownian motion
		B = sqrt(T) * StdNormal;

		// Implement strategy for portfolio
		strategy(p, S, T, T, r, 'b');

		// Implement strategy for antithetic portfolio
		strategy(antiP, antiS, T, T, r, 'b');

		// Update stock price
		S = S0 * exp (mu * T + sigma * B);

		// Update stock price
		antiS = S0 * exp (mu * T + sigma * -B);

		// Close your postion
		closePositions(p, S);

		// Close your antithetic position
		closePositions(antiP, antiS);

		// Find average
		money = (p.funds + antiP.funds) / 2;

		// Compute first and second moment
		++n;
		pbar = (pbar * (n-1) + money) / n;
		p2bar = (p2bar * (n-1) + money * money) / n;

		// Error testing every 100000 tests
		if (n % 100000 == 0){

			stdhatP = sqrt(p2bar - pbar * pbar);

			error = 1.96 * stdhatP/ sqrt(n);

			elapsed_time = Time();

			t_star = elapsed_time * pow (error / epsilon, 2);

			printf ("%10.0i   %8.4f   %8.6f %8.3f %8.3f\n", 
				n, pbar, error, elapsed_time, t_star);

			if (error < epsilon){
				done = 1;
			}

		}
	}

	Exit();
}