#include "Declarations.h"
#include "Definitions.h"
#include "Portfolio.h"

int main(){

	double r = 0.05, sigma = 0.3, S0 = 100, T = 0.5, N = 200, P0 = 0.0,
	 B, pbar, p2bar, elapsed_time, t_star, stdhatP, error, epsilon=0.1,
	done, S, antiS, StdNormal, money, dB, t;
	int n, test;
	double mu = r - 0.5*sigma*sigma;
	double dt = T/N;
	
	// Seed the RNG
	MTUniform(1);

	// Print column headings for output to execution window.
	printf ("         n       Port        +/-        t       t*\n");

	// Initialize the counter variables
	done = n = pbar = p2bar = money = test = 0;

	// Start the clock
	Time();

	// Initialize portfolios
	Portfolio p, antiP;

	while (!done){

		// More to initialize
		S = antiS = S0;
		B = 0;
		t = 0;

		// Initialize the portfolio to a money market investment
		clearAll(p, P0);

		// Initialize the antithetic portfolio to a money market investment
		clearAll(antiP, P0);

		for (int i = 1; i <= N; i++){

			// generate a standard normal
			StdNormal = PsiInv(MTUniform(0));

			// current time
			t = i * dt;

			// Brownian step
			dB = sqrt(dt) * StdNormal;

			// Brownian motion
			B += dB;

			// update the stock price
			S = S0 * exp(mu * t + sigma * B);

			// update the antithetic stock price
			antiS = S0 * exp(mu * t + sigma * -B);
			
			// implement the strategy
			strategy(p, S, dt, T-t, r, 'e');
			strategy(antiP, antiS, dt, T-t, r, 'e');
		}

		// close positions
		closePositions(p, S);
		closePositions(antiP, antiS);
		money = (p.funds + antiP.funds) / 2;

		// Compute the first and second moments
		++n;
		pbar = (pbar * (n-1) + money) / n;
		p2bar = (p2bar * (n-1) + money * money) / n;
		++test;

		// Error tests
		if (test % 100000 == 0){

			stdhatP = sqrt(p2bar - pbar * pbar);

			error = 1.96 * stdhatP/ sqrt(n);

			elapsed_time = Time();

			t_star = elapsed_time * pow (error / epsilon, 2);

			printf ("%10.0i   %8.4f   %8.6f %8.3f %8.3f\n", 
				n, pbar, error, elapsed_time, t_star);

			test = 0;

			if (error < epsilon){
				done = 1;
			}

		}
	}

	Exit();
}