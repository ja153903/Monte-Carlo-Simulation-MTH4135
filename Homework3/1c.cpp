#include "Declarations.h"
#include "Definitions.h"
#include "Portfolio.h"

int main(){

	double r = 0.05, sigma = 0.3, S0 = 100, T = 0.5, N = 200, P0 = 1000.0,
	 B, pbar, p2bar, elapsed_time, t_star, stdhatP, error, epsilon=0.1,
	done, S, antiS, StdNormal, money, dB, t, U;
	int n;
	double mu = r - 0.5*sigma*sigma;
	double dt = T/N;
	
	// Seed the RNG
	MTUniform(1);

	// Print column headings for output to execution window.
	printf ("         n       Port        +/-        t       t*\n");

	// Initialize counter variables
	done = n = pbar = p2bar = money = 0;

	// Start the timer
	Time();

	// initialize regular and antithetic portfolio
	Portfolio p, antiP;

	while (!done){

		// More values to initialize
		S = antiS = S0;
		B = 0;
		t = 0;

		// Initialize the portfolio to a money market investment
		p.funds = P0;
		p.shares = 0;
		p.mma = 1;

		// Initialize the antithetic portfolio to a money market investment
		antiP.funds = P0;
		antiP.shares = 0;
		antiP.mma = 1;

		for (int i = 1; i <= N; i++){

			// current time
			t = i * dt;

			// generate random uniform
			U = MTUniform(0);

			// generate standard normal
			StdNormal = PsiInv(U);

			// Brownian step
			dB = sqrt(dt) * StdNormal;

			// Brownian motion at that step
			B += dB;

			// update stock price
			S = S0 * exp(mu * t + sigma * B);

			// update antithetic stock price
			antiS = S0 * exp(mu * t + sigma * -B);
			
			// implement strategy
			strategy(p, S, dt, T, r, 'c');
			strategy(antiP, antiS, dt, T, r, 'c');
		}

		// close positions
		closePositions(p, S);
		closePositions(antiP, antiS);

		// find the average of antithetic and regular variable
		money = (p.funds + antiP.funds) / 2;

		// Compute first and second moment
		++n;
		pbar = (pbar * (n-1) + money) / n;
		p2bar = (p2bar * (n-1) + money * money) / n;

		// Error test
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