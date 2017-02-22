#include "Declarations.h"

main() {

   double t, S, V, Vbar, V2bar, elapsed_time, t_star, stdhatV, error, n, U, W, antiS;
   int i, done, test;

   // Time to expiration.
   double T = 0.5;

   // Number of stock price periods.
   int N = 1000;

   // Time increment per period.
   double dt = T / N;

   // Compute this oft-used value once and for all.
   double sqrtdt = sqrt(dt);

   // Risk-free interest rate.
   double r = 0.05;

   // Compute the oft-used discount factor once and for all.
   double Discount_factor = exp(-r*T);

   // Stock price volatility.
   double sigma = .30;

   // Drift term.
   double mu = r - log(cosh(sigma*sqrtdt)) / dt;

   // Initial stock price.
   double S0 = 100.0;

   // Specify the 95% error tolerance.
   double epsilon = 0.005;

   // Initialize strike price
   double K = 110.0;

   // Used for estimating a
   double varA = 0;
   double covCA = 0;
   double callPO = 0;
   double antiCallPO = 0;
   double A = 0;

   // Start the clock to time the computations.
   Time ();

   // Seed the RNG.
   MTUniform (1);

   for (int i = 1; i <= 100; i++){

      // Assign S0 to S
      S = S0;

      // Assign 0 to W
      W = 0;

      // Assign 0 to t
      t = 0;

      for (int j = 1; j <= N; j++){

         // Generate random variable
         U = MTUniform(0);

         // Advance the random walk
         W = W + (U <= 0.5 ? -1 : 1);

         // Advance the random walk
         t = t + dt;

      }

      // Compute the stock price.
      S = S * exp (mu*t + sigma*sqrtdt*W);
      
      // Compute call option payoff
      callPO = (S > K) ? S-K : 0;

      // Compute the antithetic stock price.
      antiS = S * exp (mu*t + sigma*sqrtdt*-W);

      // Compute the antithetic call option payoff
      antiCallPO = (antiS >  K) ? antiS - K : 0;

      // Create control variable
      A = W*W - N;

      // Compute sample variance of A
      varA = ((i-1)*varA + A*A)/i;

      // Compute sample covariance of C and A
      covCA = ((i-1)*covCA + ((antiCallPO + callPO)*0.5)*A)/i;
   }

   // Print column headings for output to execution window.
   printf ("         n       Vbar        +/-        t       t*\n");

   // Initialize certain values.
   V2bar = Vbar = done = n = test = 0;

   // Initialize these values to 0
   callPO = antiCallPO = 0;

   // Begin the simulation loop.
   while (!done) {

      // Initialize the stock price.
      S = S0;

      // Initialize the random walk.
      W = 0;

      // Initialize time.
      t = 0;

      // Simulate a stock price path.  Go forward period-by-period computing
      //   the next stock price.
      for (i = 1; i <= N; i++) {

         // Advance the random walk.
         U = MTUniform (0);
         W = W + (U <= 0.5 ? -1 : +1);

         // Advance time by one period.
         t = t + dt;

      }

      // Compute the stock price.
      S = S0 * exp (mu*t + sigma*sqrtdt*W);
      
      // Compute call option payoff
      callPO = (S > K) ? S-K : 0;

      // Compute the antithetic stock price.
      antiS = S0 * exp (mu*t + sigma*sqrtdt*-W);

      // Compute the antithetic call option payoff
      antiCallPO = (antiS >  K) ? antiS - K : 0;

      // Create control variable
      A = W*W - N;

      // Discount back to time 0.
      V = Discount_factor * (((callPO+antiCallPO)/2) + (-covCA/varA)*A);

      // Update the simulation counter and the sample moments of V at time T.
      ++n;
      Vbar  = ((n-1)*Vbar + V) / n;
      V2bar = ((n-1)*V2bar + V*V) / n;

      // Reset the call option payoffs
      callPO = antiCallPO = 0;

      // Update the error condition test counter.
      ++test;

      // Test the error condition periodically.
      if (test == 100000) {

         // Estimate the standard deviation of S.
         stdhatV  = sqrt (V2bar - Vbar*Vbar);

          // Estimate the error of the Vbar estimator.
         error = 1.96 * stdhatV / sqrt(n);

         // Compute the elapsed time and estimate the time of completion.
         elapsed_time = Time ();
         t_star = elapsed_time * pow (error / epsilon, 2);

         // Report.
         printf ("%10.0f   %8.4f   %8.6f %8.3f %8.3f\n", n, Vbar, error, elapsed_time, t_star);

         // Reset the "test" counter and see if error tolerance is met.
         test = 0;
         if (error < epsilon) {
            done = 1;
         }

      }

   }

   Exit ();

}

#include "Definitions.h"