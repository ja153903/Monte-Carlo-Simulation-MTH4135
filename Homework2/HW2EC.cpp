#include "Declarations.h"

double blackScholesOptionPricing(double S0, double K, double r, double sigma,
   double T){
   // d1 and d2 for Black Scholes Call Option Formula
   double d1 = (log(S0/K) + (r+sigma*sigma/2.)*(T))/(sigma*sqrt(T));
   double d2 = (log(S0/K) + (r-sigma*sigma/2.)*(T))/(sigma*sqrt(T));
   // Call Option Price via Black-Scholes Pricing Model
   return S0 * Psi(d1) - K * exp(-r * T) * Psi(d2);
}
double max(double a, double b){
   return (a > b ? a : b);
}
int main() {
   double t, V,  elapsed_time, t_star, stdhatV, error, U, antiS;
   double **S;
   int **Shout;
   int W;
   // Time to expiration.
   double T = 0.5;
   // Number of stock price periods.
   int N = 50;
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
   // Allocate space for the binomial tree
   // initialize space of multidim array
   S = (double **) calloc (N+1, sizeof (double *));
   // Loop to initialize space for each price path
   for (int i = 0; i <= N; i++){
      // initializes path with 0
      S[i] = (double *) calloc (2*i + 1, sizeof (double));
      S[i] = S[i] +1;
   }
   // Allocate space for the Shout indicator
   // initialize space for shout strategy
   Shout = (int **) calloc (N+1, sizeof (int *));
   // loop to initialize space for each price path shout strategy
   for (int i = 0; i <= N; i++){
      // initialize path with 0
      Shout[i] = (int *) calloc (2 *i + 1, sizeof (int));
      Shout[i] = Shout[i] + 1;
   }
   double value;
   // Create binomial tree and populate strategy array
   for (int i = 0; i <= N; i++){
      t = i * dt; // advance time
      // loop takes care of all possible permutations.
      for (int j = -i; j <= i; j += 2){
         S[i][j] = S0 * exp(mu * t + sigma * sqrtdt * j);
         value = blackScholesOptionPricing(S0, K, r, sigma, T-t);
         if (S[i][j] - K > value){
            // Fill shout multidimensional array with 1
            // because this is when we wish to chat
            Shout[i][j] = 1;
         }
      }
   }
   // Used for estimating a
   double varA = 0;
   double covCA = 0;
   double callPO = 0;
   double antiCallPO = 0;
   double A = 0;
   // Initialize certain values.
   double V2bar = 0;
   double Vbar = 0;
   int done = 0;
   double n = 0;
   int test = 0;
   double shout, antiShout;
   // Start the clock to time the computations.
   Time ();
   // Seed the RNG.
   MTUniform (1);
   // Estimate constant for control variable
   for (int i = 1; i <= 100; i++){
      // Set to -1 if haven't shouted yet.
      shout = antiShout = -1;
      // Assign 0 to W
      W = 0;
      for (int j = 1; j <= N; j++){
         // Generate random variable
         U = MTUniform(0);
         // Advance the random walk
         W = W + (U <= 0.5 ? -1 : 1);
         // You can only shout once if
         // nothing has been shouted and
         // if strategy tells you to
         if (Shout[j][W] && shout < 0){shout = S[j][W];}
         // You can only shout once if
         // nothing has been shouted and
         // if strategy tells you to
         if (Shout[j][-W] && antiShout < 0){antiShout = S[j][-W];}

      }
      // Compute call option payoff
      callPO = max(max(S[N][W]-K, shout - K), 0);
      // Compute the antithetic call option payoff
      antiCallPO = max(max(S[N][-W]-K, antiShout -K), 0);
      // Create control variable
      A = W*W - N;
      // Compute sample variance of A
      varA = ((i-1)*varA + A*A)/i;
      // Compute sample covariance of C and A
      covCA = ((i-1)*covCA + ((antiCallPO + callPO)*0.5)*A)/i;
   }
   double a = -(covCA/varA);
   // Print column headings for output to execution window.
   printf ("         n       Vbar        +/-        t       t*\n");
   // Initialize these values to 0
   callPO = antiCallPO = 0;
   // Begin the simulation loop.
   while (!done) {
      shout = antiShout = -1;
      // Initialize the random walk.
      W = 0;
      // Simulate a stock price path.  Go forward period-by-period computing
      //   the next stock price.
      for (int i = 1; i <= N; i++) {
         // Advance the random walk.
         U = MTUniform (0);
         W = W + (U <= 0.5 ? -1 : +1);
         if (Shout[i][W] && shout < 0){ shout = S[i][W];}
         if (Shout[i][-W] && antiShout < 0){ antiShout = S[i][-W];}
      }
      // Compute call option payoff
      callPO = max(max(S[N][W]-K, shout - K), 0);
      // Compute the antithetic call option payoff
      antiCallPO = max(max(S[N][-W]-K, antiShout -K), 0);
      // Create control variable
      A = W*W - N;
      // Discount back to time 0.
      V = Discount_factor * (((callPO+antiCallPO)/2) + a*A);
      // Update the simulation counter and the sample moments of V at time T.
      ++n;
      Vbar  = ((n-1)*Vbar + V) / n;
      V2bar = ((n-1)*V2bar + V*V) / n;
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
         printf ("%10.0f   %8.4f   %8.6f %8.3f %8.3f\n", 
            n, Vbar, error, elapsed_time, t_star);
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