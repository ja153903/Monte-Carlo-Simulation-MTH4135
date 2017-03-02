#include "Declarations.h"
int main() {
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
   // d1 and d2 for Black Scholes Call Option Formula
   double d1 = (log(S0/K) + (r+sigma*sigma/2.)*(T))/(sigma*sqrt(T));
   double d2 = (log(S0/K) + (r-sigma*sigma/2.)*(T))/(sigma*sqrt(T));
   // Call Option Price via Black-Scholes Pricing Model
   double blackScholesCallOption = S0 * Psi(d1) - K * Discount_factor * Psi(d2);
   // Print out call option price
   printf("The actual price of the call option is %8.4f\n", blackScholesCallOption);
   Exit ();
}
#include "Definitions.h"