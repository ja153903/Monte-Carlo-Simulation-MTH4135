#include "Declarations.h"


// This program simulates many ICS paths of the form ICS = mu + sigma * U,
// The ICS process starts at its long term mean of mu = 84.971.


int main () {

  double dt, U, antiU, ICS, alpha=0.050, u, mu=84.971, sigma=4.055, mu0, sigma0,
      n, tbar, t2bar, error, stdhatt, elapsed_time, t_star;

    int i, t, test, done;
    

   // In this application Delta t is 1 (month).
   dt = 1.0;

   // Specify the 95% error tolerance: one tenth of a month.
   double epsilon = 0.1;

   // Seed the RNG.
   MTUniform (1);

   // Print column headings for output to execution window.
   printf ("         n       tbar        +/-        t       t*\n");

   // Initialize certain values.
   t2bar = tbar = done = n = test = 0;

   // Conditional standard deviation.
   sigma0 = sqrt ( (1.0 - exp (-2.0*alpha*dt)) / (2.0*alpha) );

   // Start the clock to time the computations.
   Time ();

   // Begin the simulation loop 
   while (!done) {

     // Generate the O-U path, starting at...
     ICS = 85;

     // Compute the corresponding value of the underlying O-U process.
     U = (ICS - mu) / sigma;
      
     // Initialize time (in months).
     t = 0;

    // March forward in time to get each month's ICS value.
      for(i=0; i<120; ++i) {

       // Compute the mean and variance of U_(t+dt) given that U_t = u.
       u = U;

       // Conditional mean.
       mu0 = u * exp (-alpha * dt);

       // Now generate a normal RV with this mean and standard deviation.
       U = mu0 + sigma0 * PsiInv(MTUniform(0));
       antiU = mu0 + sigma0 * -PsiInv(MTUniform(0));    
       ICS = mu + sigma * ((U+antiU)/2);          

       if(ICS>=95)
        ++t;


       }
  
    // Update the error test counter, simulation counter, and sample moments of t.
       ++test;
       ++n;
       tbar = ((n-1)*tbar + t) / n;
       t2bar  = ((n-1)*t2bar + t*t) / n;

      // Test the error condition periodically.
      if (test % 100000 == 0) {

        // Estimate the standard deviation of S.
        stdhatt  = sqrt (t2bar - tbar*tbar);

        // Estimate the error of the Vbar estimator.
        error = 1.96 * stdhatt / sqrt(n);

        // Compute the elapsed time and estimate the time of completion.
        elapsed_time = Time ();
        t_star = elapsed_time * pow (error / epsilon, 2);

        // Report.
        printf ("%10.0f   %8.4f   %8.6f %8.3f %8.3f\n", n, tbar, error, elapsed_time, t_star);

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