/*
MTH 4135 Spring 2017
Group Members: Jaime Abbariao, Nadime Uddin, Abdul Haque, Michelle Tam, Xuming Shi
*/
#include "Declarations.h"
// type gcc TimedPi.cpp into command prompt to compile the code.
int main() {

   unsigned int i, n;
   double X, Y, V, W, Z, I, S, Ibar, elapsedtime;

   // Get the number of simulations you want to perform
   printf ("How many simulations should I do?... ");
   n = GetUnsignedInteger ();
   
   // We use an arbitrary seed. In this case, I stick with 25.
   MTUniform(25);

   // Start the clock
   Time ();

   S = 0;

   // we loop according to the number of simulations entered
   for (i = 1; i <= n; i++) {

      /*
         MTUniform(0) runs through the values of the generated uniform random variable on (0, 1)
         the 0 means that it's moving from value to value such that all values of MTUniform(0) used
         is distinct.
      */


      // We have to multiply this value by 2 and subtract by 1 because we need the point to to be 
      // between all values from (-1, 1)

      X = 2*MTUniform(0) - 1;

      Y = 2*MTUniform(0) - 1;

      V = 2*MTUniform(0) - 1;

      W = 2*MTUniform(0) - 1;

      Z = 2*MTUniform(0) - 1;

      // If this condition passes, then the point is on the ball
      // We increment S by 1 if the point is in the ball otherwise we increment by 0
      // S in this case works like an indicator random variable collecting binary values
      S += (X*X + Y*Y + V*V + W*W + Z*Z <= 1.0 ?  1 : 0);

   }

   // We find the probability of the event {the point is in the ball} happening
   // S denotes the number of times the event happened and n denotes the total number of simulations
   Ibar = S / n;

   // Time stops here
   elapsedtime = Time ();

   // Results:
   // Note that since we know the volume of the 5 dimensional cube C is 32
   // and that B is a subset of C, then we can estimate B as B = Vol(C)*Ibar
   printf ("The estimated volume of the 5-dimensional ball is %8.6f\n", 32*Ibar);

   // Wikipedia says that the numerical value of the 5 dimensional ball is
   // (8(Pi^2) / 15)* R^5 where R in this case is 1 as we have a unit ball
   // M_PI is defined within the math.h as a constant for PI
   printf ("The numerical volume of the 5-dimensional ball is %8.6f\n", 8 * (M_PI*M_PI) / 15);
   printf ("%d simulations took %lf seconds.\n", n, elapsedtime);
   
   
   Exit ();

}

#include "Definitions.h"
