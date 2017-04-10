#include "Declarations.h"
#include <iostream>
#include <iomanip>
// Compile with g++ -o elec ElectionSim.cpp
// Run with ./elec

int main () {


   int i;
   double E;
   int ElectoralVotes[52];
   double ClintonProbability[52];

   ElectoralVotes[ 1] =     9; ClintonProbability[ 1] =   0.1;  // Alabama
   ElectoralVotes[ 2] =     3; ClintonProbability[ 2] =  23.5;  // Alaska
   ElectoralVotes[ 3] =    11; ClintonProbability[ 3] =  33.4;  // Arizona
   ElectoralVotes[ 4] =     6; ClintonProbability[ 4] =   0.4;  // Arkansas
   ElectoralVotes[ 5] =    55; ClintonProbability[ 5] = 100.0;  // California
   ElectoralVotes[ 6] =     9; ClintonProbability[ 6] =  75.6;  // Colorado
   ElectoralVotes[ 7] =     7; ClintonProbability[ 7] =  97.3;  // Connecticut
   ElectoralVotes[ 8] =     3; ClintonProbability[ 8] =  91.5;  // Delaware
   ElectoralVotes[ 9] =     3; ClintonProbability[ 9] = 100.0;  // D.C.
   ElectoralVotes[10] =    29; ClintonProbability[10] =  55.1;  // Florida
   ElectoralVotes[11] =    16; ClintonProbability[11] =  20.9;  // Georgia
   ElectoralVotes[12] =     4; ClintonProbability[12] =  98.9;  // Hawaii
   ElectoralVotes[13] =     4; ClintonProbability[13] =   0.9;  // Idaho
   ElectoralVotes[14] =    20; ClintonProbability[14] =  98.3;  // Illinois
   ElectoralVotes[15] =    11; ClintonProbability[15] =   2.5;  // Indiana
   ElectoralVotes[16] =     6; ClintonProbability[16] =  30.2;  // Iowa
   ElectoralVotes[17] =     6; ClintonProbability[17] =   2.6;  // Kansas
   ElectoralVotes[18] =     8; ClintonProbability[18] =   0.4;  // Kentucky
   ElectoralVotes[19] =     8; ClintonProbability[19] =   0.5;  // Louisiana
   ElectoralVotes[20] =     4; ClintonProbability[20] =  82.6;  // Maine
   ElectoralVotes[21] =    10; ClintonProbability[21] = 100.0;  // Maryland
   ElectoralVotes[22] =    11; ClintonProbability[22] = 100.0;  // Massachusetts
   ElectoralVotes[23] =    16; ClintonProbability[23] =  78.9;  // Michigan
   ElectoralVotes[24] =    10; ClintonProbability[24] =  85.0;  // Minnesota
   ElectoralVotes[25] =     6; ClintonProbability[25] =   2.2;  // Mississippi
   ElectoralVotes[26] =    10; ClintonProbability[26] =   3.9;  // Missouri
   ElectoralVotes[27] =     3; ClintonProbability[27] =   4.1;  // Montana
   ElectoralVotes[28] =     5; ClintonProbability[28] =   2.3;  // Nebraska
   ElectoralVotes[29] =     6; ClintonProbability[29] =  58.3;  // Nevada
   ElectoralVotes[30] =     4; ClintonProbability[30] =  69.8;  // New Hampshire
   ElectoralVotes[31] =    14; ClintonProbability[31] =  96.9;  // New Jersey
   ElectoralVotes[32] =     5; ClintonProbability[32] =  82.6;  // New Mexico
   ElectoralVotes[33] =    29; ClintonProbability[33] =  99.8;  // New York
   ElectoralVotes[34] =    15; ClintonProbability[34] =  55.5;  // North Carolina
   ElectoralVotes[35] =     3; ClintonProbability[35] =   2.3;  // North Dakota
   ElectoralVotes[36] =    18; ClintonProbability[36] =  35.4;  // Ohio
   ElectoralVotes[37] =     7; ClintonProbability[37] =   0.0;  // Oklahoma
   ElectoralVotes[38] =     7; ClintonProbability[38] =  93.7;  // Oregon
   ElectoralVotes[39] =    20; ClintonProbability[39] =  77.0;  // Pennsylvania
   ElectoralVotes[40] =     4; ClintonProbability[40] =  93.2;  // Rhode Island
   ElectoralVotes[41] =     9; ClintonProbability[41] =  10.3;  // South Carolina
   ElectoralVotes[42] =     3; ClintonProbability[42] =   6.1;  // South Dakota
   ElectoralVotes[43] =    11; ClintonProbability[43] =   2.7;  // Tennessee
   ElectoralVotes[44] =    38; ClintonProbability[44] =   6.0;  // Texas
   ElectoralVotes[45] =     6; ClintonProbability[45] =   3.3;  // Utah
   ElectoralVotes[46] =     3; ClintonProbability[46] =  98.1;  // Vermont
   ElectoralVotes[47] =    13; ClintonProbability[47] =  85.5;  // Virginia
   ElectoralVotes[48] =    12; ClintonProbability[48] =  98.4;  // Washington
   ElectoralVotes[49] =     5; ClintonProbability[49] =   0.3;  // West Virginia
   ElectoralVotes[50] =    10; ClintonProbability[50] =  83.5;  // Wisconsin
   ElectoralVotes[51] =     3; ClintonProbability[51] =   1.1;  // Wyoming

   int done = 0; // Flag variable to end simulation
   int n = 0; // Represents the number simulations
   double success = 0, antiSuccess = 0; // Represents the number of times the simulation votes is less than 233
   int actualClintonVotes = 233; // Actual Clinton Votes from election
   double standardDeviation, error, epsilon = 0.001, elapsedTime, timeLeft, expectedVotesFirstMoment = 0, 
   expectedVotesSecondMoment = 0, C = 0, U, antiC = 0, avr;

   MTUniform(41); // Seed the RNG

   Time();

   double rho; // the parameter that governs the correlation of the model
   std::cout << "Enter in your rho (Either 0 or 1): ";
   std::cin >> rho;

   std::cout << "Simulations " << " Votes " << " Error " << " T " << " T_star" << std::endl;

   
   int clintonVotesForWin = 270; // number of votes needed to win the presidency

   while (!done){

      C = 0;
      antiC = 0;

      success = 0;
      antiSuccess = 0;

      // If rho is 1, generate 1 fresh uniform to use throughout the path
      if (rho == 1){
         U = MTUniform(0);
         for (int j = 1; j <= 51; j++){

            if (U <= ClintonProbability[j]/100.0){
               C = C + ElectoralVotes[j];
            }

            if (1- U <= ClintonProbability[j]/100.0){
               antiC = antiC + ElectoralVotes[j];
            }

         }
      }

      // If rho is 0, generate fresh uniforms every iteration of the path
      if (rho == 0){
         for (int j = 1; j <= 51; j++){
            U = MTUniform(0);

            if (U <= ClintonProbability[j]/100.0){
               C = C + ElectoralVotes[j];
            }

            if (1- U <= ClintonProbability[j]/100.0){
               antiC = antiC + ElectoralVotes[j];
            }

         }
      }

      // Increment if Clinton becomes president
      if (C >= clintonVotesForWin){ ++success; } 
      if (antiC >= clintonVotesForWin){ ++antiSuccess; }

      avr = (antiSuccess + success) * 0.5; // antithetic variance reduction

      ++n; // Increment the number of simulations

      // Sample first moment
      expectedVotesFirstMoment = (expectedVotesFirstMoment * (n - 1) + avr) / n;

      // Sample second moment
      expectedVotesSecondMoment = (expectedVotesSecondMoment * (n - 1) + avr * avr) / n;

      // Error testing
      if (n % 100000 == 0){

         standardDeviation = sqrt(expectedVotesSecondMoment - expectedVotesFirstMoment * expectedVotesFirstMoment);

         error = 1.96 * standardDeviation / sqrt(n); // 95% confidence

         elapsedTime = Time();

         timeLeft = elapsedTime * pow(error / epsilon, 2);

         std::cout << std::setprecision(8) << n << " " << expectedVotesFirstMoment << " "
         << error << " " << elapsedTime << " " << timeLeft << std::endl;

         if (error < epsilon) { done = 1; }

      }

      
   }

   
   Exit();

}





#include "Definitions.h"