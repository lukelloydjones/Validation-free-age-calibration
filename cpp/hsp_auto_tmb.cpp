#include <TMB.hpp>                                // Links in the TMB libraries
#include <cstdlib>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(deltaY);                          // Difference in capture years
  DATA_VECTOR(deltaG);                          // Difference in Ben score
  DATA_VECTOR(hspsYN);                          // Vectors of yes/no for whether individuals are HSPs
    
  PARAMETER(theta);                             // Inverse to adult pop size
  PARAMETER(psi);                               // mort/fec rate
  PARAMETER(beta);                              // Golden - ratio for age
  PARAMETER(logsi);                             // parameter for intra-chort pairs


  Type epi_llk = 0.0;                           // TMB weird declare the "objective function" (neg. log. likelihood)

  // Initial declarations
  int comptot;
  comptot = deltaY.size();
   
  Type thetaT, chk, diffi;

  thetaT = (exp(theta) / (1 + exp(theta)));

  // std::cout << "Number of comparisons " << comptot << std::endl;
    
  vector<Type> hsp_prob(comptot);
  vector<Type> deltaA(comptot);
  Type si   = exp(logsi);
  Type psiT = exp(psi);

  deltaA = (1 / beta) * deltaG;
  // hsp_prob =  4 * thetaT * exp(-psi * abs(deltaY - deltaA));

  vector<Type> chk_vec(comptot);
 
  // -------------------------
  // Compute likelihood Part 1 
  // -------------------------
  
  // This needs a little fixing to be +- or absolute
  for (int i = 0; i < comptot; i++) {

    chk   = deltaY(i) - deltaG(i) / beta;
    if (chk < 0)
    {
      chk_vec(i) = si;
    } else {
      chk_vec(i) =  1;
    }
  }

  hsp_prob = 4 * chk_vec * thetaT * exp(-psiT * abs(deltaY - (1 / beta) * deltaG));

  // ==================
  // Log-likelihood 
  // ==================
 
	epi_llk -= sum(dpois(hspsYN, hsp_prob, true));


  REPORT(hspsYN);
  REPORT(epi_llk);
  REPORT(hsp_prob);

  ADREPORT(1/thetaT);
  ADREPORT(si);

  return epi_llk;
}