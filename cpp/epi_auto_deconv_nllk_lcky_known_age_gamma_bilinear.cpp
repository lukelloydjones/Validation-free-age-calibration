#include <TMB.hpp>                                // Links in the TMB libraries
// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Known age components

  DATA_VECTOR(kwn_a); // Known age vector
  DATA_VECTOR(knw_g);


  DATA_MATRIX(comb_yg);                         // Combinations matrix over Y and G
  DATA_MATRIX(kin_yg);                          // Kin matrix over Y and G
  DATA_IMATRIX(nz_inds);                        // Indices of the matrix to be iterated over

  // DATA_MATRIX(est_pr_a_samp_y);              // This is from deconvodisc and need to be incorporated later

  // Data 

  DATA_VECTOR(g_range);                         // Range to work over for epi-score, age, year samp, birth years
  DATA_VECTOR(a_range);                          
  DATA_VECTOR(y_range);                         
  DATA_VECTOR(b_range);                         
                                                                         

  DATA_SCALAR(fng_rate);                         // Intercept for the linear model
  DATA_INTEGER(phi_prior);                       // sd prior?
  DATA_SCALAR(Mu_phi);                           // Normal prior mean for phi if required
  DATA_SCALAR(SD_phi);                           // Normal prior variance for phi if required

    
  PARAMETER(N0_ad);                              // Abundance year naught
  PARAMETER(Z1);                                 // Female/or and male survival
  PARAMETER(R1);                                 // Female/or and male ROI
  PARAMETER(log_alpha);
  PARAMETER(beta);                               // Slope term to calibration epi-age to age
  PARAMETER(bp);                                 // Slope term to calibration after break point
  PARAMETER(theta);                              // Where does break in slopes happen?

  // -------------------------
  // Deconvolution necessities
  // -------------------------

  DATA_MATRIX(count);        // Observations (one dataset per column)
  DATA_MATRIX(Desmat_fix);   // Random effect design matrix
  DATA_MATRIX(Desmat_rand);  // Fixed effect design matrix
  DATA_MATRIX(chol_S);       // full-rank bit of S(moother)
  // DATA_MATRIX( M);        // convolution - this is done inside the code 
  PARAMETER_MATRIX(u);       // Random effects vectors (one dataset per col)
  PARAMETER_MATRIX(gammas);  // Fixed effects vectors
  PARAMETER(logsdu);         // Random effect standard deviations
  PARAMETER(log_phi);     // Variance of g-score given age
  PARAMETER(log_c);          // Variance of g-score given age

  // Setup the random effect (u)
  // ---------------------------

  Type epi_llk = 0.0;

  matrix<Type> u_indept = chol_S * u;
  epi_llk += dnorm( u_indept.vec(), Type(0), exp(logsdu), true).sum();


  // --------------------
  // Initial declarations
  // --------------------

  int yrs_tot, g_tot, a_tot, y_tot, byr_1, byr_2, year, gy_tot, n_combs; //, gy_tot, byrs_tot, syrs_tot, year, byr_1, byr_2, syrs_diff, yad, bju, rel_yad, sib_mt_it;
  Type Z1_t, exptd_g, pr_hsp_b1b2gy_a1a2, pr_hsap_b1b2gy, brk, c_t, shape, scale, alpha, phi;
  
  // Pre-exponentiate elements

  Z1_t  = exp(Z1);
  c_t   = exp(log_c);
  alpha = exp(log_alpha);
  phi   = exp(log_phi);


  // ------------------------
  // Known age log-likelihood
  // ------------------------

  vector<Type> ghat(kwn_a.rows());
  for (int i = 0; i < kwn_a.rows(); i++) 
  {
    // if (kwn_a(i) <= theta) {

    //   ghat(i) = alpha + beta * kwn_a(i);  // First segment: y = beta1 * x

    // } else {

    //   ghat(i) = alpha + beta * theta + bp * (kwn_a(i) - theta);  // Second segment: y = beta1 * cp + beta2 * x

    // }

    //Type diff = kwn_a(i) - theta;
    //ghat(i) = alpha + beta * CppAD::CondExpLe(diff, Type(0), diff, Type(0)) + bp * CppAD::CondExpGe(diff, Type(0), diff, Type(0))

    Type diff = kwn_a(i) - theta;
    ghat(i) = alpha - beta * CppAD::CondExpLe(-theta, Type(0), -theta, Type(0)) 
                    - bp * CppAD::CondExpGe(-theta, Type(0), -theta, Type(0))
                    + beta * CppAD::CondExpLe(diff, Type(0), diff, Type(0)) 
                    + bp * CppAD::CondExpGe(diff, Type(0), diff, Type(0));

    shape   = 1 / phi;
    scale   = phi * ghat(i);
    epi_llk += dgamma(knw_g(i), shape, scale, true);
  }


  // Get the range sizes

  yrs_tot = b_range.size();
  g_tot   = g_range.size();
  a_tot   = a_range.size();
  y_tot   = y_range.size();
  gy_tot  = comb_yg.rows();
  n_combs = nz_inds.rows();

  std::cout << "Years total " << yrs_tot << std::endl;
  std::cout << "G-Y combos "  << gy_tot << std::endl;
  std::cout << "a_tot "       << a_tot << std::endl;


  // Set up our main vectors and matrices

  vector<Type> N_ad(yrs_tot);
  matrix<Type> hsp_prob(yrs_tot, yrs_tot);
  matrix<Type> g_a_beta_theta(g_tot, a_tot);
  array<Type>  pr_ag_y(g_tot, y_tot, a_tot);
  array<Type>  pr_g_y(g_tot, y_tot);
  array<Type>  pr_a_gy(g_tot, y_tot, a_tot);
  array<Type>  all_prob_tab(gy_tot, gy_tot);


  // ---------------
  // Pop dynamics
  // ---------------

  // R1 = exp(R1); Change of scale for optimisation?
  N_ad[0] = exp(N0_ad);
  for (year = 1; year < yrs_tot; year++)
  {
    N_ad[year] = N_ad[0] * exp(R1 * year);
    // std::cout << "Adults year "  << year << " " << N_ad[year] << std::endl;
  }

  // std::cout << "N_ad size "  << N_ad.size() << std::endl;



  // Distribution of gs|a,beta,theta
  // These need to sum to 1 so make it a difference of pnorms 
  brk = (g_range(1) - g_range(0)) / 2;
  std::cout << "break" << brk << std::endl;
  Type minval = 1e-6;
  for (int i = 0; i < a_tot; i++)
  {
    // if (a_range[i] <= theta) {

    //   exptd_g = alpha + beta * a_range[i];  // First segment: y = beta1 * x

    // } else {

    //   exptd_g = alpha + beta * theta + bp * a_range[i];  // Second segment: y = beta1 * cp + beta2 * x

    // }

    Type diff = a_range[i] - theta;
    exptd_g = alpha - beta * CppAD::CondExpLe(-theta, Type(0), -theta, Type(0)) 
                    - bp * CppAD::CondExpGe(-theta, Type(0), -theta, Type(0))
                    + beta * CppAD::CondExpLe(diff, Type(0), diff, Type(0)) 
                    + bp * CppAD::CondExpGe(diff, Type(0), diff, Type(0));

    shape   = 1 / phi;
    scale   = phi * exptd_g;

    // First bin: Full CDF up to g_range[0] + brk, including point mass
    g_a_beta_theta(0, i) = pgamma(g_range(0) + brk, shape, scale);
  
    // Middle bins
    for (int j = 1; j < g_tot - 1; j++) {
      Type prob = pgamma(g_range[j] + brk, shape, scale) - 
                  pgamma(g_range[j] - brk, shape, scale);
      g_a_beta_theta(j, i) = prob > Type(0) ? prob : Type(0);
    }
  
    // Last bin: Remaining mass
    // Compute sum from row 0 to g_tot - 2
    Type sum_so_far = Type(0);
    for (int j = 0; j < g_tot - 1; j++) {
      sum_so_far += g_a_beta_theta(j, i);
    }
    g_a_beta_theta(g_tot - 1, i) = Type(1) - sum_so_far > Type(0) ? Type(1) - sum_so_far : Type(0);

  }

  //std::cout << "g_a_beta_theta probabilities" << std::endl;
  //for (year = 0; year < g_a_beta_theta.rows(); year++)
  //{
  //  std::cout << g_a_beta_theta.row(year) << std::endl;
  //}

  // // ===============================
  // // Then the deconvolution
  // //  - Initial we are going to 
  // // plug in about the right answer
  // // ===============================

  matrix<Type> log_E_unconv = Desmat_fix * gammas + Desmat_rand * u;

    // std::cout << "log_E_unconv" << std::endl;
  // for (year = 0; year < log_E_unconv.rows(); year++)
  // {
  //   std::cout << log_E_unconv.row(year) << std::endl;
  // }

  matrix<Type> E_unconv = exp( log_E_unconv.array()); // sans array() it f***s up dims :/

  /*
  vector<int> d_E_unconv(2);
  d_E_unconv(0) = E_unconv. rows();
  d_E_unconv(1) = E_unconv. cols();
  REPORT( d_E_unconv);
  */
  
  matrix<Type> E_count = g_a_beta_theta * E_unconv;

  
  // Probabilities (the real goal here)
  matrix<Type> E_norm( E_unconv.rows(), E_unconv.cols());
  for( int j=0; j < E_unconv.cols(); j++) {
    Type inv_jsum = Type( 1.0) / E_unconv.col( j).sum();
    
    for( int i=0; i < E_unconv.rows(); i++) {
      E_norm( i, j) = E_unconv( i, j) * inv_jsum;
      //std::cout << "E_norm "  <<  E_norm( i, j) << std::endl; 
    };
  };
  
  // std::cout << "E_norm probabilities" << std::endl;
  // for (year = 0; year < E_unconv.rows(); year++)
  // {
  //   std::cout << E_norm.row(year) << std::endl;
  // }

  epi_llk += dpois( count.vec(), E_count.vec(), true).sum();


  // ===============================
  // Compute the ideal HSP probabilities
  // ===============================

  for (byr_1 = 0; byr_1 < yrs_tot; byr_1++)
  {
    for (byr_2 = 0; byr_2 < yrs_tot; byr_2++)
    {
      if (byr_2 > byr_1)
      {
        // Female survival
        Type cumul_surv = exp(-Z1_t * (byr_2 - byr_1)); 
        hsp_prob(byr_1, byr_2) = 4 * fng_rate * cumul_surv / N_ad[byr_2] ;
        // std::cout << "hsp_prob "  <<  hsp_prob(byr_1, byr_2) << std::endl;  
      } else if (byr_2 == byr_1)
      {
        Type cumul_surv = exp(-Z1_t * (byr_2 - byr_1)); 
        hsp_prob(byr_1, byr_2) = c_t * fng_rate * 4 * cumul_surv / N_ad[byr_2];
      }
    }
  }

  // --------------------------------------------------
  // The components to compute P(a|g,y,beta,theta,samp)
  // --------------------------------------------------


  for (int i = 0; i < g_tot; i++)
  {
    for (int j = 0; j < y_tot; j++)
    {
      for (int k = 0; k < a_tot; k++)
      {

        pr_ag_y(i, j, k) = g_a_beta_theta(i, k) * E_norm(k, j);

        //std::cout << "i,j,k "  <<  i << " " << j << " " << k << std::endl;
        //std::cout << "pr_ag_y "  <<  pr_ag_y(i, j, k) << std::endl;  
      }
    }
  }
  
  // Can only slice the last elements of the arrays.

  for (int i = 0; i < a_tot; i++)
  { 
    pr_g_y += pr_ag_y.col(i);
  }
  

  for (int i = 0; i < g_tot; i++)
  {
    for (int j = 0; j < y_tot; j++)
    {
      for (int k = 0; k < a_tot; k++)
      {

        pr_a_gy(i, j, k) = pr_ag_y(i, j, k) / pr_g_y(i, j);

        //std::cout << "i,j,k "  <<  i << " " << j << " " << k << std::endl;
        //std::cout << "pr_a_gy "  <<  pr_a_gy(i, j, k) << std::endl; 
      }
    }
  }
  

  // ----------------------
  // Do the summing over As
  // ----------------------

  //for (int i = 0; i < n_combs; i++)
  //for (int i = 0; i < 2; i++)
  for (int i = 0; i < n_combs; i++)
  {
 
    // if ((i % 2000) == 0) 
    // { 
    //   std::cout << i << std::endl; 
    // }
    
    int G1 = nz_inds(i, 2); 
    int G2 = nz_inds(i, 3);

    int Y1 = nz_inds(i, 4);  
    int Y2 = nz_inds(i, 5); 

    pr_hsap_b1b2gy = 0.0;
    for (int a1 = 0; a1 < a_tot; a1++)
    {
      for (int a2 = 0; a2 < a_tot; a2++)
      {

        int B1 = Y1 - a1 + yrs_tot - y_tot;
        int B2 = Y2 - a2 + yrs_tot - y_tot;

        //std::cout << "a1, a2 "  <<  a1 << " " << a2  << std::endl;
        //std::cout << "B1, B2 "  <<  B1 << " " << B2  << std::endl;
        //std::cout << "G1, G2 "  <<  G1 << " " << G2  << std::endl;
        //std::cout << "Y1, Y2 "  <<  Y1 << " " << Y2  << std::endl;
   

        if (B2 >= B1)
        {
          pr_hsp_b1b2gy_a1a2 = hsp_prob(B1, B2) * pr_a_gy(G1,  Y1, a1) * pr_a_gy(G2, Y2, a2);

          //std::cout << "pr_hsp_b1b2gy_a1a2 "  <<  pr_hsp_b1b2gy_a1a2  << std::endl;

          pr_hsap_b1b2gy += pr_hsp_b1b2gy_a1a2;
        } 

      }
    }
    all_prob_tab(nz_inds(i, 0), nz_inds(i, 1)) = pr_hsap_b1b2gy;
    //std::cout << "prob_matrix "  <<  all_prob_tab(nz_inds(i, 0), nz_inds(i, 1))  << std::endl;

    // Compute the likelihood
    epi_llk += dbinom( kin_yg(nz_inds(i, 0),       nz_inds(i, 1)), 
                       comb_yg(nz_inds(i, 0),      nz_inds(i, 1)), 
                       all_prob_tab(nz_inds(i, 0), nz_inds(i, 1)), true);

  }
  
  // //std::cout << "log-likelihood "  <<  epi_llk  << std::endl;

  // // // ================================
  // // // Log-likelihood HSP - Broad scale
  // // // ================================

  // // std::cout << "Starting HSP llk" << std::endl;
  // // for (byr_1 = 0; byr_1 < byrs_tot; byr_1++)
  // // {
  // //   for (byr_2 = 0; byr_2 < byrs_tot; byr_2++)
  // //   {
  // //     // if (byr_1 <= byr_2)
  // //     if (byr_1 < byr_2)
  // //     {
  // //       // std::cout << "HSP llk " << ciff_llk << std::endl;
  // //       // std::cout << "by1 by2 " << byr_1 << " " << byr_2 << std::endl;
  // //       ciff_llk -= dbinom( kin_sib(byr_1, byr_2),
  // //                           comb_sib(byr_1, byr_2),
  // //                           df_sib_prob(byr_1, byr_2), true);
  // //     }
  // //   }
  // // }

  // // std::cout << "Total llk " << ciff_llk << std::endl;

  if (phi_prior)
  {
    std::cout << "HEEELLP phi PRIOR" << std::endl;
    epi_llk -= pow(log_phi - Mu_phi, 2.0) / (2.0 * pow(SD_phi, 2.0));
  }


  vector<Type> ghat2(31);
  for (int i = 0; i < 31; i++) {

    // if (i <= theta) {

    //   ghat2(i) = alpha + beta * i;  // First segment: y = beta1 * x

    // } else {

    //   ghat2(i) = alpha + beta * theta + bp * (i- theta);  // Second segment: y = beta1 * cp + beta2 * x

    // }

    Type diff = i - theta;
    ghat2(i) = alpha - beta * CppAD::CondExpLe(-theta, Type(0), -theta, Type(0)) 
                    - bp * CppAD::CondExpGe(-theta, Type(0), -theta, Type(0))
                    + beta * CppAD::CondExpLe(diff, Type(0), diff, Type(0)) 
                    + bp * CppAD::CondExpGe(diff, Type(0), diff, Type(0));
    std::cout << "ghat2 i " << ghat2(i) << std::endl;
  }

  //REPORT(hsp_prob);
  //REPORT(all_prob_tab);
  REPORT(epi_llk);
  REPORT(g_a_beta_theta);

  ADREPORT(N0_ad);      // Qbundance year naught
  ADREPORT(Z1);         // Mortality
  ADREPORT(R1);         // PopROI

  ADREPORT(log(N_ad));
  ADREPORT(Z1_t);
  //ADREPORT(E_norm);
  ADREPORT(c_t);
  ADREPORT(ghat2);

  return -epi_llk;
}
