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
                                                                         

  DATA_SCALAR(fng_rate);                            // Intercept for the linear model
  DATA_INTEGER(sig_prior);                           // sd prior?
  DATA_SCALAR(Mu_log_sig_ga);                            // Normal prior mean for ROI if required
  DATA_SCALAR(SD_log_sig_ga);                            // Normal prior variance for ROI if required

  PARAMETER(N0_ad);                              // Abundance year naught
  PARAMETER(Z1);                                 // Female/or and male survival
  PARAMETER(R1);                                 // Female/or and male ROI
  PARAMETER(alpha);
  PARAMETER(Beta);                                // Slope term to calibration epi-age to age
  

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
  PARAMETER(log_sig_ga);         // Sd of g-score given age
  PARAMETER(log_c);          // Variance of g-score given age

  Type Z1_t = exp(Z1);
  Type c_t  = exp(log_c);
  Type sig_ga = exp(log_sig_ga);

  // Setup the random effect (u)
  // ---------------------------

  Type epi_llk = 0.0;

  matrix<Type> u_indept = chol_S * u;
  epi_llk += dnorm( u_indept.vec(), Type(0), exp(logsdu), true).sum();

  // Known age likelihood
  // --------------------

  vector<Type> ghat = alpha + Beta * kwn_a;
  epi_llk += dnorm(knw_g, ghat, sig_ga, true).sum();

  
  // --------------------
  // Initial declarations
  // --------------------

  int yrs_tot, g_tot, a_tot, y_tot, byr_1, byr_2, year, gy_tot, n_combs; //, gy_tot, byrs_tot, syrs_tot, year, byr_1, byr_2, syrs_diff, yad, bju, rel_yad, sib_mt_it;
  Type exptd_g, pr_hsp_b1b2gy_a1a2, pr_hsap_b1b2gy, brk;

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


  // ===============
  // Pop dynamics
  // ===============

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
  for (int i = 0; i < a_tot; i++)
  {
     exptd_g = alpha + a_range[i] * Beta;
    
    // First bin: Full CDF up to g_range[0] + brk, including point mass
    g_a_beta_theta(0, i) = pnorm(g_range(0) + brk, exptd_g, sig_ga);
  
    // Middle bins
    for (int j = 1; j < g_tot - 1; j++) {
      Type prob = pnorm(g_range[j] + brk, exptd_g, sig_ga) - 
                  pnorm(g_range[j] - brk, exptd_g, sig_ga);
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

  std::cout << "g_a_beta_theta probabilities" << std::endl;
  for (year = 0; year < g_a_beta_theta.rows(); year++)
  {
   std::cout << g_a_beta_theta.row(year) << std::endl;
  }

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
  

  if (sig_prior)
  {
    std::cout << "HEEELLP Sig PRIOR" << std::endl;
    epi_llk -= pow(log_sig_ga - Mu_log_sig_ga, 2.0) / (2.0 * pow(SD_log_sig_ga, 2.0));
    // ciff_llk += pow(Z2 - Mu_Z1, 2.0) / (2.0 * pow(SD_Z1, 2.0));
  }


  vector<Type> ghat2(a_tot);
  for (int i = 0; i < a_tot; i++) {

    ghat2(i) = alpha + Beta * i;
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
