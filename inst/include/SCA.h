

//template<class Type>
//Type objective_function<Type>::operator() ()
//{
  using namespace SCA;

  DATA_VECTOR(C_hist);    // Total catch
  DATA_VECTOR(I_hist);    // Index

  DATA_MATRIX(CAA_hist);  // Catch-at-length proportions
  DATA_VECTOR(CAA_n);     // Annual samples in CAL

  DATA_MATRIX(CAL_hist);  // Catch-at-length proportions
  DATA_VECTOR(CAL_n);     // Annual samples in CAL
  DATA_INTEGER(n_y);      // Number of years in model
  DATA_INTEGER(max_age);  // Maximum age (plus-group)
  DATA_VECTOR(M);         // Natural mortality at age

  DATA_VECTOR(mat);       // Maturity-at-age at the beginning of the year

  DATA_STRING(I_type);    // String whether index surveys B, VB, or SSB
  DATA_STRING(SR_type);   // String indicating whether Beverton-Holt or Ricker stock-recruit is used

  DATA_VECTOR(est_early_rec_dev);
  DATA_VECTOR(est_rec_dev); // Indicator of whether rec_dev is estimated in model or fixed at zero

  DATA_INTEGER(Nbins);
  DATA_VECTOR(LenBins);   // Length bins
  DATA_VECTOR(LenMids);   // Mid of length bins

  DATA_SCALAR(Linf);      // Linf
  DATA_SCALAR(min_LAA);   // Smallest LAA for selectivity
  DATA_INTEGER(yind_F);
  
  DATA_INTEGER(use_LeesEffect);
  
  // Only used if use_LeesEffect = TRUE
  DATA_INTEGER(ngtg);     // Number of growth type groups
  DATA_VECTOR(distGTG);
  DATA_VECTOR(rdist);
  
  DATA_MATRIX(WAA);       // Weight-at-age-and-GTG at the beginning of the year
  DATA_MATRIX(LAA);       // Length-at-age-and-GTG at the beginning of the year
  
  DATA_MATRIX(xout);      // Length bins and Length-at-age sorted by row
  DATA_IMATRIX(interp_check);
  DATA_IMATRIX(interp_check2);
  DATA_IMATRIX(integ_check);
  DATA_IVECTOR(integ_fac);
  DATA_IVECTOR(integ_ind);
  
  // Only used if use_LeesEffect = FALSE
  DATA_VECTOR(mean_LAA);
  DATA_VECTOR(mean_WAA);
  DATA_SCALAR(CV_LAA);

  PARAMETER(log_R0);
  PARAMETER(transformed_h);
  PARAMETER(F_equilibrium);
  PARAMETER_VECTOR(vul_par);

  PARAMETER_VECTOR(logF);

  PARAMETER(log_sigma);
  PARAMETER(log_omega);
  PARAMETER(log_tau);
  PARAMETER_VECTOR(log_early_rec_dev);
  PARAMETER_VECTOR(log_rec_dev);

  Type R0 = exp(log_R0);
  Type h;
  if(SR_type == "BH") {
    h = 0.8 * invlogit(transformed_h);
  } else {
    h = exp(transformed_h);
  }
  h += 0.2;

  Type sigma = exp(log_sigma);
  Type omega = exp(log_omega);
  Type tau = exp(log_tau);

  Type penalty = 0;
  Type prior = 0.;

  // Calculate selectivity-at-length
  Type LFS = invlogit(vul_par(0)) * (0.9 * Linf - min_LAA) + min_LAA;
  Type L5 = LFS - exp(vul_par(1)) * Linf;
  Type Vmaxlen = invlogit(vul_par(2));

  Type sl = (LFS - L5) / pow(-log2(0.05), 0.5);
  Type sr = (Linf - LFS) / pow(-log2(Vmaxlen), 0.5);

  vector<Type> Select_at_length = s_dnormal(LenMids, LFS, sl, sr, Vmaxlen);
  
  // Define F
  vector<Type> F(n_y);
  F(yind_F) = exp(logF(yind_F));
  for(int y=0;y<n_y;y++) if(y != yind_F) F(y) = F(yind_F) * exp(logF(y));
  
  ////// Equilibrium reference points and per-recruit quantities
  vector<Type> NPR_virgin(max_age);
  vector<Type> Weight_virgin(max_age);
  Weight_virgin.setZero();
  
  vector<Type> NPR_equilibrium(max_age);
  vector<Type> Weight_equilibrium(max_age);
  Weight_equilibrium.setZero();

  // Probability of length given age
  vector<matrix<Type> > probLA(n_y+1);
  vector<matrix<Type> > NPR(n_y+1);
  vector<matrix<Type> > probGTGA(n_y+1);
  matrix<Type> SAA(max_age, ngtg);
  
  matrix<Type> Select_at_age(n_y+1, max_age);
  matrix<Type> Weight_at_age(n_y+1, max_age);
  Select_at_age.setZero();
  Weight_at_age.setZero();
  
  // Setup probLA
  // Equilibrium quantities (leading into first year of model) 
  if(use_LeesEffect) { // GTG submodel
    vector<vector<int> > integ_index = split(integ_ind, integ_fac);
    SAA = s_dnormal(LAA, LFS, sl, sr, Vmaxlen); // Selectivity at age and GTG
    
    NPR_virgin = LeesApp_fn(Type(0), rdist, M, SAA, WAA, Weight_virgin, max_age, ngtg);
    NPR_equilibrium = LeesApp_fn(F_equilibrium, rdist, M, SAA, WAA, Weight_equilibrium, max_age, ngtg);
    
    for(int y=0;y<=n_y;y++) { // Also updates Weight_at_age, Select_age, NPR, probGTGA
      probLA(y) = LeesApp_fn(F, F_equilibrium, rdist, M, SAA, LenBins, LAA, WAA, xout, Select_at_length, 
             Select_at_age, Weight_at_age, NPR, probGTGA, Nbins, max_age, ngtg, y, interp_check, interp_check2,
             integ_check, integ_index);
    }
    
  } else { // no GTGs
    probLA(0) = LenAge_matrix(LenMids, mean_LAA, CV_LAA, max_age, LenMids.size(), LenMids(1) - LenMids(0)); //get Select_at_age now
    
    vector<Type> Select_at_age2(max_age);
    Select_at_age2.setZero();
    for(int a=0;a<max_age;a++) {
      for(int len=0;len<Nbins;len++) Select_at_age2(a) += probLA(0)(a,len) * Select_at_length(len);
    }
    
    NPR_virgin = calc_NPR(Type(0), Select_at_age2, M, max_age);
    Weight_virgin = mean_WAA;
    
    NPR_equilibrium = calc_NPR(F_equilibrium, Select_at_age2, M, max_age);
    Weight_equilibrium = mean_WAA;
    
    for(int y=0;y<=n_y;y++) {
      if(y>0) probLA(y) = probLA(0);
      Select_at_age.row(y) = Select_at_age2;
      Weight_at_age.row(y) = mean_WAA;
    }
  }

  Type EPR0 = 0;
  Type B0 = 0;
  for(int a=0;a<max_age;a++) {
    EPR0 += NPR_virgin(a) * Weight_virgin(a) * mat(a);
    B0 += NPR_virgin(a) * Weight_virgin(a);
  }
  
  B0 *= R0;
  Type N0 = R0 * NPR_virgin.sum();
  Type E0 = R0 * EPR0;

  Type Arec;
  Type Brec;

  if(SR_type == "BH") {
    Arec = 4 *h;
    Arec /= 1-h;
    Arec /= EPR0;
    Brec = 5*h - 1;
    Brec /= (1-h) * E0;
  } else {
    Arec = pow(5*h, 1.25);
    Arec /= EPR0;
    Brec = 1.25;
    Brec *= log(5*h);
    Brec /= E0;
  }
  Type CR = Arec * EPR0;
  
  Type EPR_eq = 0;
  for(int a=0;a<max_age;a++) {
    EPR_eq += NPR_equilibrium(a) * Weight_equilibrium(a) * mat(a);
  }
  
  Type R_eq;
  if(SR_type == "BH") {
    R_eq = Arec * EPR_eq - 1;
  } else {
    R_eq = log(Arec * EPR_eq);
  }
  R_eq /= Brec * EPR_eq;

  ////// During time series year = 1, 2, ..., n_y
  matrix<Type> N(n_y+1, max_age);   // Numbers at year and age
  matrix<Type> CAApred(n_y, max_age);   // Catch (in numbers) at year and age at the mid-point of the season
  matrix<Type> CALpred(n_y, Nbins);
  vector<Type> CN(n_y);             // Catch in numbers
  vector<Type> Cpred(n_y);
  vector<Type> Ipred(n_y);          // Predicted index at year
  vector<Type> R(n_y+1);            // Recruitment at year
  vector<Type> R_early(max_age-1);
  vector<Type> VB(n_y+1);           // Vulnerable biomass at year
  vector<Type> B(n_y+1);            // Total biomass at year
  vector<Type> E(n_y+1);            // Spawning biomass at year

  CN.setZero();
  CALpred.setZero();
  Cpred.setZero();
  VB.setZero();
  B.setZero();
  E.setZero();
  
  // Year one
  R(0) = R_eq;
  if(!R_IsNA(asDouble(est_rec_dev(0)))) R(0) *= exp(log_rec_dev(0) - 0.5 * pow(tau, 2));

  for(int a=0;a<max_age;a++) {
    if(a == 0) {
      N(0,a) = R(0) * NPR_equilibrium(a);
    } else {
      R_early(a-1) = R_eq;
      if(!R_IsNA(asDouble(est_early_rec_dev(a-1)))) R_early(a-1) *= exp(log_early_rec_dev(a-1) - 0.5 * pow(tau, 2));
      N(0,a) = R_early(a-1) * NPR_equilibrium(a);
    }
    
    B(0) += N(0,a) * Weight_at_age(0,a);
    VB(0) += N(0,a) * Weight_at_age(0,a) * Select_at_age(0,a);
    E(0) += N(0,a) * Weight_at_age(0,a) * mat(a);
    
  }

  // Loop over all other years
  for(int y=0;y<n_y;y++) {
    if(SR_type == "BH") {
      R(y+1) = BH_SR(E(y), h, R0, E0);
    } else {
      R(y+1) = Ricker_SR(E(y), h, R0, E0);
    }

    if(y<n_y-1) {
      if(!R_IsNA(asDouble(est_rec_dev(y+1)))) R(y+1) *= exp(log_rec_dev(y+1) - 0.5 * pow(tau, 2));
    }
    N(y+1,0) = R(y+1);
    
    for(int a=0;a<max_age;a++) {
      if(a<max_age-1) N(y+1,a+1) = N(y,a) * exp(-Select_at_age(y,a) * F(y) - M(a));
      Type meanN = N(y,a) * (1 - exp(-Select_at_age(y,a) * F(y) - M(a))) / (Select_at_age(y,a) * F(y) + M(a));
      CAApred(y,a) = Select_at_age(y,a) * F(y) * meanN;
      CN(y) += CAApred(y,a);
      for(int len=0;len<Nbins;len++) CALpred(y,len) += probLA(y)(a,len) * CAApred(y,a);
   
      Cpred(y) += CAApred(y,a) * Weight_at_age(y,a);
      B(y+1) += N(y+1,a) * Weight_at_age(y+1,a);
      VB(y+1) += N(y+1,a) * Weight_at_age(y+1,a) * Select_at_age(y+1,a);
      E(y+1) += N(y+1,a) * Weight_at_age(y+1,a) * mat(a);
      
    }
  }

  // Calculate nuisance parameters and likelihood
  Type q;
  if(I_type == "B") {
    q = calc_q(I_hist, B);
    for(int y=0;y<n_y;y++) Ipred(y) = q * B(y);
  } else if(I_type == "VB") {
    q = calc_q(I_hist, VB);
    for(int y=0;y<n_y;y++) Ipred(y) = q * VB(y);
  } else {
    q = calc_q(I_hist, E);
    for(int y=0;y<n_y;y++) Ipred(y) = q * E(y);
  }

  vector<Type> nll_comp(5);
  nll_comp.setZero();
  for(int y=0;y<n_y;y++) {
    if(!R_IsNA(asDouble(I_hist(y)))) nll_comp(0) -= dnorm(log(I_hist(y)), log(Ipred(y)), sigma, true);
    if(C_hist(y) > 0) {
      if(!R_IsNA(asDouble(CAA_n(y))) && CAA_n(y) > 0) {
        vector<Type> loglike_CAAobs = CAA_n(y) * CAA_hist.row(y);
        vector<Type> loglike_CAApred = CAApred.row(y)/CN(y);
        nll_comp(1) -= dmultinom_robust(loglike_CAAobs, loglike_CAApred, true);
      }

      if(!R_IsNA(asDouble(CAL_n(y))) && CAL_n(y) > 0) {
        vector<Type> loglike_CALobs = CAL_n(y) * CAL_hist.row(y);
        vector<Type> loglike_CALpred = CALpred.row(y)/CN(y);
        nll_comp(2) -= dmultinom_robust(loglike_CALobs, loglike_CALpred, true);
      }
      nll_comp(3) -= dnorm(log(C_hist(y)), log(Cpred(y)), omega, true);
    }
    if(!R_IsNA(asDouble(est_rec_dev(y)))) nll_comp(4) -= dnorm(log_rec_dev(y), Type(0), tau, true);
  }
  for(int a=0;a<max_age-1;a++) {
    if(!R_IsNA(asDouble(est_early_rec_dev(a)))) nll_comp(4) -= dnorm(log_early_rec_dev(a), Type(0), tau, true);
  }

  Type nll = nll_comp.sum() + penalty + prior;

  ADREPORT(R0);
  ADREPORT(h);
  ADREPORT(omega);
  ADREPORT(sigma);
  ADREPORT(tau);
  ADREPORT(q);

  REPORT(omega);
  REPORT(sigma);
  REPORT(tau);

  REPORT(NPR_virgin);
  REPORT(Arec);
  REPORT(Brec);
  REPORT(EPR0);
  REPORT(CR);
  REPORT(h);
  REPORT(R0);
  REPORT(B0);
  REPORT(N0);
  REPORT(E0);
  
  REPORT(NPR_equilibrium);
  REPORT(EPR_eq);
  REPORT(R_eq);

  REPORT(vul_par);
  REPORT(L5);
  REPORT(LFS);
  REPORT(Vmaxlen);
  REPORT(sl);
  REPORT(sr);

  REPORT(F);

  REPORT(N);
  REPORT(CN);
  REPORT(Cpred);
  REPORT(CAApred);
  REPORT(CALpred);
  REPORT(Ipred);
  REPORT(R);
  REPORT(R_early);
  REPORT(VB);
  REPORT(B);
  REPORT(E);

  REPORT(probLA);
  if(use_LeesEffect) {
    REPORT(probGTGA);
    REPORT(NPR);
  }
  REPORT(SAA);
  REPORT(Select_at_length);
  REPORT(Select_at_age);
  REPORT(Weight_at_age);
  
  REPORT(Weight_virgin);
  REPORT(Weight_equilibrium);

  REPORT(log_early_rec_dev);
  REPORT(log_rec_dev);
  REPORT(nll_comp);
  REPORT(nll);
  REPORT(penalty);
  REPORT(prior);

  return nll;
//}
