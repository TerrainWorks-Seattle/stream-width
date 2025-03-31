
#include <TMB.hpp>
#include "ssn_prec.hpp"

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}



template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  
  // Data
  DATA_VECTOR( y_n );  // Y at node i 
  DATA_VECTOR( x_n );
  
  DATA_IVECTOR( from_e ); // length n-1
  DATA_IVECTOR( to_e ); // length n-1
  DATA_VECTOR( dist_e ); // length n-1
  DATA_VECTOR( flow_n ); // length n
  
  // Parameters
  PARAMETER( logtheta_up );
  PARAMETER( logtheta_dn );
  PARAMETER( logphi );
  PARAMETER( beta0 );
  PARAMETER( beta1 );
  PARAMETER( gamma_up );
  PARAMETER( gamma_dn );
  
  // Random effects
  PARAMETER_VECTOR( psi_up_n );
  PARAMETER_VECTOR( psi_dn_n );
  
  // Objective function
  vector<Type> jnll(3);
  jnll.setZero();
  
  
  int N = flow_n.size(); 
  
  Type theta_up = exp(logtheta_up);
  Type theta_dn = exp(logtheta_dn);
  Type phi = exp(logphi);
  
  Eigen::SparseMatrix<Type> Q_up(N,N);
  Q_up = ssn_prec_tailup(from_e,
                         to_e,
                         dist_e,
                         flow_n,
                         theta_up );
  
  Eigen::SparseMatrix<Type> Q_dn(N,N);
  Q_dn = ssn_prec_tailup(to_e,
                         from_e,
                         dist_e,
                         flow_n,
                         theta_dn );
  
  // Probability of random effects
  jnll(0) += GMRF(Q_up)( psi_up_n ); 
  jnll(1) += GMRF(Q_dn)( psi_dn_n );
  
  // Probability of data conditional on random effects
  vector<Type> log_mu_n( N );
  vector<Type> mu_n(N);
  for( int n=0; n<N; n++) {
    log_mu_n(n) = beta0 + beta1*x_n(n) + gamma_up*psi_up_n(n) + gamma_dn*psi_dn_n(n);
    mu_n(n) = exp(log_mu_n(n));
    
    if ( !R_IsNA(asDouble(y_n(n))) ){
      jnll(2) -= dgamma( y_n(n), 1/phi, mu_n(n)*phi, true );
    }
  }
  
  // Total jnll
  Type j = jnll(0) + jnll(1) + jnll(2);
  
  // Reporting
  //REPORT( VTriplets );
  REPORT( jnll ); 
  REPORT( Q_up );
  REPORT( Q_dn );
  REPORT( beta0 );
  REPORT( beta1 );
  REPORT( gamma_up );
  REPORT( gamma_dn );
  REPORT( phi );
  REPORT( psi_up_n );
  REPORT( psi_dn_n );
  REPORT( theta_up );
  REPORT( theta_dn );
  REPORT( mu_n );
  
  ADREPORT( mu_n );
  
  return j;
}
