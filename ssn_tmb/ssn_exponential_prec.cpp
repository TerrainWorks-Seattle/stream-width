
#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// // Function for getting precision matrices
// template<class Type>
// Eigen::SparseMatrix<double> ssn_prec_tailup(
//     vector<int> from_e, 
//     vector<int> to_e,
//     vector<double> dist_e,
//     vector<double> flow_n,
//     double theta ){
//   
//   
//   typedef Eigen::Triplet<Type> T;
//   int N = flow_n.size(); 
//   int E = from_e.size();
//   
//   
//   vector<Type> weight(E);
//   weight.fill(1.0);
//   vector<Type> rho(E);
//   vector<Type> var(E);
//   
//   // get confluence nodes
//   vector<int> confl;
//   // nodes with two incoming nodes
//   for (int n=0; n<N; n++) {
//     int i = 0;
//     Type flow_total = 0.0; 
//     for (int e=0; e<E; e++) {
//       if (to_e(e) == n) {
//         i ++;
//         flow_total += flow_n(from_e(e));
//       } 
//     }
//     if (i == 2) {
//       // two nodes flowing into one. 
//       // weight upstream nodes
//       for (int e=0; e<E; e++) {
//         if (to_e(e) == n) {
//           weight(e) = flow_n(from_e(e))/flow_total;
//         }
//       }
//     }
//   }
//   
//   for (int e=0; e<E; e++ ) {
//     // autocorrelation with each upstream node
//     rho(e) = exp(-theta*dist_e(e));
//     // variance contribution from each upstream node
//     var(e) = (1-exp(-2*theta*dist_e(e)));
//   }
//   
//   vector<Type> v_n(N);
//   v_n.fill(1.0);
//   std::vector<T> GammaTriplets(E);
//   for(int e=0; e<E; e++) {
//     // fill Gamma
//     GammaTriplets.push_back(T(to_e(e), from_e(e), weight(e) * rho(e)));
//     // Diagonal variance matrix
//     v_n(to_e(e)) += -1.0 + weight(e) * var(e);
//   }
//   Eigen::SparseMatrix<Type> Gamma(N, N);
//   Gamma.setFromTriplets(GammaTriplets.begin(), GammaTriplets.end());
//   
//   std::vector<T> VTriplets(N); 
//   for (int n=0; n<N; n++) {
//     // this will break if any 0
//     VTriplets.push_back(T(n, n, pow(v_n(n), -1)));
//   }
//   Eigen::SparseMatrix<Type> V(N, N);
//   
//   
//   
//   // Precision
//   Eigen::SparseMatrix<Type> I(N, N);
//   I.setIdentity();
//   Eigen::SparseMatrix<Type> Q = (I-Gamma).transpose() * V * (I-Gamma);
//   return Q;
// }


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
  PARAMETER( logtheta1 );
  //PARAMETER( logtheta2 );
  PARAMETER( logphi );
  PARAMETER( beta0 );
  PARAMETER( beta1 );
  PARAMETER( gamma1 );
  //PARAMETER( gamma2 );
  
  // Random effects
  PARAMETER_VECTOR( psi1_n );
  //PARAMETER_VECTOR( psi2_n );
  
  // Objective function
  vector<Type> jnll(3);
  jnll.setZero();
  
  Type theta1 = exp(logtheta1);
  //Type theta2 = exp(logtheta2);
  Type phi = exp(logphi);
  
  // Make TAILUP precision matrix
  typedef Eigen::Triplet<Type> T;
  int N = flow_n.size(); 
  int E = from_e.size();
  
  
  vector<Type> weight(E);
  weight.fill(1.0);
  vector<Type> rho(E);
  vector<Type> var(E);
  vector<Type> v_n(N);
  v_n.fill(0.0);
  
  // handle confluences and sources. 
  for (int n=0; n<N; n++) {
    int i = 0;
    // sum of upstream contributing areas
    Type flow_total = 0.0; 
    for (int e=0; e<E; e++) {
      if (to_e(e) == n) {
        i ++;
        flow_total += flow_n(from_e(e));
      } 
    }
    // confluence node: two nodes flowing into one. 
    if (i == 2) {
      // set weights for upstream nodes.
      for (int e=0; e<E; e++) {
        if (to_e(e) == n) {
          weight(e) = flow_n(from_e(e))/flow_total;
        }
      }
    }
    // Source node
    // Set endogeneous variance to 1
    if (i == 0) {
      v_n(n) = Type(1.0); 
    }
  }
  
  
  std::vector<T> GammaTriplets;
  GammaTriplets.reserve(E);
  for(int e=0; e<E; e++) {
    // autocorrelation with each upstream node
    rho(e) = exp(-theta1*dist_e(e));
    
    // variance contribution from each upstream node
    var(e) = (1-exp(-2*theta1*dist_e(e)));
    // fill Gamma
    GammaTriplets.push_back(T(to_e(e), from_e(e), weight(e) * rho(e)));
    // Diagonal variance matrix
    v_n(to_e(e)) += weight(e) * var(e);
  }
  Eigen::SparseMatrix<Type> Gamma(N, N);
  Gamma.setFromTriplets(GammaTriplets.begin(), GammaTriplets.end());
  
  std::vector<T> VTriplets; 
  VTriplets.reserve(N);
  for (int n=0; n<N; n++) {
    // this will break if any 0
    VTriplets.push_back(T(n, n, pow(v_n(n), -1)));
  }
  Eigen::SparseMatrix<Type> V(N, N);
  V.setFromTriplets(VTriplets.begin(), VTriplets.end());
  
  
  // Precision
  Eigen::SparseMatrix<Type> I(N, N);
  I.setIdentity();
  Eigen::SparseMatrix<Type> Q1 = (I-Gamma) * V * (I-Gamma).transpose();
  
  
  
  
  // Probability of random effects
  jnll(0) += GMRF(Q1)( psi1_n ); // get density of MVN at omega_s
  
  
  // Probability of data conditional on random effects
  vector<Type> log_mu_n( N );
  vector<Type> mu_n(N);
  for( int n=0; n<N; n++) {
    log_mu_n(n) = beta0 + beta1*x_n(n) + gamma1*psi1_n(n);
    mu_n(n) = exp(log_mu_n(n));
    
    if ( !R_IsNA(asDouble(y_n(n))) ){
      jnll(1) -= dgamma( y_n(n), 1/phi, mu_n(n)*phi, true );
    }
  }
  
  // Total jnll
  Type j = jnll(0) + jnll(1) + jnll(2);
  
  // Reporting
  REPORT( V ); 
  //REPORT( VTriplets );
  REPORT( jnll ); 
  REPORT( Gamma );
  REPORT( I ); 
  REPORT( Q1 );
  REPORT( weight );
  REPORT( rho );
  REPORT( var );
  REPORT( v_n );
  REPORT( beta0 );
  REPORT( beta1 );
  REPORT( gamma1 );
  REPORT( phi );
  REPORT( psi1_n );
  REPORT( theta1 );
  REPORT( mu_n );
  
  ADREPORT( mu_n );
  
  return j;
}
