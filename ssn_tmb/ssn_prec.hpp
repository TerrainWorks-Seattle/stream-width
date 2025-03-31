
// Function for getting precision matrices
template<class Type>
Eigen::SparseMatrix<Type> ssn_prec_tailup(
    vector<int> from_e,
    vector<int> to_e,
    vector<Type> dist_e,
    vector<Type> flow_n,
    Type theta ){
  
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
    rho(e) = exp(-theta*dist_e(e));
    
    // variance contribution from each upstream node
    var(e) = (1-exp(-2*theta*dist_e(e)));
    // fill Gamma
    GammaTriplets.push_back(T(from_e(e), to_e(e), weight(e) * rho(e)));
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
  Eigen::SparseMatrix<Type> Q = (I-Gamma) * V * (I-Gamma).transpose();
  
  return Q;
}