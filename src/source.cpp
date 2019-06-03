#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// #######################################################
// ###### The logLikelihood functions
// #######################################################

// [[Rcpp::export]]
double logLikMk_C(arma::uword n_nodes, arma::uword n_tips, arma::uword n_states, arma::vec edge_len, arma::mat edge_mat, arma::vec parents, arma::mat X, arma::mat Q, int root_node, int root_type) {
  // This is the log-lik function for a simple Mk model fitted to the tree. This will return the likelihood of the transition matrix for the rate regimes.

  // n_nodes = number of nodes in phy
  // n_tips = number of tips in phy
  // n_states = number of states in the regime.
  // edge_len = vector with the edge length
  // edge_mat = the edge matrix for the tree.
  // parents = vector with the unique( edge.matrix[,1] )
  // X = data matrix. number of columns equal to the number of states. rows in the same order as the tip.labels. a 1 marks state is present and a 0 mark state absent.
  // Q = the transition matrix.
  // root_node = the node number for the root node
  // root_type = the type to compute the root probabilities. 0 = equal and 1 = madfitz

  arma::mat append_mat = mat(n_nodes, n_states, fill::zeros);
  arma::mat liks = join_vert(X, append_mat);
  arma::vec comp = vec(n_nodes + n_tips);
  arma::uword anc;
  arma::uvec ii; // This is a vector of indexes.
  arma::mat v = mat(n_states, 2); // Two descendant nodes.

  // Loop to traverse the tree.
  for(uword i=0; i < n_nodes; i++) {

    // Need to check the usage of 'anc'. Is it an index or a vector test?
    anc = parents[i] - 1; // This is an index. C++ starts from 0.
    ii = find( parents[i] == edge_mat.col(0) ); // More than one entry.

    uword des;
    arma::vec v_root = vec(n_states, fill::ones);
    for(uword j=0; j < 2; j++) {
      des = as_scalar( edge_mat(ii[j], 1) ) - 1; // This is an index
      v.col(j) = expmat(Q * edge_len[ ii[j] ]) * trans( liks.row( des ) );
      v_root = v_root % v.col(j);
    }
	  
    if( parents[i] == root_node ){ // The computations at the root
      if( root_type == 0 ){
	// This is the equal root probability model:
	arma::vec equal_pi = vec(n_states, fill::ones);
	equal_pi = equal_pi / n_states;
	comp[ anc ] = sum( v_root % equal_pi );
      } else{
	// if( root_type == 1) // This needs to be TRUE
	// This is the Maddison and Fitzjohn method.
	// arma::vec comp_unscaled = sum( v_root.col(0) );
	arma::vec liks_root = v_root / sum( v_root );
	arma::vec root_p = liks_root / sum( liks_root );
	comp[ anc ] = sum( v_root % root_p );
      }
    } else{
      comp[ anc ] = sum( v_root );
    }

    liks.row(anc) = trans( v_root / comp[ anc ] ); // Need row vector.
  }

  // Get the log-lik for the model:
  return sum( log( comp.subvec(0 + n_tips, n_nodes + n_tips - 1) ) );
}

double priorQ(arma::vec vec_Q, arma::vec par_prior_Q, std::string den_Q){
  // Compute the prior probability for the Q matrix rates.
  // Can be uniform or exponential prior.
  // The 'par_prior_Q' can be a vector with 1 or 2 elements, depeding on the density.
  double pp = 0.0;
  if( den_Q == "uniform" ){
    for( arma::uword i=0; i < vec_Q.n_rows; i++ ){ // vec_Q is a column.
      // If "unif" then 'par_prior_Q' is a vector with 2 elements.
      pp = pp + R::dunif(vec_Q[i], par_prior_Q[0], par_prior_Q[1], true);
    }
  } else{ // Then it is an exponential prior.
    for( arma::uword i=0; i < vec_Q.n_rows; i++ ){
      // In Rcpp the exponential is '1/rate' whereas in R it is simply 'rate' .
      // Here only one parameter is needed.
      pp = pp + R::dexp(vec_Q[i], 1/par_prior_Q[0], true);
    }
  }
  
  return pp;
}

arma::vec extractQ(arma::mat Q, arma::uword size, std::string model_Q){
  // Function to extract a column vector from the Q matrix.
  // Length of the vector will depend on the type of the model for the Q matrix.
  // Need to use the same pattern to extract and rebuild the matrix.
  arma::vec vec_Q;
  
  if( model_Q == "ER" ){
    vec_Q = vec(1, fill::zeros);
    vec_Q[0] = Q(0,1); // All off-diagonals are the same.
  } else if( model_Q == "SYM" ){
    arma::uword count = 0;
    int size_vec = ( ( size * size ) - size ) / 2;
    vec_Q = vec(size_vec, fill::zeros);
    for( arma::uword i=0; i < size; i++ ){
      for( arma::uword j=0; j < size; j++ ){
	if( i >= j ) continue;
	vec_Q[count] = Q(i,j);
	count++;
      }
    }
  } else{ // model_Q == "ARD"
    arma::uword count = 0;
    int size_vec = ( size * size ) - size;
    vec_Q = vec(size_vec, fill::zeros);
    for( arma::uword i=0; i < size; i++ ){
      for( arma::uword j=0; j < size; j++ ){
	if( i == j ) continue;
	vec_Q[count] = Q(i,j);
	count++;
      }
    }
  }

  return vec_Q;
}

arma::mat buildQ(arma::vec vec_Q, arma::uword size, std::string model_Q){
  // Function to re-build the Q matrix.
  // Need to follow the same pattern used to extract the vector.
  arma::mat Q = mat(size, size, fill::zeros);
  
  if( model_Q == "ER" ){
    Q.fill(vec_Q[0]);
    // Now fill the diagonal.
    for( arma::uword i=0; i < size; i++ ){
      Q(i,i) = -1.0 * ( sum( Q.row(i) ) - vec_Q[0] );
    }
  } else if( model_Q == "SYM" ){
    Q.fill(0); // Fill the matrix with 0.
    arma::uword count = 0;
    // Go over the matrix and fill the upper and lower-tri.
    for( arma::uword i=0; i < size; i++ ){
      for( arma::uword j=0; j < size; j++ ){
	if( i >= j ) continue;
	Q(i,j) = vec_Q[count];
	Q(j,i) = vec_Q[count]; // The trick to fill the lower-tri.
	count++;
      }
    }
    // Now fill the diagonal.
    for( arma::uword i=0; i < size; i++ ){
      Q(i,i) = -1.0 * sum( Q.row(i) );
    }
  } else{ // model_Q == "ARD"
    Q.fill(0); // Fill the matrix with 0.
    arma::uword count = 0;
    // Go over the matrix and fill the upper and lower-tri.
    for( arma::uword i=0; i < size; i++ ){
      for( arma::uword j=0; j < size; j++ ){
	if( i == j ) continue;
	Q(i,j) = vec_Q[count];
	count++;
      }
    }
    // Now fill the diagonal.
    for( arma::uword i=0; i < size; i++ ){
      Q(i,i) = -1.0 * sum( Q.row(i) );
    }
  }

  return Q;
}

void writeQToFile(std::ostream& Q_mcmc_stream, arma::vec vec_Q, int k, std::string model_Q){
  // Note the 'std::ostream&' argument here is the use of a reference.
  if( model_Q == "ER" ){
    Q_mcmc_stream << vec_Q;
  } else{
    arma::uword print_size = vec_Q.n_rows;
    for( arma::uword i=0; i < (print_size-1); i++ ){
      Q_mcmc_stream << vec_Q[i];
      Q_mcmc_stream << "; ";
    }
    Q_mcmc_stream << vec_Q[print_size-1];
    Q_mcmc_stream << "\n";
  }
}

