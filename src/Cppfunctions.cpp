


# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// This file contains all the C++ functions used


//Generate Random variable form t distribution
//Generate a random of size m by n with degree of freedom k
// [[Rcpp::export]]
arma::mat randt(int M, int N, int K)
{
  using namespace arma;
  int i, j;
  float v1, v2;
  mat T;       T.zeros(M, N);
  vec X;       X.zeros(K+1);

  for(i=0; i<M; i++){
    for(j=0; j<N; j++){
      X=randn(K+1);
      v1=mean(X);  v2=var(X);
      T(i,j)=sqrt(static_cast<double>(K+1))*v1/v2;
    }
  }

  return T;

}
//Generate Fourier basis
// [[Rcpp::export]]
arma::vec Fourier_basis(float z, int n)
{
  using namespace arma;
  int i;
  float v1,v2;
  float w=0.04;
  v1=sqrt(static_cast<double>(2)); v2=1.0*z*3.1415926*w;//v2=1.0*datum::pi*z;
  vec cita; cita.zeros(n);
  for(i=1;i<=n-1;i++){
    if(i%2==0)cita(i)=v1*cos(i*v2);
    else      cita(i)=v1*sin((i+1)*v2);
  }
  cita(0)=1;
  return cita;

}



//Calculate Huber Loss
// [[Rcpp::export]]
arma::mat Huber_loss (arma::mat X, arma::mat phi, arma::mat B, float CT, int T)
{
  using namespace arma;
  int t;       float v1;
  mat Loss;    Loss.zeros(1,1);
  mat Huber;   Huber.zeros(1,1);
  mat M1;      M1.zeros(1,1);

  for(t=0; t<T; t++){
    M1=X(t)-phi.row(t)*B;
    v1=fabs(M1(0));
    if(v1> fabs(CT)) {Huber(0)=CT*(2*M1(0)-CT);}
    else {Huber(0)= M1(0)*M1(0);}
    //printf("\n x=%f  x^2=%f   Huber=%f", M1(0), M1(0)*M1(0), Huber(0));
    Loss+=Huber/T;
  }

  return Loss;
}



//Calculate the Gradient of Huber Loss
// [[Rcpp::export]]
arma::mat Huber_gradient (arma::mat X, arma::mat phi, arma::mat B, float CT, int T)
{
  using namespace arma;
  int t, J=B.n_rows;       float v1;
  mat Grad;        Grad.zeros(J,1);
  mat Huber_dot;   Huber_dot.zeros(1,1);
  mat M1;          M1.zeros(1,1);
  mat phi_t;       phi_t=phi.t();

  for(t=0; t<T; t++){
    M1=X(t)-phi.row(t)*B;
    v1=M1(0);
    if (v1> fabs(CT))     {Huber_dot(0)=2*CT;}
    else if (v1 < -1*fabs(CT)) {Huber_dot(0)= -2*CT;}
    else {Huber_dot(0)=2*M1(0);}
    //printf("\n x=%f   Huber_dot=%f", M1(0), Huber_dot(0));
    Grad-= 1*Huber_dot(0)*phi_t.col(t)/T;
  }

  return Grad;
}

//Minimize Huber loss with gradient descent
// [[Rcpp::export]]
arma::mat Huber_descent (arma::mat X, arma::mat phi, arma::mat B, float CT)
{
  using namespace arma;

  int k, T=phi.n_rows, J=phi.n_cols;
  float v1=0, v2=0;
  mat b_1;  b_1.zeros(J,1);
  mat b_2;  b_2.zeros(J,1);
  mat test; test.zeros(J,1);
  b_1=B;

  for(k=1; k<500; k++){
    v1= as_scalar(Huber_loss (X, phi, b_1, CT, T));
    test=Huber_gradient (X, phi, b_1, CT, T);
    b_2=b_1;      b_1-=0.5* test/sqrt(static_cast<double>(k));
    v2=as_scalar(Huber_loss (X, phi, b_1, CT, T));
    //printf("\n %dth v1=%f    v2=%f   \n", k, v1, v2);
    if(fabs(v1-v2)<1.0e-10 || v1<v2+1.0e-8)k=500;
  }

  return b_2;

}




//Robust estimation with Huber loss
// [[Rcpp::export]]
float Robust_CV (arma::mat vx, arma::mat phi)
{
  using namespace arma;
  int i,k,T=phi.n_rows, J=phi.n_cols,  T_vali=0;
  float MSE_vali, MSE_small, ct_o=5, ct;
  T_vali=T/5;

  mat vx_1;      mat vx_2;
  mat vx_train;  mat vx_vali;    mat vx_hat;
  mat phi_1;     mat phi_2;
  mat phi_vali;  mat phi_train;
  vec b_2;       vec b_sol;
  vec b_1;       b_1.zeros(J,1);

  b_1=solve(phi, vx);

  for(i=1,MSE_small=1.0e8, ct_o=1.0e8; i<=21; i++){

    ct=(abs(vx).max())*i/20;
    if(i==21)ct=1.0e7;


    for(k=0, MSE_vali=0; k<5; k++){
      //printf("\n---------------  %dth ------------------\n",k);
      vx_1.resize(0,0);    vx_2.resize(0,0);
      vx_vali=vx.rows(span(k*T_vali, (k+1)*T_vali-1));
      if(k > 0)vx_1=vx.rows(span(0,k*T_vali-1));
      if(k < 4)vx_2=vx.rows(span((k+1)*T_vali,T-1));
      vx_train=join_cols(vx_1,vx_2);
      phi_1.resize(0,0);    phi_2.resize(0,0);
      phi_vali=phi.rows(span(k*T_vali, (k+1)*T_vali-1));
      if(k > 0)phi_1=phi.rows(span(0,k*T_vali-1));
      if(k < 4)phi_2=phi.rows(span((k+1)*T_vali,T-1));
      phi_train=join_cols(phi_1,phi_2);

      b_2 = Huber_descent (vx_train, phi_train, b_1, ct);
      //b_sol=solve(phi_train, vx_train);
      vx_hat=phi_vali*b_2;

      MSE_vali+=as_scalar((vx_hat-vx_vali).t()*(vx_hat-vx_vali)/T_vali);

      //cout << b_2 << endl;
      //cout << b_sol << endl;
    }
    if(MSE_vali<MSE_small){MSE_small=MSE_vali; ct_o=ct;}
    //printf("\n  CT=%f,   MSE_vali=%f \n", ct, MSE_vali);

  }

  return ct_o;

}






//Robust estimation with Huber loss
// [[Rcpp::export]]
arma::mat Robust_estimate (arma::mat X, arma::mat phi, arma::mat B, float CT)
{
  using namespace arma;

  int i, T=phi.n_rows, J=phi.n_cols, N=X.n_rows;
  mat b_1;    b_1.zeros(J,1);
  mat b_2;    b_2.zeros(J,1);
  mat b_sol;  b_sol.zeros(J,1);
  mat vx;     vx.zeros(T,1);
  mat B_hat;  B_hat.zeros(N,J);
  b_1=B;

  for(i=0; i<N; i++){
    //printf("\n---------- i=%d --------------------\n", i);
    vx=X.row(i).t();
    b_2 = Huber_descent (vx, phi, b_1, CT);
    b_sol=solve(phi, vx);
    B_hat.row(i)=b_2.t();
    //cout << b_2 << endl;
    //cout << b_sol << endl;
    //printf("\n------------------------------------\n");
  }

  return B_hat;

}


// [[Rcpp::export]]
arma::mat mu_robust_F( arma::mat X, arma::mat phi)
{
  using namespace arma;
  int i, P, K;
  P=X.n_rows;
  K=phi.n_cols;

  float Tau;
  mat F_H_0; F_H_0.ones(K);
  mat Xi;
  mat mu_hat; mu_hat.zeros(K,P);

  for(i=0;i<P;i++){
    Rcpp::checkUserInterrupt();
    Xi=X.row(i);
    F_H_0=solve(phi,trans(Xi));
    Tau= Robust_CV (trans(Xi),(phi));
    mu_hat.col(i)=Huber_descent(Xi, phi, F_H_0, Tau);
  }
  return mu_hat;

}


//Output: Estimated cov matrix Sigma_hat
// [[Rcpp::export]]
arma::mat Cov_Huber( arma::mat X, arma::mat phi)
{
  using namespace arma;
  int i, j, P=X.n_rows;

  //Tuning parameter
  float Tau;

  //Define the matrices
  mat Xi, Xj;
  mat Sigma_hat; Sigma_hat.zeros(P,P);
  mat F_H_0; F_H_0.ones(1);

  //Entry-wise Huber method
  for(i=0;i<P;i++){
    Rcpp::checkUserInterrupt();
    for(j=0;j<=i;j++){
      Xi=X.row(i); Xj=X.row(j);
      F_H_0=solve(phi,trans(Xi%Xj));
      Tau= Robust_CV (trans(Xi%Xj), phi);
      Sigma_hat(i,j)=arma::conv_to<double>::from(Huber_descent (Xi%Xj,phi, F_H_0, Tau));
      Sigma_hat(j,i)=Sigma_hat(i,j);
    }
  }


  return Sigma_hat;

}




//Find factors
// [[Rcpp::export]]
arma::mat Find_factors (arma::mat Sigma, arma::mat X, int N, int P, int K)
{
  using namespace arma;


//mat XX;             XX.zeros(N, N);
  vec eigval_cov;     eigval_cov.zeros(N);
  mat eigvec_cov;     eigvec_cov.zeros(N,N);
  mat F_hat;          F_hat.zeros(N,K);
  mat Lambda_hat;     Lambda_hat.zeros(P,K);


  //PCA on XPX, eatimate loadings by the first K eigenvectors
  //XX= X *  X.t();

  eig_sym(eigval_cov, eigvec_cov, Sigma);
  eigval_cov=sort(eigval_cov,"descend");
  eigvec_cov=fliplr(eigvec_cov);
  Lambda_hat=eigvec_cov.cols(0,K-1)* sqrt(static_cast<double>(P));

  //Estimate Factors
  F_hat=X.t()*Lambda_hat/P;

  //Generate Y_star and U_hat
  mat I_n;  I_n.eye(N,N);
  mat P_F; P_F=F_hat.t()*F_hat;
  P_F=I_n-F_hat*P_F.i()*F_hat.t();
  return P_F;
}



// [[Rcpp::export]]
arma::mat Find_PF(arma::mat F_hat, int N)
{
  using namespace arma;

  //Generate X_star
  mat I_n;  I_n.eye(N,N);
  mat P_F; P_F=F_hat.t()*F_hat;
  P_F=I_n-F_hat*P_F.i()*F_hat.t();
  return P_F;
}


// [[Rcpp::export]]
arma::mat Find_lambda_class (arma::mat Sigma, arma::mat X, int N, int P, int K)
{
  using namespace arma;


  //mat XX;             XX.zeros(N, N);
  vec eigval_cov;     eigval_cov.zeros(N);
  mat eigvec_cov;     eigvec_cov.zeros(N,N);
  mat Lambda_hat;     Lambda_hat.zeros(P,K);


  //PCA on XPX, eatimate loadings by the first K eigenvectors
  //XX= X *  X.t();

  eig_sym(eigval_cov, eigvec_cov, Sigma);
  eigval_cov=sort(eigval_cov,"descend");
  eigvec_cov=fliplr(eigvec_cov);
  Lambda_hat=eigvec_cov.cols(0,K-1)* sqrt(static_cast<double>(P));
  return Lambda_hat;
}

// [[Rcpp::export]]
arma::mat Find_factors_class (arma::mat Lambda_hat, arma::mat X, int N, int P, int K)
{
  using namespace arma;
  //Estimate Factors
  mat F_hat;          F_hat.zeros(N,K);
  F_hat=X.t()*Lambda_hat/P;
  return F_hat;
}

// [[Rcpp::export]]
arma::mat Find_X_star_class (arma::mat F_hat,arma::mat Lambda_hat,  arma::mat X)
{
  using namespace arma;
  mat U_star;
  U_star= (X-Lambda_hat*F_hat.t());
  mat UT;      UT=U_star.t();
  return UT;

}



//Find factors
// [[Rcpp::export]]
arma::mat Find_X_star(arma::mat P_F, arma::mat X)
{   using namespace arma;
 mat U_star;  U_star=X*P_F;
  //U_star= (X-Lambda_hat*F_hat.t())*P_F;
  mat UT;      UT=U_star.t();
  return UT;
}



// [[Rcpp::export]]
arma::mat Find_Y_star (arma::mat P_F,arma::mat Y)
{   using namespace arma;
  mat Y_star;
  //Y_star=Y;
  Y_star=P_F*Y;
  return Y_star;
}

//Input:  Covariance matrix M
//Output: all
// [[Rcpp::export]]
arma::mat Eigen_Decomp( arma::mat M)
{
  using namespace arma;
  int  P=M.n_rows;

  //Define matrices for eigenvalues and eigen-vectors
  vec eigval_cov;  eigval_cov.zeros(P);
  mat eigvec_cov;  eigvec_cov.zeros(P,P);
  mat eigall_cov;  eigall_cov.zeros(P,P+1);

  eig_sym(eigval_cov, eigvec_cov, M);
  eigval_cov=sort(eigval_cov,"descend");
  eigvec_cov=fliplr(eigvec_cov);
  eigall_cov = join_rows(eigvec_cov, eigval_cov);

  return eigall_cov;

}


