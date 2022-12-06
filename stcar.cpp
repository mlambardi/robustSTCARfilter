#include <Rcpp.h>
using namespace Rcpp;

// This is the export of a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar).

// Function to subset a data vector in a R-like fashion
// It creates a vector with the only components of v indexed by i
NumericVector getdata(NumericVector v, IntegerVector i) {
  NumericVector w(i.length());
  for (int j = 0; j < w.length(); j++) {
    w[j] = v[i[j]];
  }
  return w;
}

// Function to overwrite a vector and avoid exceeding with instances
// It replaces all components of w with c(v1, v2, v3)
NumericVector overwr(NumericVector w, NumericVector v1, double v2, double v3, double v4) {
  for (int i = 0; i < v1.length(); i++) {
    w[i] = v1[i];
  }
  w[v1.length()] = v2;
  w[v1.length()+1] = v3;
  w[v1.length()+2] = v4;
  return w;
}

// Function to update the Hessian matrix in filtering, diagonal version
// It adds v to the diagonal of m^{-1}
NumericMatrix updatediag(NumericMatrix m, NumericVector v, bool seas) {
  for (int i = 0; i < m.ncol() - !seas; i++) {
    m(i,i) = 1.0/(1.0/m(i,i) + v[i]);
  }
  return m;
}

// It computes the matrix-times-vector product m%*%v
NumericVector prodmatvec(NumericMatrix m, NumericVector v) {
  NumericVector pr(m.nrow());
  for (int i = 0; i < m.nrow(); i++) {
    pr[i] = sum(m(i,_)*v);
  }
  return pr;
}

// Function to update the Hessian matrix in filtering, non-diagonal version
// It computes the inverse of m^{-1} + l*v%*%t(v)
NumericMatrix updatematr(NumericMatrix m, NumericVector v, double l) {
  NumericVector pr = prodmatvec(m, v);
  double den = 1.0/l + sum(v*pr);
  for (int i = 0; i < m.nrow(); i++) {
    m(i,_) = m(i,_) - pr*pr[i]/den;
  }
  return m;
}

// Sign function
double sign(double x) {
  if (x > 0) {
    return +1;
  } else {
    if (x < 0) {
      return -1;
    } else {
      return 0;
    }
  }
}

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
// [[Rcpp::export]]
List updater(NumericVector par, double s1, double s2, double dens, double lambda, int nloc, int tlag, double teff, IntegerVector indexy, IntegerMatrix indexx, NumericMatrix dat, NumericMatrix datspat, NumericMatrix infox, NumericMatrix errmat, NumericMatrix parmat, NumericMatrix predmat, NumericVector uvec, NumericVector tracevec, int starttime, bool is_median, bool is_diagonal, bool verbose, int tseas, bool seas) {
  List output;
  int p = par.length();
  NumericVector score(p);
  NumericVector xvec(p);
  double news1, news2, newdens;
  int i, time;
  int ry;
  IntegerVector rx;
  NumericVector objxst(p-3);
  NumericVector objxsT(p-3);
  NumericVector objxSt(p-3);
  NumericVector objxST(p-3);
  NumericVector objxstseas(p-3);
  NumericVector objxsTseas(p-3);
  NumericVector objxStseas(p-3);
  NumericVector objxSTseas(p-3);
  double objyst, objysT, objySt, objyST;
  double objystseas, objysTseas, objyStseas, objySTseas;
  double rho = par[p-3], phi = par[p-2], phiseas = par[p-1];
  NumericVector beta = par[Range(0, p-4)];
  double err, cnst;
  starttime = starttime + tseas;
  Progress progr(dat.nrow()-starttime, verbose);
  for (time = starttime; time < dat.nrow(); time++) {
    score.fill(0.0);
    news1 = (1.0 - lambda*nloc)*s1;
    news2 = (1.0 - lambda*nloc)*s2;
    newdens = (1.0 - lambda*nloc)*dens;
    for (i = 0; i < nloc; i++) {
      ry = indexy[i];
      rx = indexx(i,_);
      objyst =     dat(time-0*tlag,ry);
      objysT =     dat(time-1*tlag,ry);
      objySt = datspat(time-1*tlag,ry);
      objyST = datspat(time-2*tlag,ry);
      objxst = getdata(dat(    time-1*tlag,_), rx);
      objxsT = getdata(dat(    time-2*tlag,_), rx);
      objxSt = getdata(datspat(time-2*tlag,_), rx);
      objxST = getdata(datspat(time-3*tlag,_), rx);
      objystseas =     dat(time-0*tlag-tseas,ry);
      objysTseas =     dat(time-1*tlag-tseas,ry);
      objyStseas = datspat(time-1*tlag-tseas,ry);
      objySTseas = datspat(time-2*tlag-tseas,ry);
      objxstseas = getdata(dat(    time-1*tlag-tseas,_), rx);
      objxsTseas = getdata(dat(    time-2*tlag-tseas,_), rx);
      objxStseas = getdata(datspat(time-2*tlag-tseas,_), rx);
      objxSTseas = getdata(datspat(time-3*tlag-tseas,_), rx);
      
      xvec = overwr(
        xvec,
        (objxst - rho*objxSt) - phi*(objxsT - rho*objxST) - phiseas*((objxstseas - rho*objxStseas) - phi*(objxsTseas - rho*objxSTseas)),
        (objySt - phi*objyST) - sum(beta*(objxSt - phi*objxST)) - phiseas*((objyStseas - phi*objySTseas) - sum(beta*(objxStseas - phi*objxSTseas))),
        (objysT - rho*objyST) - sum(beta*(objxsT - rho*objxST)) - phiseas*((objysTseas - rho*objySTseas) - sum(beta*(objxsTseas - rho*objxSTseas))),
        ((objystseas - phi*objysTseas) - rho*(objyStseas - phi*objySTseas) - sum(beta*((objxstseas - phi*objxsTseas) - rho*(objxStseas - phi*objxSTseas))))*seas
      );
      
      err = objyst - rho*objySt - phi*(objysT - rho*objyST) - sum(beta*((objxst - rho*objxSt) - phi*(objxsT - rho*objxST)))
        -phiseas*(objystseas - rho*objyStseas - phi*(objysTseas - rho*objySTseas) - sum(beta*((objxstseas - rho*objxStseas) - phi*(objxsTseas - rho*objxSTseas))));
      
      news1 = news1 + lambda*fabs(err);
      news2 = news2 + lambda*err*err;
      newdens = newdens + lambda*R::dnorm4(err/s1, 0.0, 1.33*pow(lambda, 0.2), FALSE);
      
      if (is_median) {
        err = sign(err);
        cnst = 2.0*dens/s1;
      } else {
        // err = err / s2;
        // cnst = 1.0 / s2;
        cnst = 1.0;
      }
      
      score = score + xvec*err*lambda;
      
      if (is_diagonal) {
        infox = updatediag(infox, xvec*xvec * cnst*lambda/(1.0 - lambda*nloc), seas);
      } else {
        infox = updatematr(infox, xvec, cnst*lambda/(1.0 - lambda*nloc));
      }
    }
    infox = infox / (1.0 - lambda*nloc);
    s1 = news1;
    s2 = news2;
    dens = newdens;
    for (i = 0; i < p; i++) {
      parmat(time, i) = par[i];
    }
    parmat(time, p) = s1;
    parmat(time, p+1) = sqrtf(s2);
    parmat(time, p+2) = dens;
    uvec = prodmatvec(infox, uvec);
    tracevec[time] = 1.0/sqrtf(sum(uvec*uvec));
    uvec = uvec * tracevec[time];
    par = par + prodmatvec(infox, score);
    
    beta = par[Range(0,p-4)];
    rho = par[p-3];
    phi = par[p-2];
    phiseas = par[p-1];
    
    if (time + 1*tlag < dat.nrow()) {
      for (i = 0; i < nloc; i++) {
        ry = indexy[i];
        rx = indexx(i,_);
        objyst =     dat(time+1*tlag,ry);
        objysT =     dat(time-0*tlag,ry);
        objySt = datspat(time-0*tlag,ry);
        objyST = datspat(time-1*tlag,ry);
        objxst = getdata(dat(    time-0*tlag,_), rx);
        objxsT = getdata(dat(    time-1*tlag,_), rx);
        objxSt = getdata(datspat(time-1*tlag,_), rx);
        objxST = getdata(datspat(time-2*tlag,_), rx);
        objystseas =     dat(time+1*tlag-tseas,ry);
        objysTseas =     dat(time-0*tlag-tseas,ry);
        objyStseas = datspat(time-0*tlag-tseas,ry);
        objySTseas = datspat(time-1*tlag-tseas,ry);
        objxstseas = getdata(dat(    time-0*tlag-tseas,_), rx);
        objxsTseas = getdata(dat(    time-1*tlag-tseas,_), rx);
        objxStseas = getdata(datspat(time-1*tlag-tseas,_), rx);
        objxSTseas = getdata(datspat(time-2*tlag-tseas,_), rx);
        errmat(time + 1*tlag,i) = (objyst - phi*objysT) - rho*(objySt - phi*objyST) - sum(beta*((objxst - phi*objxsT) - rho*(objxSt - phi*objxST)))
          -phiseas*((objystseas - phi*objysTseas) - rho*(objyStseas - phi*objySTseas) - sum(beta*((objxstseas - phi*objxsTseas) - rho*(objxStseas - phi*objxSTseas))));
        predmat(time + 1*tlag,i) = objyst - errmat(time + 1*tlag,i);
      }
    }
    progr.increment();
  }
  
  output["pars"] = parmat;
  output["mineigen"] = tracevec;
  output["pred"] = predmat;  
  output["resid"]=errmat;
  return output;
}

/*** R
stcar <- function(f, dat, tlag=1, tseas=0, teff, is.median=T, is.diagonal=F, wmat, verbose = F) {
  dat <- as.matrix(dat)
  seas <- tseas > 0
  
  varresp <- gsub("`", "", formula.tools::lhs.vars(f))
  varcova <- gsub("`", "", formula.tools::rhs.vars(f))
  varalll <- c(varresp, varcova)
  varspat <- varalll[grepl("%d", varalll)]
  
  suff <- unique(drop(as.vector(sapply(varspat, function(i) {
    fmt <- paste0("^", gsub("%d", "(\\\\d+)", i), "$")
    sub(fmt, "\\1", colnames(dat)[grepl(fmt, colnames(dat))])
  }))))
  suff <- suff[order(as.numeric(suff))]
  nloc <- length(suff)
  lambda <- 1/(teff * nloc)
  wmat <- wmat[suff,suff]
  
  datspat <- dat
  for (i in varspat) {
    i <- sprintf(sub("%d", "%s", i), suff)
    datspat[,i] <- as.matrix(dat[,i]) %*% wmat
  }
  
  indexy <- as.vector(sapply(sprintf(sub("%d", "%s", varresp), suff), function(i) which(colnames(dat) == i)))
  indexx <- sapply(varcova, function(v) sapply(if (grepl("%d", v)) sprintf(sub("%d", "%s", v), suff) else rep(v, length(suff)), function(j) which(colnames(dat)==j)))
  rownames(indexx) <- suff
  
  inity <- head(dat, ceiling(teff))[,indexy]
  
  beta <- sapply(varcova, function(i) 0)
  if ("intercept%d" %in% varcova) {
    drift <- beta[which(varcova=="intercept%d")] <- mean(inity, na.rm=T)
  } else {
    drift <- 0
  }
  rho <- 0
  phi <- 0
  phiseas <- 0
  s1 <- mean(abs(inity - drift), na.rm=T)
  s2 <- mean((inity - drift)^2, na.rm=T)
  dens <- mean(dnorm((inity-drift)/s1, sd = 1.33*lambda^.2), na.rm=T)
  
  infox <- apply(indexx, 2, function(i) {
    mean(head(dat, ceiling(teff))[,i]^2, na.rm=T)
  })
  infox <- c(infox, rep(s2, 3))
  infox <- infox * ifelse(is.median, 2*dens/s1, 1)
  infox <- 1/infox
  if (!seas) infox[length(infox)] <- 0
  infox <- diag(infox)
  uvec <- rep(1/sqrt(ncol(infox)), length=ncol(infox))
  
  parmat <- matrix(NA, ncol = ncol(infox) + 3, nrow = nrow(dat))
  colnames(parmat) <- c(colnames(indexx), "rho", "phi", "phiseas", "gamma", "RMSE", "f(0)")
  predmat <- errmat <- matrix(NA, ncol = nloc, nrow = nrow(dat))
  tracevec <- rep(NA, length=nrow(dat))
  
  par <- c(beta, rho=rho, phi=phi, phiseas=phiseas)
  
  updater(
    par = par, s1 = s1, s2 = s2, dens = dens, lambda = lambda, nloc = nloc, tlag = tlag, teff = teff,
    indexy = indexy-1, indexx = indexx-1, dat = dat, datspat = datspat, errmat = errmat, parmat = parmat, predmat = predmat, uvec = uvec, tracevec = tracevec,
    infox = infox, starttime = max(ceiling(teff) + 1, 3*tlag + 1)-1, is_median = is.median, is_diagonal = is.diagonal, verbose = verbose,
    tseas = tseas, seas = seas
  ) %>%
    purrr::list_modify(
      formula = as.character(f),
      tlag = tlag, tseas = if (seas) tseas else NA, teff = teff,
      is.median = is.median, is.diagonal = is.diagonal,
      locs = suff,
      varalll = varalll,
      varspat = varspat,
      dat = dat,
      datspat = datspat
    )
}
*/
