#include <RcppArmadillo.h>
#include <iostream>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List achr(const Rcpp::List& model, const Rcpp::List& state, const arma::mat& warmupPnts, const int nPnts, const int stepsPerPnt) {
    
    std::cout << "Preparing for ACHR sampling..." << std::endl;

    const double maxMinTol = 1e-09, uTol = 1e-09, dTol = 1e-14;
    arma::vec lb = Rcpp::as<arma::vec>(model["lb"]);
    arma::vec ub = Rcpp::as<arma::vec>(model["ub"]);
    arma::sp_mat S = Rcpp::as<arma::sp_mat>(model["S"]);
    arma::mat* Sd = new arma::mat(S);
    arma::mat N = arma::null(*Sd);
    delete Sd;
    unsigned int nRxns = S.n_cols;
    unsigned int nWarmupPnts = warmupPnts.n_cols;
    arma::vec centerPnt = Rcpp::as<arma::vec>(state["center.pnt"]);
    arma::vec prevPnt = Rcpp::as<arma::vec>(state["prev.pnt"]);
    unsigned int totalStepCount = Rcpp::as<Rcpp::IntegerVector>(state["n.tot.steps"])[0];

    std::cout << "Begin ACHR sampling..." << std::endl;

    arma::mat points(nRxns, nPnts);
    arma::vec curPnt;
    arma::vec randVector;
    unsigned int pointCount = 0;
    while (pointCount < nPnts) {
        
        // create the random step size vector
        randVector = arma::randu(stepsPerPnt);
    
        unsigned int stepCount = 0;
        while (stepCount < stepsPerPnt) {
    
            // pick a random warmup point
            int randPntID = std::floor(nWarmupPnts*arma::randu());
            arma::vec randPnt = warmupPnts.col(randPntID);
            // get a direction from the center point to the warmup point
            arma::vec u = randPnt - centerPnt;
            u = u / norm(u);
            // figure out the distances to upper and lower bounds
            arma::vec distUb = ub - prevPnt;
            arma::vec distLb = prevPnt - lb;
            // figure out if we are too close to a boundary
            arma::uvec validDir = arma::find((distUb > dTol) % (distLb > dTol)); // element-wise multiplication (%) as a workaround for &
            // figure out positive and negative directions
            arma::vec validU = u(validDir);
            arma::uvec posDirn = arma::find(validU > uTol);
            arma::uvec negDirn = arma::find(validU < -uTol);
            // figure out all the possible maximum and minimum step sizes
            arma::vec inverseValidU = 1 / validU;
            arma::vec maxStepTemp = distUb(validDir) % inverseValidU;
            arma::vec minStepTemp = -distLb(validDir) % inverseValidU;
            arma::vec maxStepVec = arma::join_vert(maxStepTemp(posDirn), minStepTemp(negDirn));
            arma::vec minStepVec = arma::join_vert(minStepTemp(posDirn), maxStepTemp(negDirn));
            // figure out the true max & min step sizes
            double maxStep = maxStepVec.min();
            double minStep = minStepVec.max();
            // find new direction if we're getting too close to a constraint
            if ((std::abs(minStep) < maxMinTol && std::abs(maxStep) < maxMinTol) || (minStep > maxStep)) continue;
            // pick a rand out of list_of_rands and use it to get a random step distance
            double stepDist = randVector(stepCount)*(maxStep-minStep)+minStep;
            // advance to the next point
            curPnt = prevPnt + stepDist*u;
            // reproject the current point
            if (totalStepCount % 10 == 0) {
                if (arma::abs(S*curPnt).max() > 1e-9) {
                  curPnt = N*(N.t()*curPnt);
                }
            }
            arma::uvec overInd = ub < curPnt;
            arma::uvec underInd = lb > curPnt;
            if (any(overInd) || any(underInd)) {
                arma::uvec overIndi = arma::find(overInd);
                arma::uvec underIndi = arma::find(underInd);
                curPnt(overIndi) = ub(overIndi);
                curPnt(underIndi) = lb(underIndi);
            }
            // recalculate the center point
            centerPnt = ((nWarmupPnts+totalStepCount)*centerPnt + curPnt)/(nWarmupPnts+totalStepCount+1);
            // next
            prevPnt = curPnt;
            totalStepCount = totalStepCount + 1;
            stepCount = stepCount + 1;
        }
    
        // add the current point to points
        points.col(pointCount) = curPnt;
        pointCount = pointCount + 1;
        // print progress
        if (pointCount % (nPnts/10) == 0) {
            std::cout << "Sampling progress: " << 100*pointCount/nPnts << "%" << std::endl;
        }
    }

    std::cout << "Finished ACHR sampling." << std::endl;

    Rcpp::List endStat = Rcpp::List::create(Rcpp::Named("center.pnt") = centerPnt,
                                            Rcpp::Named("prev.pnt") = prevPnt,
                                            Rcpp::Named("n.tot.steps") = totalStepCount);
    return Rcpp::List::create(Rcpp::Named("pnts") = points,
                              Rcpp::Named("stat") = endStat);
}