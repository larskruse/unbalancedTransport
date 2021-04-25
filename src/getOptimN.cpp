#include "Rcpp.h"
#include "getOptimN.h"



void getNtSW(Rcpp::NumericVector ukSW, Rcpp::NumericVector vkSW, Rcpp::List indexGroupsSWi, Rcpp::List indexGroupsSWj, Rcpp::NumericVector iList, Rcpp::NumericVector jList, int &nt, int &k0, int &k1, double &shiftTransportCost){

    // indexes of the supply and demand vectors with positive mass
    Rcpp::NumericVector uniqueI = Rcpp::unique(iList).sort();
    Rcpp::NumericVector uniqueJ = Rcpp::unique(jList).sort();

    // vector to store the minimum transport costs for each "shift group"
    Rcpp::NumericVector minKSW;

    // calculate best shift in each "shift group"
    for(int i = 0; i < indexGroupsSWi.length(); i++){
        // add new element to minKSW
        minKSW.push_back(
            // minimum value in the ukSW values indicated by the "shift group"
            Rcpp::min(
                Rcpp::as<Rcpp::NumericVector>(ukSW[
                                                  Rcpp::in(uniqueJ, Rcpp::as<Rcpp::NumericVector>(indexGroupsSWj[i]))]
                )
            )
        // maximum value in the vkSW values indicated by the "shift group"
        - Rcpp::max(
                Rcpp::as<Rcpp::NumericVector>(vkSW[
                                                  Rcpp::in(uniqueI, Rcpp::as<Rcpp::NumericVector>(indexGroupsSWi[i]))]
                )
        )
        );

    }


    // index of best shift value
    nt = Rcpp::which_min(minKSW);
    // best shift value
    shiftTransportCost = minKSW[nt];


    // k1 and k0 are the values in the indexGroups... lists that are used to
    // compute the optimal shift
    k1 = Rcpp::which_min(
        Rcpp::as<Rcpp::NumericVector>(ukSW[
                                          Rcpp::in(uniqueJ, Rcpp::as<Rcpp::NumericVector>(indexGroupsSWj[nt]))]
        ));
    k1 = Rcpp::as<Rcpp::NumericVector>(indexGroupsSWj[nt])[k1]; //-1];


    k0 = Rcpp::which_max(
        Rcpp::as<Rcpp::NumericVector>(vkSW[
                                          Rcpp::in(uniqueI, Rcpp::as<Rcpp::NumericVector>(indexGroupsSWi[nt]))]
        ));
    k0 = Rcpp::as<Rcpp::NumericVector>(indexGroupsSWi[nt])[k0]; //-1];
}



void getNtNO(Rcpp::NumericVector ukNO, Rcpp::NumericVector vkNO, Rcpp::List indexGroupsNOi, Rcpp::List indexGroupsNOj, Rcpp::NumericVector iList, Rcpp::NumericVector jList, int &nt, int &k0, int &k1, double &shiftTransportCost){

    // indexes of the supply and demand vectors with positive mass
    Rcpp::NumericVector uniqueI = Rcpp::unique(iList).sort();
    Rcpp::NumericVector uniqueJ = Rcpp::unique(jList).sort();

    // vector to store the minimum transport costs for each "shift group"
    Rcpp::NumericVector minKNO;


    // calculate best shift in each "shift group"
    for(int i = 0; i < indexGroupsNOi.length(); i++){
        // add new element to minKNO


        minKNO.push_back(
            // minimum value in the ukNO values indicated by the "shift group"
            Rcpp::min(
                Rcpp::as<Rcpp::NumericVector>(ukNO[
                                                  Rcpp::in(uniqueI, Rcpp::as<Rcpp::NumericVector>(indexGroupsNOi[i]))]
                ))
        // maximum value in the vkNO values indicated by the "shift group"
        -Rcpp::max(
                Rcpp::as<Rcpp::NumericVector>(vkNO[
                                                  Rcpp::in(uniqueJ, Rcpp::as<Rcpp::NumericVector>(indexGroupsNOj[i]))]
                )));

    }


    // index of best shift value
    nt = Rcpp::which_min(minKNO);
    shiftTransportCost = minKNO[nt];


    // k1 and k0 are the values in the indexGroups... lists that are used to
    // compute the optimal shift
    k0 = Rcpp::which_min(
        Rcpp::as<Rcpp::NumericVector>(ukNO[
                                          Rcpp::in(uniqueI, Rcpp::as<Rcpp::NumericVector>(indexGroupsNOi[nt]))]
        ));
    k0 = Rcpp::as<Rcpp::NumericVector>(indexGroupsNOi[nt])[k0]; //-1];

    k1 = Rcpp::which_max(
        Rcpp::as<Rcpp::NumericVector>(vkNO[
                                          Rcpp::in(uniqueJ, Rcpp::as<Rcpp::NumericVector>(indexGroupsNOj[nt]))]
        ));
    k1 = Rcpp::as<Rcpp::NumericVector>(indexGroupsNOj[nt])[k1]; //-1];


}



void getNtSWNO(Rcpp::NumericVector ukSW, Rcpp::NumericVector vkSW, Rcpp::List indexGroupsSWi, Rcpp::List indexGroupsSWj, Rcpp::NumericVector ukNO, Rcpp::NumericVector vkNO, Rcpp::List indexGroupsNOi, Rcpp::List indexGroupsNOj, Rcpp::NumericVector iList, Rcpp::NumericVector jList, int &nt, int &k0, int &k1, double &shiftTransportCost){

    // variables to store the results of getNtNO and getNtSW
    int k0SW;
    int k1SW;

    int k0NO;
    int k1NO;

    int ntSW;
    int ntNO;

    double transportCostSW;
    double transportCostNO;



    // calculate the optimal north-east and south-west shift
    getNtSW(ukSW, vkSW, indexGroupsSWi, indexGroupsSWj, iList, jList, ntSW, k0SW, k1SW, transportCostSW);

    getNtNO(ukNO, vkNO, indexGroupsNOi, indexGroupsNOj, iList, jList, ntNO, k0NO, k1NO, transportCostNO);

    // the shift with the lower cost is the overall optimal shift
    // change reference values accordingly
    if(transportCostSW > transportCostNO){
        shiftTransportCost = transportCostNO;
        k0 = k0NO;
        k1 = k1NO;
        nt = ntNO;
    }else{
        shiftTransportCost = transportCostSW;
        k0 = k0SW;
        k1 = k1SW;
        nt = ntSW;

    }

    return;
}
