
#include "Rcpp.h"
#include "getOptimN.h"
#include "utilMonge.h"
#include "shifts.h"

//' C++ implementation of the Monge algorithm
//'
//' This function calculates the optimal transport cost and transport
//' plan for unbalanced optimal transport problems with a Monge cost matrix
//'
//'
//' @param costMatrix A numeric matrix that fulfills the Monge property
//' @param supply A numeric vector. Its length must be equal to the numbers of rows in the cost matrix.
//' @param demand A numeric vector. Its length must be equal to the number of column in the cost matrix.
//' @param constructionCost The cost of creating mass at any point.
//' @param destructionCost The cost of destructing mass at any point.
//' @return It returns a list of the transport cost, three numeric vectors that indicate the optimal transport path, and two vectors that indicate
//'     where mass is created and where it is destructed.
//' @noRd
//[[Rcpp::export]]
Rcpp::List Monge_Rcpp(Rcpp::NumericMatrix &costMatrix, Rcpp::NumericVector &supply,
                      Rcpp::NumericVector &demand, double &constructionCost,
                      double &destructionCost){

    // storing the original supply and demand vectors
    // the variables p and q will be altered to give a balanced transport problem
    Rcpp::NumericVector supplyTrue = Rcpp::clone(supply);
    Rcpp::NumericVector demandTrue = Rcpp::clone(demand);

    double supplySum = Rcpp::sum(supply);
    double demandSum = Rcpp::sum(demand);

    // expand the cost matrix, supply, and demand vectors to an balanced transport problem
    UnbalancedToBalanced(costMatrix, supply, demand);


    //difference between mass supply and demand of the balanced and unbalanced problem
    double difCreateMass = Rcpp::sum(demand) - demandSum;
    double difDestructMass =  Rcpp::sum(supply) - supplySum;

    // amount of mass that is created and destructed
    // will be updated in each iteration depending on the the optimal shift
    double createdMass = 0.0;
    double destructedMass = 0.0;

    // the optimal transport plan is similar the the path given by the north west corner rule
    // it is represented by three vectors:
    //  iList: the row index of the path elements
    //  jList: the column index of the path elements
    //  weightList: the weight transported at the position given by the iList and jList entries
    Rcpp::NumericVector iList;
    Rcpp::NumericVector jList;
    Rcpp::NumericVector weightList;

    // calculate initial transport cost of the balanced problem
    double transportCost;
    transportCost = Nw_Corner_Rule(costMatrix, supply, demand, iList, jList, weightList);

    // Variables to store values of the previous iteration
    // The main loop is repeated until the transport cost does not decline anymore.
    // Then, the values from the previous interation give to optimal transport solution
    double transportCostPrev = transportCost + 1;

    Rcpp::NumericVector supplyPrev = Rcpp::clone(supply);
    Rcpp::NumericVector demandPrev = Rcpp::clone(demand);
    Rcpp::NumericVector weightListPrev = Rcpp::clone(weightList);
    Rcpp::NumericVector iListPrev = Rcpp::clone(iList);
    Rcpp::NumericVector jListPrev = Rcpp::clone(jList);

    // Indicates which element of the current path has the highest cost
    // Eliminating one element from the path is an alternative to a shift
    int pathMinIndex;
    double pathMinTransportCost;


    // variables to store the information about the optimal shift
    // the variables are named according to the paper
    double shiftMinTransportCost;
    int nt;
    int k0;
    int k1;
    double shiftTransportCost;


    // storing information on which shifts in the south-west and north-east
    // directions are possible
    Rcpp::List indexGroupsSWi;
    Rcpp::List indexGroupsSWj;
    Rcpp::List indexGroupsNOi;
    Rcpp::List indexGroupsNOj;

    // u(k) and v(k) vectors according to the paper
    Rcpp::NumericVector ukSW;
    Rcpp::NumericVector vkSW;
    Rcpp::NumericVector ukNO;
    Rcpp::NumericVector vkNO;

    // amount of mass that is removed from the supply and demand vectors by the shfit
    double pathChangeValue;

    // import and export vectors
    Rcpp::NumericVector importVec;
    Rcpp::NumericVector exportVec;



    // main iteration of the function
    // each iteration computes a new optimal shift and changes the supply/demand vectors accordingly
    // the new transport cost is than compared to the cost computed in the previous iteration
    // if the new cost is higher, the optimal transport plan was the one given by the last iteration

    // since the original problem is an unbalanced transport, at least the mass added by
    // making it a balanced transport problem has to be removed, even if that increases the
    // transport cost

    while(transportCostPrev >= transportCost || Rcpp::sum(supplyPrev) > supplySum || Rcpp::sum(demandPrev) > demandSum){


        // save the variables obtained by the last iteration
        iListPrev = Rcpp::clone(iList);
        jListPrev = Rcpp::clone(jList);
        weightListPrev = Rcpp::clone(weightList);
        supplyPrev = Rcpp::clone(supply);
        demandPrev = Rcpp::clone(demand);
        transportCostPrev = transportCost;


        // if the transport path is not NULL, removing a single element from the path is also a
        // valid "shift"
        // the cost for that is computed independently of the shifts between distinct points in the path
        if(iList.length() > 0){


            pathMinTransportCost = -costMatrix(iList[0], jList[0]);
            pathMinIndex = 0;

            for(int i = 1; i < iList.length(); i++){

                if(-costMatrix(iList[i], jList[i]) < pathMinTransportCost){
                    pathMinTransportCost = -costMatrix(iList[i], jList[i]);
                    pathMinIndex = i;
                }

            }
        }else{

            break;
        }




        // if there is more than 1 element in the transport path, a shift is possible
        if(iList.length() > 1){

            // calculate between which points a shift can be applied
            GetIndexGroups(iList, jList, indexGroupsSWi, indexGroupsSWj, indexGroupsNOi, indexGroupsNOj);



            // calculate the u(k), v(k) vectors according to the paper
            CalcSWShift(costMatrix, iList, jList, ukSW, vkSW);
            CalcNOShift(costMatrix, iList, jList, ukNO, vkNO);


            // if south-west AND north-east shifts are possible find the overall best shift
            if(indexGroupsSWi.length() > 0 && indexGroupsNOi.length() > 0){

                getNtSWNO(ukSW, vkSW, indexGroupsSWi, indexGroupsSWj, ukNO, vkNO, indexGroupsNOi, indexGroupsNOj, iList, jList, nt, k0, k1, shiftTransportCost);

                // if north-east shifts are not possible, calculate the best south-west shift
            }else if(indexGroupsSWi.length() > 0){


                getNtSW(ukSW, vkSW, indexGroupsSWi, indexGroupsSWj, iList, jList, nt, k0, k1, shiftTransportCost);

                // if south-west shifts are not possible, calculate the best north-east shift
            }else if(indexGroupsNOi.length() > 0){


                getNtNO(ukNO, vkNO, indexGroupsNOi, indexGroupsNOj, iList, jList, nt, k0, k1, shiftTransportCost);

                // else set nt to -1 to indicate, that no shifts are possible
            }else{
                nt = -1;
            }

            // no shift is possible if there is only one element in the transport path
        }else{
            nt = -1;
        }



        // the mass shifted is equal to the minimum mass in any element of the transport path
        // this makes sure, that the shifts add at most one element at every turning of the path
        pathChangeValue = Rcpp::min(weightList);

        // if no shift is possible, or the the cost induced by the shift is larger than the cost
        // of only removing one element from the given path, remove only the one element

        if((nt >= 0 && shiftTransportCost > pathMinTransportCost) || nt < 0 ){
            k0 = iList[pathMinIndex];
            k1 = jList[pathMinIndex];
            pathChangeValue = weightList[pathMinIndex];
        }


        // update the supply and demand vectors
        supply[k0] = supply[k0] - pathChangeValue;
        demand[k1] = demand[k1] - pathChangeValue;

        // set to 0 if the fall under the threshold. This eliminates floating
        // point errors
        if(supply[k0] < 1e-10){
            supply[k0] = 0;
        }

        if(demand[k1] < 1e-10){
            demand[k1] = 0;
        }



        createdMass += pathChangeValue;
        destructedMass += pathChangeValue;

        // the new transport cost after the shift is calculated using the
        // north-west-corner rule and adding the  cost for creating and
        // destructing the mass
        // note:
        //  This is not the real transport cost in case of an unbalanced input
        //  since the mass that is added by constructing a balanced problem has not been taken into
        //  account.
        transportCost = Nw_Corner_Rule(costMatrix, supply, demand, iList, jList, weightList);

        transportCost +=createdMass*constructionCost+destructedMass*destructionCost;


        // if the supply and demand vector are both zero and the new transport cost is smaller
        // than the previous transport cost, the new transport plan and cost is optimal
        if(weightList[0] == -1){

            if(transportCost < transportCostPrev){
                transportCostPrev = transportCost;
                weightListPrev.erase(0, weightListPrev.length());
                iListPrev.erase(0, iListPrev.length());
                jListPrev.erase(0, jListPrev.length());

            }
            break;

        }

    }

    // take into account the mass added by constructing a balanced problem
    transportCostPrev -= destructionCost*difDestructMass + constructionCost*difCreateMass;

    exportVec = supplyTrue - supplyPrev;
    importVec = demandTrue - demandPrev;

return Rcpp::List::create(
    Rcpp::Named("transportCost") = transportCostPrev,
    Rcpp::Named("iList") = iListPrev,
    Rcpp::Named("jList") = jListPrev,
    Rcpp::Named("weightList") = weightListPrev,
    Rcpp::Named("import") = importVec,
    Rcpp::Named("export") = exportVec);

}
