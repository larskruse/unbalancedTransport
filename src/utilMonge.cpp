#include "Rcpp.h"
#include "utilMonge.h"

// C++ implantation of the North-West Corner Rule
double Nw_Corner_Rule(Rcpp::NumericMatrix costMatrix, Rcpp::NumericVector supply, Rcpp::NumericVector demand, Rcpp::NumericVector &iList, Rcpp::NumericVector &jList, Rcpp::NumericVector &weightList){


    double checkThreshold;
    double threasholdValue = 1e-10;

    double transport = 0.0;

    // Copying the original supply and demand vectors to preserve the original values
    Rcpp::NumericVector supplyAlt = Rcpp::clone(supply);
    Rcpp::NumericVector demandAlt = Rcpp::clone(demand);


    // Setting index variables to 0
    int i = 0;
    int j = 0;

    // Deleting the old vector entries
    weightList.erase(0, weightList.length());
    iList.erase(0, iList.length());
    jList.erase(0, jList.length());


    // If the supply and demand are equal to 0, no mass is transported
    if(Rcpp::is_true(Rcpp::all(supply == 0)) && Rcpp::is_true(Rcpp::all(demand == 0))){

        // Setting the first value in weightList to -1
        // This indicates that no mass is transported
        weightList.push_front(-1);
        return 0;

    }

    // Main loop of the function.
    // Repeat until no mass is left to transport, e.g. supply and demand are 0
    while(Rcpp::is_true(Rcpp::any(supplyAlt != 0)) || Rcpp::is_true(Rcpp::any(demandAlt != 0))){

        if(i == supplyAlt.length() || j == demandAlt.length()){
            break;
        }


        // Getting the next non zero entries of the supply and demand vectors
        while(supplyAlt[i] == 0){
            i++;
        }

        while(demandAlt[j] == 0){
            j++;
        }

        if(i >= supplyAlt.length() || j >= demandAlt.length()){
            break;
        }


        // If the current supply and demand mass are equal,
        // set both to zero and push the mass and indexes to the
        // vectors
        if(supplyAlt[i] == demandAlt[j]){


            weightList.push_back(supplyAlt[i]);
            iList.push_back(i);
            jList.push_back(j);
            supplyAlt[i] = 0;
            demandAlt[j] = 0;

            i++;
            j++;

            if(i == supplyAlt.length() || j == demandAlt.length()){

                break;
            }

        // If the current supply is larger than the demand,
        // push the position, and mass of the demand vector to the vectors.
        // Set the current demand mass to zero and update the supply mass.
        }else if(supplyAlt[i] > demandAlt[j]){

            weightList.push_back(demandAlt[j]);
            iList.push_back(i);
            jList.push_back(j);
            supplyAlt[i] -= demandAlt[j];
            checkThreshold = supplyAlt[i];

            if(std::abs(checkThreshold) <  threasholdValue){
                supplyAlt[i] = 0;
                i++;
            }
            demandAlt[j] = 0;
            j++;
            if(i == supplyAlt.length() || j == demandAlt.length()){
                break;
            }


        // If the current demand is larger than the supply,
        // push the position, and mass of the supply vector to the vectors.
        // Set the current supply mass to zero and update the demand mass.
        }else{
            weightList.push_back(supplyAlt[i]);
            iList.push_back(i);
            jList.push_back(j);
            demandAlt[j] -= supplyAlt[i];
            checkThreshold = demandAlt[j];
            if(std::abs(checkThreshold) <  threasholdValue){
                demandAlt[j] = 0;
                j++;
            }
            supplyAlt[i] = 0;
            i++;
            if(i == supplyAlt.length() || j == demandAlt.length()){
                break;
            }

        }

    }

    // The transport cost is calculated by multiplying all weights with the transport
    // cost at the corresponding position in the cost matrix.
    for(int k = 0; k < weightList.length(); k++){
        transport = transport + weightList[k] * costMatrix(iList[k], jList[k]);
    }


    return transport;

}





// Expanding an unbalanced to a balanced optimal transport problem

void UnbalancedToBalanced(Rcpp::NumericMatrix &costMatrix, Rcpp::NumericVector &supply, Rcpp::NumericVector &demand){


    // Mass difference between supply and demand vector
    double dif;

    // Maximum of the cost matrix
    double Kc = Rcpp::max(costMatrix);

    // Dimensions of the cost matrix
    int N = costMatrix.nrow();
    int M = costMatrix.ncol();

    // Next value to add to the cost matrix
    double newValue;

    // New row or column for the cost matrix
    Rcpp::NumericVector newRowCol = Rcpp::NumericVector::create(Kc);


    // If the supply and demand mass is equal, the problem is already balanced.
    // Otherwise, a new entry in the supply or demand vector and a corresponding
    // row or column in the cost matrix is added.
    if(Rcpp::sum(supply) != Rcpp::sum(demand)){


        if(Rcpp::sum(supply) > Rcpp::sum(demand)){

            // Adding missing mass to the demand vector
            dif = Rcpp::sum(supply)-Rcpp::sum(demand);
            demand.push_back(dif);



            // Calculate new cost matrix column according to the method described
            // in the paper

            for(int i = 0; i < (N-1); i++){
                // compute the new value such that the matrix maintains the Monge property
                newValue = costMatrix(N-i-2,M-1) + newRowCol[0] - costMatrix(N-i-1, M-1);



                // Adding the new value to the beginning of the new column
                newRowCol.push_front(std::max(newValue, Kc));
            }

            costMatrix = Rcpp::cbind(costMatrix,newRowCol);

        }else{

            // Adding missing mass to the supply vector
            dif = Rcpp::sum(demand) - Rcpp::sum(supply);
            supply.push_back(dif);

            // Calculate new cost matrix column according to the method described
            // in the paper
            costMatrix = transpose(costMatrix);
            for(int i = 0; i < (M-1); i++){
                // compute the new value such that the matrix maintains the Monge property
                newValue = costMatrix(M-i-2,N-1) + newRowCol[0] - costMatrix(M-i-1, N-1);
                newRowCol.push_front(std::max(newValue, Kc));
            }

            costMatrix = Rcpp::cbind(costMatrix, newRowCol);
            costMatrix = transpose(costMatrix);

        }

    }

}


// Calculating which shifts are possible and if they are in the north-east
// or south-west direction

void GetIndexGroups(Rcpp::NumericVector iList, Rcpp::NumericVector jList, Rcpp::List &indexGroupsSWi, Rcpp::List &indexGroupsSWj, Rcpp::List &indexGroupsNOi, Rcpp::List &indexGroupsNOj){

    // Indexes of the supply and demand vectors with positive mass
    Rcpp::NumericVector uniqueI = Rcpp::unique(iList).sort();
    Rcpp::NumericVector uniqueJ = Rcpp::unique(jList).sort();

    // Temp variables to hold all indexes that belong to the same shift group
    std::vector<double> groupSWi;
    Rcpp::NumericVector groupSWj;

    Rcpp::NumericVector groupNOi;
    std::vector<double> groupNOj;



    // creating new lists
    indexGroupsSWi = Rcpp::List::create();
    indexGroupsSWj = Rcpp::List::create();
    indexGroupsNOi = Rcpp::List::create();
    indexGroupsNOj = Rcpp::List::create();


    // Shifts between two points are only possible if iList has more than one element.
    if(iList.length() > 1){

        // computing indexes for south-west shifts
        for(int i = 0; i < iList.length()-1; i++){

            // all points that have the same j index, belong in the same shift group and can be
            // shifted to the j indexes larger than that.
            // The corresponding i indexes are stored in groupSWi
            if(jList[i] == jList[i+1]){
                groupSWi.push_back(iList[i]);


                // If the next element in the path has a different j index,
                // a shift from the indexes in groupSWi to all j indexes larger than their
                // j index is possible.
                // These j indexes are stored in groupSWj
                // both vectors are added to the corresponding lists
            }else if(groupSWi.size()>0){
                if(iList[i] != iList[i+1]){
                    groupSWi.push_back(iList[i]);
                }

                indexGroupsSWi.push_back(Rcpp::wrap(groupSWi));
                groupSWi.clear();

                groupSWj = uniqueJ[uniqueJ > jList[i]];
                indexGroupsSWj.push_back(Rcpp::wrap(groupSWj));
                groupSWj.erase(0, groupSWj.length());


                // If neither the i nor j index of the next element in the path
                // is equal to the indexes of the  current element,
                // it is added as a single point.
            }else if(iList[i] != iList[i+1]){
                groupSWi.push_back(iList[i]);
                indexGroupsSWi.push_back(Rcpp::wrap(groupSWi));
                groupSWi.clear();

                groupSWj = uniqueJ[uniqueJ > jList[i]];
                indexGroupsSWj.push_back(Rcpp::wrap(groupSWj));
                groupSWj.erase(0, groupSWj.length());


            }

            // Computing indexes for north-east shifts:
            // All points that have the same i index, belong in the same shift group
            // and can be shifted to all i indexes larger than theirs.
            if(iList[i] == iList[i+1]){
                groupNOj.push_back(jList[i]);

            }else if(groupNOj.size() > 0){
                if(jList[i] != jList[i+1]){
                    groupNOj.push_back(jList[i]);
                }

                indexGroupsNOj.push_back(Rcpp::wrap(groupNOj));
                groupNOj.clear();

                groupNOi = uniqueI[uniqueI > iList[i]];
                indexGroupsNOi.push_back(Rcpp::wrap(groupNOi));
                groupNOi.erase(0, groupNOi.length());

                // If neither the i nor j index of the next element in the
                // transport path is equal to the current indexes, it is added
                // as a single point.
            }else if(jList[i] != jList[i+1]){
                groupNOj.push_back(jList[i]);
                indexGroupsNOj.push_back(Rcpp::wrap(groupNOj));
                groupNOj.clear();

                groupNOi = uniqueI[uniqueI > iList[i]];
                indexGroupsNOi.push_back(Rcpp::wrap(groupNOi));
                groupNOi.erase(0, groupNOi.length());

            }

        }

    }


    return;

}
