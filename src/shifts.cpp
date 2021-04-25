#include "Rcpp.h"


// Calculating u(k) and v(k) for the south-west-shifts
void CalcSWShift(Rcpp::NumericVector costMatrix, Rcpp::NumericVector iList, Rcpp::NumericVector jList, Rcpp::NumericVector &ukSW, Rcpp::NumericVector &vkSW){

    // Largest index in supply and demand vector that holds mass
    int maxI = Rcpp::max(iList);
    int maxJ = Rcpp::max(jList);

    // remove all old values
    ukSW.erase(0, ukSW.length());
    vkSW.erase(0, vkSW.length());

    // first entry of v(k) is always 0
    vkSW.push_back(0);

    // stepwise u(k) and v(k) shift cost
    double cost = 0.0;

    // loop index
    int i = 0;

    // calculating ukSW
    // all values are computed while going through the transport path once
    while(i < jList.length() && jList[i] <= maxJ){

        // every element in the path where the j index increases is a new end point of
        // a "u(k) shift".
        while(i < jList.length()-1 && iList[i+1] == iList[i] && jList[i+1] <= maxJ){

            ukSW.push_back(cost-costMatrix(iList[i], jList[i]));
            i++;

            if(i == iList.length()){
                break;
            }
        }

        ukSW.push_back(cost-costMatrix(iList[i], jList[i]));

        // Update the cost if maxJ is not reached yet
        if(jList[i] != maxJ){

            cost -= costMatrix(iList[i],jList[i]);

            // find next south-west corner
            while(jList[i] == jList[i+1]){
                i++;

                if(i == jList.length()){
                    break;
                }

            }



            // increase cost by the cost of the south-west corner
            cost += costMatrix(iList[i+1],jList[i]);
            i++;

        }else{
            break;
        }

    }

    // resetting loop index and cost
    i = 0;
    cost = 0;

    // calculating vkSW
    // all values are computed while going through the transport path once
    while(i <= iList.length() && iList[i] < maxI){

        // increase loop index until the next north-east corner is found
        while(i < iList.length()-1 && iList[i+1] == iList[i]){
            i++;
        }
        // update cost according to the corners cost
        cost -= costMatrix(iList[i], jList[i]);

        // finding the next south-west corner
        // very new element in the column is an end point to the v(k) shift
        while(i < jList.length()-1 && jList[i] == jList[i+1] && iList[i+1] < maxI){
            i++;
            vkSW.push_back(cost + costMatrix(iList[i], jList[i]));

        }
        // If the next element differs in i and j, add a new v(k) shift value
        //  in the same column but new row, to keep the overall mass in the
        //  column equal to before the shift
        if(iList[i+1] != iList[i]){
            vkSW.push_back(cost+costMatrix(iList[i+1], jList[i]));

        }

        // update the cost and go to next element
        cost += costMatrix(iList[i+1], jList[i]);
        i++;

    }
    return;
}



// Calculating u(k) and v(k) for the north-east shifts
void CalcNOShift(Rcpp::NumericVector costMatrix, Rcpp::NumericVector iList, Rcpp::NumericVector jList, Rcpp::NumericVector &ukNO, Rcpp::NumericVector &vkNO){

    // Largest index in supply and demand vector that holds mass
    int maxI = max(iList);
    int maxJ = max(jList);

    // delete all old values
    ukNO.erase(0, ukNO.length());
    vkNO.erase(0, vkNO.length());

    // the first value in vk is always 0 and
    // the first in uk is always -(cost for first element in the transport path)
    vkNO.push_back(0);
    ukNO.push_back(-costMatrix(iList[0], jList[0]));

    // loop index
    int i = 0;

    // stepwise transport cost
    double cost = 0;


    // calculating ukNO
    // all values are computed while going through the transport path once

    // The first south-west corner is computed before the main loop starts.
    // Each new element in the path is the end point of a shift
    while(iList[i+1] <= maxI && jList[i+1] == jList[i]){
        i++;
        ukNO.push_back(-costMatrix(iList[i], jList[i]));
    }

    // update cost with cost of the first south-west corner
    cost = -costMatrix(iList[i], jList[i]);

    // repeat until the first instance of the larges i value is reached
    while(i < iList.length() && iList[i] <= maxI){

        // find the next north-east corner
        while(i < iList.length()-1 && iList[i] == iList[i+1]){
            i++;

        }

        if(i+1 >= jList.length()){
            break;
        }

        // update the shift cost with the cost of the north-east corner
        cost += costMatrix(iList[i], jList[i+1]);
        i++;


        // Find the next south-west corner
        // Every element in the column is the end point of a shift
        while(i < jList.length()-1 && jList[i] == jList[i+1] && iList[i+1] <= maxI){

            ukNO.push_back(cost-costMatrix(iList[i], jList[i+1]));
            i++;
        }

        // update shift cost using the cost of the corner
        cost -= costMatrix(iList[i], jList[i]);
        if(i < iList.length()){
            ukNO.push_back(cost);
        }

    }

    // resetting loop index and cost
    i = 0;
    cost = 0;

    // calculating vkNO
    // all values are computed while going through the transport path once
    while(i < jList.length()-1 && jList[i] < maxJ){

        // find the first south-west corner
        while(i < iList.length()-1 && jList[i] == jList[i+1]){
            i++;
        }

        // update shift cost using the cost of the south-west corner
        cost -= costMatrix(iList[i], jList[i]);

        // find the next north-east corner
        // every element until the corner gives a new value for vkNO
        while(i < iList.length()-1 && iList[i] == iList[i+1] && jList[i+1] <= maxJ){
            i++;
            vkNO.push_back(cost+costMatrix(iList[i], jList[i]));

        }

        // If the next element is in the same column just update the cost.
        // No v(k) value is added here.
        if(i < iList.length()-1 && jList[i] == jList[i+1]){
            cost += costMatrix(iList[i], jList[i]);

            // If the next value is not in the same column, update the cost and
            // add a new v(k) value
        }else if(i < iList.length()-1){
            vkNO.push_back(cost + costMatrix(iList[i], jList[i+1]));
            cost += costMatrix(iList[i], jList[i+1]);

            i++;

        }

    }

    return ;

}
