#include <bits/stdc++.h>

#include "Rcpp.h"
#include "utilTree.cpp"

// Nodes for the binary search tree (splay tree) that is used to store the linear function segments
struct node{
    using np = node*;

    // p: parent node
    // chr: right child node
    // chl: left child node
    np p, chl, chr;

    // val: slope value
    // add: value that is added to the slope by moving up the tree
    // len: length of this element
    // size: sum of the lengths of all elements in the subtree of this node
    // nb: number of elements in the subtree of this node
    double val, add, len, size, nb;

    // key: index of the node in the transport tree, where this segment node was generated
    // it is used to determine where mass is created or destructed
    int key;


    node(){}

    node(int key, double val, double add, double len, double size, double nb, np p, np chl, np chr) : key(key), val(val), add(add), len(len), size(size), nb(nb), p(p), chl(chl), chr(chr){}


    bool is_root(){
        return !p;
    }

    // splay tree functions
    void rot();
    np splay();



    // determines whether this node is the left or right child of its parent, or if it is the root node
    int pos(){
        if(p){
            if(p->chl == this) return -1;
            if(p->chr == this) return 1;
        }
        return 0;
    }
};


// used for storing the keys and lengths of segments that vanish to the left of right side
Rcpp::NumericVector keysLeft;
Rcpp::NumericVector keysRight;

Rcpp::NumericVector lengthLeft;
Rcpp::NumericVector lengthRight;

using np = node*;

node dummy(-1,0, 0, 0, 0, 0, NULL, NULL, NULL);

std::vector<node> nodes;

// iterator for the nodes vector
int it;

inline int key(np t){return t->key; }
inline double size(np t){ return t->size; }
inline int nb(np t){ return t->nb; }

// Creates a new leaf node with given values.
// Since it is a leaf node, length and size are equal
np newnode(int key, double val, double len, double add, np p){
    nodes[it] = node(key, val, add, len, len, 1, p, &dummy, &dummy);
    return &nodes[it++];
}


// Calculating the new slope of a node and updating the values of its subtree accordingly
void push(np t){

    t->val += t->add;
    t->chl->add += t->add;
    t->chr->add += t->add;
    t->add = 0;
}


// Recursively updates the length and number of elements in the subtree of a node
np update(np t){
    if(!t || t == &dummy) return t;
    t->size = size(t->chl) + size(t->chr) + t->len;
    t->nb = nb(t->chl) + nb(t->chr) + 1;
    return t;
}

// rotates the tree at a node depending on its position under its parent
void node::rot(){
    np q = p->p;
    int pps = p->pos();
    if(pps == -1) q->chl = this;
    if(pps == 1 ) q->chr = this;
    int ps = pos();
    if(ps == -1){
        p->chl = chr;
        chr->p = p;
        chr = p;
    } else if(ps == 1){
        p->chr = chl;
        chl->p = p;
        chl = p;
    }
    p->p = this;
    update(p);
    p = q;
    update(this);
    update(q);
}

//rotates the node to the root of the tree
np node::splay(){
    while(!is_root()){
        int ps = pos();
        int pps = p->pos();
        if(pps == 0){
            rot();
        } else if(ps == pps){
            p->rot();
            rot();
        } else {
            rot();
            rot();
        }
    }
    return this;
}


// Updating all slopes to the left of the node and rotating the tree such that
// the leftmost node is the new root
np le(np t){
    assert(t != &dummy);
    push(t);
    if(t->chl != &dummy) return le(t->chl);
    return t->splay();
}

// Updating all slopes to the right of the node and rotating the tree such that
// the rightmost node is the new root
np re(np t){
    assert(t != &dummy);
    push(t);
    if(t->chr != &dummy) return re(t->chr);
    return t->splay();
}



// Updating the segment node values in x < 0
//
// This function is used to calculate the step from e_{v,x} to t_{p,x} by subtracting
// the edge weight from all segment node values in x < 0. If a segment passes through
// x = 0, it is split in two parts and a new segment node is inserted in to the tree.
//
//
// @param t A segment node.
// @param k A numeric value. Typically, the distance from the leftmost element to 0.
// @param c A numeric value. The value that is subtracted from all segment node slopes in x < 0.
//      If any segment goes through x = 0, the segment is divided in two parts.
// @return The first segment node in x < 0
//
np find(np t, double k, double c){
    if(t == &dummy) return t;

    // Checks if all elements in the subtree of  t are in x < 0.
    // In that case, c must be subtracted from all segment nodes in the subree.
    if(t->size <= k){
        t->add -= c;
        return t;
    }
    push(t);
    // If k is smaller then the length of all elements in the left subtree of t,
    // only those nodes are used.
        if(k < size(t->chl)) return find(t->chl, k, c);

        // If the element in t goes through 0, it is split and a new segment nod
        // is added.
        else if(k < size(t->chl) + t->len){
            // length of the segment in x < 0
            double ck = k - size(t->chl);
            // Rotating t to the root of the segment tree
            np res = t->splay();
            // Create new node in x < 0
            np l = newnode(t->key, t->val, ck, -c, res);
            // Shorten the current node to its length in x > 0
            res->len -= ck;
            // Insert the new node l into the tree. Since, res is the first segment
            // in x > 0 and l is the first in x < 0, l is inserted as the left
            // child.
            l->chl = res->chl;
            res->chl->p = l;
            res->chl = l;
            // Since res is fully in x > 0, c has to be added.
            res->val += c;


            // since res is the root of the tree, all other segments in x > 0
            // are in its right subtree and c has to be added.
            res->chr->add += c;
            update(l);
            return update(res);
        }
        // in any other case, at least one more slope element is needed to reach 0.
        // t-> chr: the next element with a larger slope than t is its right child.
        // the new value for k is the remaining distance to cover when the elements
        //in the left subtree of t and t itself is not enough
        return find(t->chr, k - size(t->chl) - t->len, c);
}


// Structure to hold the segment tree and additional information.
// One is used for each node in the transport tree.
struct TH{
    np root;

    // ll: slope at the left edge with length infinity
    // rr: slope at the right edge with length infinity
    double ll, rr;


    // keys of the transport tree nodes where the elements with slope
    // ll and rr were generated
    int llkey;
    int rrkey;

    // coordinate and function value of the end of the leftmost segment with finite length
    // pos: x-coordinate
    // ini: function value at pos
    double pos;
    double ini;
    TH(){}
};


// Inserts a new node with given parameters
// @param p parent of the new node
// @param val slope value of the new node
// @param len length of the new node
// @param key key of the new node

np insert(np &t, np p, double val, double len, int key, double k){
    // if t is a dummy node, create a new node in its position
    if(t == &dummy){
        np cur = newnode(key, val, len, 0, p);
        t = cur;
        return cur->splay();
    }
    push(t);

    // Otherwise it is inserted depending on its slope. If its value is smaller than
    // the value of node t it is inserted to the left and otherwise to the right.
    if(val != 0 || t->val != 0){

        if(val <= t->val){
            return insert(t->chl, t, val, len, key,k);
        } else{ // if(x->val > t->val){
            return insert(t->chr, t, val, len, key,k);
        }
    }else{

        if(k <= size(t->chl) + t->len){
            return insert(t->chr, t, val, len, key,k);
        }else{
            return insert(t->chl, t, val, len, key,k);
        }

        return t;
    }
}




// The same as above but with a given node to insert instead of values.
np insert(np &t, np p, np x, double k){
    // if t is a dummy node, create a new node in its position
    if(t == &dummy){
        x->p = p;
        x->chl = x->chr = &dummy;
        t = x;
        return update(x)->splay();
    }
    push(t);

    // Otherwise it is inserted depending on its slope. If its value is smaller than
    // the value of node t it is inserted to the left and otherwise to the right.

    if(x->val != 0 || t->val != 0){

        if(x->val <= t->val){
            return insert(t->chl, t, x,k);
        } else{ // if(x->val > t->val){
            return insert(t->chr, t, x,k);
        }
    }else{
        if(k <= size(t->chl) + t->len){
            return insert(t->chr, t, x,k);
        }else{
            return insert(t->chl, t, x,k);
        }

        return t;
    }

}


// Merging the subtree of s into the subtree of t
// @param t A segment node
// @param s A segment node
void merge(np &t, np s, double tpos, double spos){
    if(s == &dummy){
        return;
    }
    push(s);
    merge(t, s->chl, tpos, spos);
    merge(t, s->chr, tpos, spos);
    t = insert(t, NULL, s, tpos);
}




// The main function of the algorithm.
//
// It calculates the segment trees for all nodes in the given transport tree.
//
// @param T The transport tree as 'Tree' structure.
// @param a The position of the functions minimum in the initial state.
// @param creation A numeric vector. Give the cost for creating mass in each node.
// @param destruction A numeric vector. Give the cost for creating mass in each node.
// @param dp A vector of segment trees in TH structure.
// @param v The key of the current transport tree node.
// @param p The key of the current transport key nodes parent node.
// @return None
//
void treegkr_dfs(Tree &T, std::vector<double> &a, Rcpp::NumericVector &creation, Rcpp::NumericVector &destruction, std::vector<TH> &dp, int v, int p){

    // Setting the slope and keys of the leftmost and rightmost segments. These have length
    // infinity and are not stored in the segment tree.
    dp[v].ll = -creation[v];
    dp[v].llkey = v;
    dp[v].rr = destruction[v];
    dp[v].rrkey = v;
    // Creating the segment between 0 and the initial minimum of the function.
    if(a[v] >= 0){
        dp[v].root = newnode(v, -creation[v], a[v], 0, NULL);
        // Since the only segment node is in x > 0, the position of the left end of
        // the leftmost segment is 0.
        dp[v].ini = a[v] * creation[v];
        dp[v].pos = 0;
    } else {
        dp[v].root = newnode(v, destruction[v], -a[v], 0, NULL);
        // The leftmost element ends at the minimum 0.
        dp[v].ini = 0;
        // and is at x = a[v].
        dp[v].pos = a[v];
    }

    // Recursively calculating the segment trees for all nodes in the subtree of
    // the current node v.
    // These trees are then merged with the current transport tree nodes segment tree.
    for(std::pair<int, double> it: T.G[v]){
        // i: key of new transport tree node
        // c: edge cost from the current node v to the new node i
        int i = it.first;
        double c = it.second;
        // If the new node is the node that was visited before the current node, it is skipped.
        // (Nodes a and b are linked by an edge. After going from a to b, there is an edge
        // back from b to a)
        if(i == p){
            continue;
        }

        // calculating the transport on the subtree of i
        treegkr_dfs(T, a, creation, destruction, dp, i, v);

        // Adding and subtracting the edge cost from v to i from all segment nodes
        // depending on its position regarding x = 0.
        dp[i].root = find(dp[i].root, -dp[i].pos, c);
        // The function value change at 0 only depends on the edge cost and distance from 0 to
        // the end of the leftmost segment
        dp[i].ini -= dp[i].pos * c;


        dp[i].ll -= c;
        dp[i].rr += c;

        // Merging the segment trees of the current node with the one of the new node:

        // Updating ini and pos
        dp[v].ini += dp[i].ini;
        dp[v].pos += dp[i].pos;


        if (dp[v].ll < dp[i].ll){
            dp[v].llkey = dp[i].llkey;
            dp[v].ll = dp[i].ll;
        }

        if(dp[v].rr > dp[i].rr){
            dp[v].rrkey = dp[i].rrkey;
            dp[v].rr = dp[i].rr;
        }


        // The segment tree with fewer segment nodes is merged into the
        // other tree. This reduces the number of operations.
        if(dp[v].root == &dummy){
            dp[v].root = dp[i].root;
        } else if(dp[i].root != &dummy){
            if(dp[v].root->nb < dp[i].root->nb){
                std::swap(dp[i].root, dp[v].root);
            }
            merge(dp[v].root, dp[i].root, dp[v].pos, dp[i].pos);
        }

        // The merged tree can be simplified:

        // All elements with a slope smaller than ll are removed.
        while(dp[v].root != &dummy){
            dp[v].root = le(dp[v].root);
            if(dp[v].root->val >= dp[v].ll){
                break;
            }
            // The keys and length of the removed segment nodes are saved.
            // These will be used to compute the import and export vectors.
            if(dp[v].root->len > 0){
                keysLeft.push_back(dp[v].root->key);
                lengthLeft.push_back(dp[v].root->len);
            }

            dp[v].ini += dp[v].root->val * dp[v].root->len;
            dp[v].pos += dp[v].root->len;
            dp[v].root = dp[v].root->chr;
            dp[v].root->p = NULL;
        }

        // If one of the deleted segments with slope smaller than ll
        // crossed x = 0 a new element with slope ll is added to the tree.
        if(dp[v].pos > 0){
            dp[v].root = insert(dp[v].root, NULL, dp[v].ll, dp[v].pos, dp[v].llkey, dp[v].pos);
            dp[v].ini -= dp[v].ll * dp[v].pos;
            dp[v].pos = 0;
        }


        // All elements with a slope larger than rr are removed.
        while(dp[v].root != &dummy){
            dp[v].root = re(dp[v].root);
            if(dp[v].root->val <= dp[v].rr){
                break;
            }
            // The keys and length of the removed segment nodes are saved.
            // These will be used to compute the import and export vectors.
            if(dp[v].root->len > 0){
                keysRight.push_back(dp[v].root->key);
                lengthRight.push_back(dp[v].root->len);
            }

            dp[v].root = dp[v].root->chl;
            dp[v].root->p = NULL;
        }

        // If pos + length of all elements in the segment tree is smaller than 0,
        // a new element with slope rr is inserted to cover the distance from the
        // end of the rightmost element to 0.
        if(dp[v].pos + size(dp[v].root) < 0){
            dp[v].root = insert(dp[v].root, NULL, dp[v].rr, -(dp[v].pos + size(dp[v].root)), dp[v].rrkey, dp[v].pos);
        }


        }
    }



// A function to print the segment tree structure to the console
// Used for Debugging.
// void printTree(np t){
//
//     if(t->chl != &dummy){
//         Rcpp::Rcout << "left \n";
//         printTree(t->chl);
//     }
//
//     if(t->len > 0){
//         Rcpp::Rcout << "value: " << t->val << "  length: " << t->len <<
//             "  key :" << 1+t->key << "  add :" << t->add << "  size :" << t->size <<
//                 "  nb :" << t->nb  << "\n";
//     }
//     if( t->chr != &dummy){
//         Rcpp::Rcout << "right \n";
//         printTree(t->chr);
//     }
//
//     Rcpp::Rcout << "up \n";
//
//
// }





// Calcualting the keys and length of all segments with length > 0 in the subtree
//  of t
//
// @param t A node in a segment tree
// @param keys A vector to store the keys
// @param length A vector to store the keys
// @return by reference: Key and length vectors
//
void getKeysLengths(np t, Rcpp::NumericVector &keys, Rcpp::NumericVector &lengths){

    if(t->chl != &dummy){
        getKeysLengths(t->chl, keys, lengths);
    }

    if(t->len > 0){
        keys.push_back(t->key);
        lengths.push_back(t->len);
    }
    if( t->chr != &dummy){
        getKeysLengths(t->chr, keys, lengths);
    }
}



// function to call the algorithm from R
// TreeAd: nx2 Matrix containing the data tree structure. Each row contains the keys of
// two nodes, that are connected by an edge
// costMatrix: the Trees cost matrix. Only the costs of the edges are used
// r: suplly vector
// s: demand vector
// creation: cost of creation of mass in the data tree nodes
// destruction: cost of destruction of mass in the data tree nodes

// The tree metric unbalanced optimal transport algorithm
//
// This function makes the unbalanced optimal transport algorithm for tree metrics
// accessible from R. It calculated the optimal transport cost and the import vector.
//
// @param tree A tree structure given in list form. Each entry in the list represents
//          an edge: (first node, second node, edge weight)
// @param supply The supply vector.
// @param demand The demand vector.
// @param creationCost A numeric vector giving the creation cost at each node.
// @param destructionCost A numeric vector giving the destruction cost at each node.
// @return A list containing the optimal transport cost and the import vector.
// @export
// [[Rcpp::export]]
Rcpp::List treegkr_Rcpp (Rcpp::List &tree, Rcpp::NumericVector &supply, Rcpp::NumericVector &demand,
                         Rcpp::NumericVector &creation, Rcpp::NumericVector &destruction){
    // Number of nodes
    int n = supply.length();

    Rcpp::NumericVector import (n);


    // creating tree structure
    Tree T(n);
    Rcpp::NumericVector currentEdge;
    for (int i = 0; i < n-1; i++){
        currentEdge = tree[i];
        T.add_edge((int)currentEdge[0]-1, (int) currentEdge[1]-1,  (double) currentEdge[2]);
    }


    nodes.resize(4 * n + 1);
    it = 0;
    std::vector<TH> dp(n);

    // Calculating the positions of the initial minimum for each transport tree node
    std::vector<double> a(n);
    for(int i = 0; i < n; i++){
        a[i] = demand[i] - supply[i];
    }

    // Calling the main function
    treegkr_dfs(T, a, creation, destruction, dp, 0, -1);

    double pos = dp[0].pos;

    // Calculating the import vector

    // Calling find with 0, does not change the values of the segment tree nodes,
    // but splits any segment that crosses x = 0. Therefore, one element in the
    // segment tree ends at x = 0.
    dp[0].root = find(dp[0].root, -pos, 0);
    Rcpp::NumericVector finalKeys;
    Rcpp::NumericVector finalLengths;

    getKeysLengths(dp[0].root, finalKeys, finalLengths);

    // The import of mass depends on the position of the segment tree nodes.
    // If a segment generated at a transport tree node with initial minimum position
    // in x > 0 ends up in x < 0 in the finale segment tree, it indicates destruction
    // of mass at its transport tree node.
    // A segment from x < 0 that ends up at x > 0 indicates mass creation.
    // The length of the segment determines the amount of mass that is created
    // or destructed.
    // The import vector is updated accordingly.
    for(int i = 0; i < finalKeys.length(); i++){
        pos += finalLengths[i];
        if(pos <= 0   && a[finalKeys[i]] < 0){
            import[finalKeys[i]] -= finalLengths[i];
        }else if(pos > 0 && a[finalKeys[i]] > 0){
            import[finalKeys[i]] += finalLengths[i];
        }


    }

    // The same holds for segments that were deleted from the tree because its
    // slope was larger than rr or smaller then ll.
    for(int i = 0; i < keysLeft.length(); i++){

        if(a[keysLeft[i]] < 0){
            import[keysLeft[i]] -= lengthLeft[i];
        }
    }

    for(int i = 0; i < keysRight.length(); i++){
        if(a[keysRight[i]] > 0){
            import[keysRight[i]] += lengthRight[i];
        }

    }


    pos = dp[0].pos;
    double ans = dp[0].ini;

    // If the segment tree is empty, the function value in ans is the optimal
    // transport cost
    if(dp[0].root == &dummy){

        return Rcpp::List::create(
            Rcpp::Named("import") = import,
            Rcpp::Named("cost") = ans);
    }


    // Otherwise, the function value at x = 0 has to be caluclated.

    // Starting at the leftmost segment tree node, all nodes are traversed in
    // order of their values until x = 0 is reached.
    while(1){

        // Rotate the leftmost node until it is the root node.
        dp[0].root = le(dp[0].root);

        // If this segment is long enough to reach zero, calculate the final
        // transport cost
        if(pos + dp[0].root->len >= 0){

            return Rcpp::List::create(
                Rcpp::Named("import") = import,
                Rcpp::Named("cost") = (ans - dp[0].root->val * pos));
        // If not, update the position and ans values
        }else {

            pos += dp[0].root->len;
            ans += dp[0].root->val * dp[0].root->len;

        }
        // Since the leftmost element is the root, its right child is the the
        // node with the second smallest slope value.
        dp[0].root = dp[0].root->chr;
    }

}
