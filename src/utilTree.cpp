#include <vector>
#include <tuple>


struct Tree{

    // Every node in the tree is assigned a vector. These vectors hold pairs of values.
    // The first value is the key of node it is connected to with an edge.
    // The second value is the weight of that edge.
    std::vector<std::vector<std::pair<int, double>>> G;
    // the number of nodes in the tree
    int n;

    // the constructor sets the size of the tree and initializes the tree structure.
    Tree(int size){
        n = size;
        G.assign(n, std::vector<std::pair<int, double>>(0));
    }

    // Every new edge adds a new pair to two vectors. One pair for each of the nodes
    // the edge connects.
    void add_edge(int a, int b, double c){
        G[a].emplace_back(b, c);
        G[b].emplace_back(a, c);
    }

};

