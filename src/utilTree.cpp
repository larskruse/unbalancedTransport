/*
MIT License

Copyright (c) 2020 joisino

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

//See https://github.com/joisino/treegkr


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

