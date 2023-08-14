#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include"Simbody.h"

#define WEIGHTED   0
#define UNWEIGHTED 1

using namespace std;
class Community {
  struct graph G;
  unsigned int nb_nodes;
  unsigned long nb_links;
  double total_weight;  

  vector<unsigned long> degrees;
  vector<vector<int> > nodes;
  vector<unsigned int> links;
  vector<float> weights;
  // constructors:
  // reads graph from file using graph constructor
  // type defined the weighted/unweighted status of the graph file
  Community (struct graph);

  // return the number of neighbors (degree) of the node
 int nb_neighbors(int node);

  // return the number of self loops of the node
  double nb_selfloops(int node);

  // return the weighted degree of the node
  double weighted_degree(int node);

  // return pointers to the first neighbor and first weight of the node
  pair<vector<int>::iterator, vector<float>::iterator > neighbors(int node);
};
