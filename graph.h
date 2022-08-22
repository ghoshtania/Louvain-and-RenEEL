// File: graph_binary.h
// -- graph handling header file
//-----------------------------------------------------------------------------
// Community detection 
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// This program must not be distributed without agreement of the above mentionned authors.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time	    : February 2008
//-----------------------------------------------------------------------------
// see readme.txt for more details

#ifndef GRAPH_H
#define GRAPH_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>

#define WEIGHTED   0
#define UNWEIGHTED 1

using namespace std;

class Graph {
public:
  class Node
  {
    public:
    vector<int> actual_nodes;
  };
  
  unsigned int nb_nodes;
  unsigned int nb_links;
  double total_weight;  
  // unsigned int com;
  //double *nb_selfloop; //an array(of N elements) of which  return the number of self loops of the node
  //double *weighted_degree; //an array (of N elements) of which each element returns the weighted degree of the node

  vector<unsigned int> degrees;
  vector<vector<int> > nodes;
  vector<unsigned int> links;
  vector<double> weights;

  Graph();
  Graph(int size, int links);

  // binary file format is
  // 4 bytes for the number of nodes in the graph
  // 8*(nb_nodes) bytes for the cumulative degree for each node:
  //    deg(0)=degrees[0]
  //    deg(k)=degrees[k]-degrees[k-1]
  // 4*(sum_degrees) bytes for the links
  // IF WEIGHTED 4*(sum_degrees) bytes for the weights in a separate file
  //Graph(char *filename, char *filename_w, int type);
  Graph(vector<vector<int> >& c_nodes);
  Graph(int nb_nodes, int nb_links, double total_weight, int *degrees, int *links, double *weights);
  
  void add_node(vector<int>& n);
  void display(void);
  void display_for_check(void);
  void display_reverse(void);
 // void display_binary(char *outfile);
  void write_community(const char *filename_w);
  bool check_symmetry();
  //double nb_selfloops(unsigned int node);
  //double calculate_weighted_degree(unsigned int node);
//inline unsigned int cumulative_degree(unsigned int i);
  // return the number of neighbors (degree) of the node
  inline unsigned int nb_neighbors(unsigned int i);

  //return the number of self loops of the node
  inline double nb_selfloops(unsigned int node);

  //return the weighted degree of the node
  inline double weighted_degree(unsigned int node);
  //inline int Links(int i);
  //inline double get_Total_Weight(void);

  // return pointers to the first neighbor and first weight of the node
  inline pair<vector<unsigned int>::iterator, vector<double>::iterator > neighbors(unsigned int node);
};
 


inline unsigned int
Graph::nb_neighbors(unsigned int node) {
  assert(node>=0 && node<nb_nodes);

  if (node==0)
    return degrees[0];
  else
    return degrees[node]-degrees[node-1];
}



inline double
Graph::nb_selfloops(unsigned int node) {
  assert(node>=0 && node<nb_nodes);
  double self = 0;
  //cerr<<"HIIII"<<endl;
  pair<vector<unsigned int>::iterator, vector<double>::iterator > p = neighbors(node);
  for (unsigned int i=0 ; i<nb_neighbors(node) ; i++) {
    if (*(p.first+i)==node) {
    //cerr<<"Hello"<<endl;
      if (weights.size()!=0)
      { self = (double)*(p.second+i);
  return self; }
      else 
  {self = 1.;
    return self;}
    //cout<<"self"<<self<<endl;
    }
  
} 
  //cerr<<"HIIII"<<endl;
  return 0.;
}

    


inline double
Graph::weighted_degree(unsigned int node) {

  assert(node>=0 && node<nb_nodes);
    if (weights.size()==0){
        return (double)nb_neighbors(node);}
  else {
    pair<vector<unsigned int>::iterator, vector<double>::iterator > p = neighbors(node);
    // cout<< "Hey2"<<endl;
    double res = 0;
    for (unsigned int i=0 ; i<nb_neighbors(node) ; i++) {
      res += (double)*(p.second+i);
    }
    cout<<"weighted_degree"<<res<<endl; 
    return res;
  }
    
}

inline pair<vector<unsigned int>::iterator, vector<double>::iterator >
Graph::neighbors(unsigned int node) { //return pointers to the first neighbor and first weight of the node

  assert(node>=0 && node<nb_nodes);


  if (node==0)
    return make_pair(links.begin(), weights.begin());
  else if (weights.size()!=0)
    return make_pair(links.begin()+degrees[node-1], weights.begin()+degrees[node-1]);
  else

    return make_pair(links.begin()+degrees[node-1], weights.begin());
}


#endif // GRAPH_H
