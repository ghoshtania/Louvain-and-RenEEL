// File: graph_binary.cpp
// -- graph handling source
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

#include <sys/mman.h>
#include <fstream>
#include "graph.h"
#include "math.h"
#include "body.h"



Graph::Graph() {
  nb_nodes     = 0;
  nb_links     = 0;
  total_weight = 0;
  denominator =  0;
  nb_com = 0;
}

Graph::Graph(int size, int link){
  nb_nodes = size+1;
  ///denominator = link;
  total_weight = 0;
  // com = size;
  //degrees.resize(nb_nodes);
 

  
  //nb_selfloop = (double*)malloc(sizeof(int)*size);
  //weighted_degree = (double*)malloc(sizeof(int)*size);
  // Compute total weight
  // for (unsigned int i=0 ; i<nb_nodes ; i++) {
  //   total_weight += (double)weighted_degree(i);
  // }
  for (unsigned int i=0 ; i<nb_nodes ; i++) {
    vector<int> n;
    n.push_back(i);
    nodes.push_back(n);

  }
  // for (unsigned int i=0 ; i<nb_nodes ; i++) {
  //   total_weight += (double)weighted_degree(i);
  // }
}

// Graph::Graph(char *filename, char *filename_w, int type) {
//   ifstream finput;
//   finput.open(filename,fstream::in | fstream::binary);

//   // Read number of nodes on 4 bytes
//   // finput.read((char *)&nb_nodes, 4);
//   // cout<<"Nodes"<<nb_nodes<<endl;
//   // assert(finput.rdstate() == ios::goodbit);

//   // // Read cumulative degree sequence: 8 bytes for each node
//   // // cum_degree[0]=degree(0); cum_degree[1]=degree(0)+degree(1), etc.
//   // degrees.resize(nb_nodes);
//   // finput.read((char *)&degrees[0], nb_nodes*8);
//   // cout<<"degree"<<degrees[2]<<endl;

//   // // Read links: 4 bytes for each link (each link is counted twice)
//   // nb_links=degrees[nb_nodes-1];
//   // links.resize(nb_links);
//   // finput.read((char *)(&links[0]), (long)nb_links*8);
//   // cout<<"links"<<links[1]<<endl;

//   // IF WEIGHTED : read weights: 4 bytes for each link (each link is counted twice)
//   weights.resize(0);
//   total_weight=0;
//   if (type==WEIGHTED) {
//     ifstream finput_w;
//     finput_w.open(filename_w,fstream::in | fstream::binary);
//     weights.resize(nb_links);
//     finput_w.read((char *)&weights[0], (long)nb_links*4);  
//   }

//   //initilize first graph without contraction
//   for (unsigned int i=0 ; i<nb_nodes ; i++) {
//     vector<int> n;
//     n.push_back(i);
//     nodes.push_back(n);
//   }
  
  
  
  
// }


Graph::Graph(int n, int m, double t, int *d, int *l, double *w) {
/*  nb_nodes     = n;
  nb_links     = m;
  total_weight = t;
  degrees      = d;
  links        = l;
  weights      = w;*/
}

Graph::Graph(vector<vector<int> >& c_nodes) {
  nb_nodes     = 0;
  nb_links     = 0;
  total_weight = 0;
  denominator = 0;
  nodes.reserve(c_nodes.size());
  for (size_t i=0; i<c_nodes.size(); i++) {
    nodes.push_back(c_nodes[i]);
  }
}

void 
Graph::add_node(vector<int>& n) {
    //Graph::Node node;
    //node.actual_nodes = n;
    nodes.push_back(n);
}
void
Graph::display_for_check()
{
  cout<<"No of Nodes"<<nb_nodes <<endl;
  cout<<"No of Links"<<nb_links <<endl;
  
  for (int i =0; i<nb_nodes; i++)
  {
    cout<<degrees[i]<<" ";
  }

  cout << endl;
  cout<<"total weight"<< total_weight<<endl;
  cout<<"nb_links"<<nb_links<<endl;
// for (int i =0; i<nb_nodes; i++)
//   {
//     for(int j=0;j<nb_nodes;j++){
//       cout << nl[i][j]<<"  ";
//     }
//     cout << endl;
//   }

// for (int i = 0; i<nb_links; i++)
// cout<<"Links no" << Links(i)<<endl;

// for (int i = 0; i<nb_nodes; i++)
//   cout<<"Weighteddeg"<<weighted_degree(i)<<endl;
  for (int i=0; i<nb_links; i++)
    {
      cout<<"Weights "<<weights[i]<<' '<<"Links "<<links[i]<<endl;
    }

  
}
void
Graph::display() {
  /*
  for (unsigned int node=0 ; node<nb_nodes ; node++) {
    pair<vector<unsigned int>::iterator, vector<float>::iterator > p = neighbors(node);
    cout << node << ":" ;
    for (unsigned int i=0 ; i<nb_neighbors(node) ; i++) {
      if (true) {
	if (weights.size()!=0)
	  cout << " (" << *(p.first+i) << " " << *(p.second+i) << ")";
	else
	  cout << " " << *(p.first+i);
      }
    }
    cout << endl;
  }
  */
  for (size_t i=0; i<nodes.size(); i++) {
    for (size_t j=0; j<nodes[i].size(); j++) {
      cout << nodes[i][j] << " ";  
    }
    cout << "\n";  
  }
}

void
Graph::write_community(const char *filename_w) {
  FILE* fOut = fopen(filename_w, "wt");
  for (size_t i=0; i<nodes.size(); i++) {
    for (size_t j=0; j<nodes[i].size(); j++) {
      fprintf(fOut, "%d ", nodes[i][j]);  
    }
    fprintf(fOut, "\n");  
  }
}

// void
// Graph::display_reverse() {
//   for (unsigned int node=0 ; node<nb_nodes ; node++) {
//     pair<vector<unsigned int>::iterator, vector<double>::iterator > p = neighbors(node);
//     for (unsigned int i=0 ; i<nb_neighbors(node) ; i++) {
//       if (node>*(p.first+i)) {
// 	if (weights.size()!=0)
// 	  cout << *(p.first+i) << " " << node << " " << *(p.second+i) << endl;
// 	else
// 	  cout << *(p.first+i) << " " << node << endl;
//       }
//     }   
//   }
// }


// bool
// Graph::check_symmetry() {
//   int error=0;
//   for (unsigned int node=0 ; node<nb_nodes ; node++) {
//     pair<vector<unsigned int>::iterator, vector<double>::iterator > p = neighbors(node);
//     for (unsigned int i=0 ; i<nb_neighbors(node) ; i++) {
//       unsigned int neigh = *(p.first+i);
//       double weight = *(p.second+i);
      
//       pair<vector<unsigned int>::iterator, vector<double>::iterator > p_neigh = neighbors(neigh);
//       for (unsigned int j=0 ; j<nb_neighbors(neigh) ; j++) {
// 	unsigned int neigh_neigh = *(p_neigh.first+j);
// 	double neigh_weight = *(p_neigh.second+j);

// 	if (node==neigh_neigh && weight!=neigh_weight) {
// 	  cout << node << " " << neigh << " " << weight << " " << neigh_weight << endl;
// 	  if (error++==10)
// 	    exit(0);
// 	}
//       }
//     }
//   }
//   return (error==0);
// }


// void
// Graph::display_binary(char *outfile) {
//   ofstream foutput;
//   foutput.open(outfile ,fstream::out | fstream::binary);

//   foutput.write((char *)(&nb_nodes),4);
//   foutput.write((char *)(&degrees[0]),4*nb_nodes);
//   foutput.write((char *)(&links[0]),8*nb_links);
// }
