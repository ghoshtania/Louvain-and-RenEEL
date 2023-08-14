#include <sys/mman.h>
#include <fstream>
#include "math.h"
#include "SimCommunity.h"
// #include"Simbody.h"

using namespace std;

Community::Community(graph G){

nb_nodes = G.N;
nb_links = G.n;


//creat links vector
links.resize(nb_links);
for (unsigned int i=0 ; i<nb_nodes ; i++)
	for (int j = 1; j<=G.nl[i][0];j++)
		links.push_back(G.nl[i][j]);

// Compute total weight
for (int i=0 ; i<nb_nodes ; i++) {
    total_weight += (double)weighted_degree[i];
  }
};

inline int
Community::nb_neighbors(int node){
	return G.nl[node][0];
}
inline double
Community ::nb_selfloops(int node){
	return 2*G.ml[node][0];
}
inline double
Community ::weighted_degree(int node){
	return G.dl[node];
}