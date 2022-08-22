// File: community.h
// -- community detection source file
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

#include "community.h"

using namespace std;

// Community::Community(char * filename, char * filename_w, int type, int nbp, double minm) {
//   g = Graph(filename, filename_w, type);
//   size = g.nb_nodes;   //size is the number of nodes present in the network
    
//   neigh_weight.resize(size,-1);
//   neigh_pos.resize(size);
//   neigh_last=0;
    

//   n2c.resize(size);
//   ini.resize(size);
//   tot.resize(size);

//   for (int i=0 ; i<size ; i++) {
    
//     n2c[i] = i;
//     tot[i] = g.weighted_degree(i); //
//     // cout<<"Hello1"<<endl;
//     ini[i]  = g.nb_selfloops(i);
//     // cout<<"Hello2"<<endl;
//       }
//     // cout<<"Hello"<<endl;
//   nb_pass = nbp;
//   min_modularity = minm;
//     //cout<<"Hello11"<<endl;
// }

Community::Community(struct Graph gc, int nbp) {
  g = gc; //
  size = g.nb_nodes;
  cout<<"size"<<size;
  neigh_weight.resize(size,-1);
  neigh_pos.resize(size);
  neigh_last=0;


  n2c.resize(size);
  ini.resize(size);
  cout<<"ini.size"<<ini.size()<<endl;
  tot.resize(size);

  for (int i=0 ; i<size ; i++) {
  //cout<<"check3"<<endl;
    n2c[i] = i;
    ini[i]  = g.nb_selfloops(i);
    tot[i] = g.weighted_degree(i);

  }
 
  nb_pass = nbp;
  
}

// void
// Community::init_partition(char * filename) {
//   ifstream finput;
//   finput.open(filename,fstream::in);

//   // read partition
//   while (!finput.eof()) {
//     unsigned int node, comm;
//     finput >> node >> comm;
    
//     if (finput) {
//       int old_comm = n2c[node];
//       neigh_comm(node); //produces all communities of that node

//       remove(node, old_comm, neigh_weight[old_comm]); //

//       unsigned int i=0;
//       for ( i=0 ; i<neigh_last ; i++) {
// 	unsigned int best_comm     = neigh_pos[i];
// 	float best_nblinks  = neigh_weight[neigh_pos[i]];
// 	if (best_comm==comm) {
// 	  insert(node, best_comm, best_nblinks);
// 	  break;
// 	}
//       }
//       if (i==neigh_last)
// 	insert(node, comm, 0);
//     }
//   }
//   finput.close();
// }

// inline void
// Community::remove(int node, int comm, double dnodecomm) {
//   assert(node>=0 && node<size);

//   tot[comm] -= g.weighted_degree(node);
//   in[comm]  -= 2*dnodecomm + g.nb_selfloops(node);
//   n2c[node]  = -1;
// }

// inline void
// Community::insert(int node, int comm, double dnodecomm) {
//   assert(node>=0 && node<size);

//   tot[comm] += g.weighted_degree(node);
//   in[comm]  += 2*dnodecomm + g.nb_selfloops(node);
//   n2c[node]=comm;
// }

void
Community::display() {
  for (int i=0 ; i<size ; i++)
    cerr << " " << i << "/" << n2c[i] << "/" << ini[i] << "/" << tot[i] ;
  cerr << endl;
}


double
Community::modularity() {
  double q  = 0.;
  double m2 = (double)g.total_weight;
  //cout <<"m2 value "<<m2<<endl;

  for (int i=0 ; i<size ; i++) {
    if (tot[i]>0)
      q += (double)ini[i]/m2 - ((double)tot[i]/m2)*((double)tot[i]/m2);
  }
 //cout<<"check"<<endl;
  return q;
}

void
Community::neigh_comm(unsigned int node) { //compute all neighbour communities of a particular node
  for (unsigned int i=0 ; i<neigh_last ; i++)
    neigh_weight[neigh_pos[i]]=-1;
  neigh_last=0;

  pair<vector<unsigned int>::iterator, vector<double>::iterator> p = g.neighbors(node);

  unsigned int deg = g.nb_neighbors(node);

  neigh_pos[0]=n2c[node];
  neigh_weight[neigh_pos[0]]=0;
  neigh_last=1;

  for (unsigned int i=0 ; i<deg ; i++) {
    unsigned int neigh        = *(p.first+i);
    unsigned int neigh_comm   = n2c[neigh];
    double neigh_w = (g.weights.size()==0)?1.:*(p.second+i);
    
    if (neigh!=node) {
      if (neigh_weight[neigh_comm]==-1) {
	    neigh_weight[neigh_comm]=0.;
	    neigh_pos[neigh_last++]=neigh_comm;
      }
      neigh_weight[neigh_comm]+=neigh_w;
    }
  }
}

void
Community::partition2graph() {
  vector<int> renumber(size, -1);
  for (int node=0 ; node<size ; node++) {
    renumber[n2c[node]]++;
  }

  int final=0;
  for (int i=0 ; i<size ; i++)
    if (renumber[i]!=-1)
      renumber[i]=final++;


  for (int i=0 ; i<size ; i++) {
    pair<vector<unsigned int>::iterator, vector<double>::iterator> p = g.neighbors(i);

    int deg = g.nb_neighbors(i);
    for (int j=0 ; j<deg ; j++) {
      int neigh = *(p.first+j);
      cout << renumber[n2c[i]] << " " << renumber[n2c[neigh]] << endl;
    }
  }
}

void
Community::display_partition() {
  vector<int> renumber(size, -1);
  for (int node=0 ; node<size ; node++) {
    renumber[n2c[node]]++;
  }

  int final=0;
  for (int i=0 ; i<size ; i++)
    if (renumber[i]!=-1)
      renumber[i]=final++;

  
  ofstream outfile1;
  outfile1.open("Com.txt",ofstream::app);
  //outfile1<<"1234 1234"<<endl;
  for (int i=0 ; i<size ; i++)
    outfile1 << i << " " << renumber[n2c[i]] << endl;
  outfile1.close();
}

double
Community :: partitioning(int which_node)
{
    {
      vector<int> renumber1(size, -1);
      for (int node=0 ; node<size ; node++) {
        renumber1[n2c[node]]++;
      }

      int final=0;
      for (int i=0 ; i<size ; i++)
        if (renumber1[i]!=-1)
          renumber1[i]=final++;

      for (int i=0 ; i<size ; i++)
          //cout<<"Here it is" << renumber1[n2c[i]]<<endl;
         return (renumber1[which_node]) ;
    }
}
Graph
Community::partition2graph_binary() {
  // Renumber communities
  vector<int> renumber(size, -1);
  for (int node=0 ; node<size ; node++) {
    renumber[n2c[node]]++;
  }

  int final=0;
  for (int i=0 ; i<size ; i++)
    if (renumber[i]!=-1)
      renumber[i]=final++;

  // Compute communities
  vector<vector<int> > comm_nodes(final);
  vector<vector<int> > communities(final);
  //printf("%s %d \n", __FILE__, __LINE__);
  for (int node=0 ; node<size ; node++) {
    //TODO add node handling
    vector<int>& comm = communities[renumber[n2c[node]]];  
    comm.insert(comm.end(), g.nodes[node].begin(), g.nodes[node].end());
    comm_nodes[renumber[n2c[node]]].push_back(node);
  }
  /*
  for (size_t i=0; i<g.nodes.size(); i++) {
    for (size_t j=0; j<g.nodes[i].size(); j++) {
      printf("%d n ", g.nodes[i][j]);
    }
    printf("n \n");
  }
  */
  /*
  for (size_t i=0; i<communities.size(); i++) {
    for (size_t j=0; j<communities[i].size(); j++) {
      printf("%d c ", communities[i][j]);
    }
    printf("i: %d \n", i);
  }
  printf("%s %d \n", __FILE__, __LINE__);
  */
  Graph g2(communities);  
  
  g2.nb_nodes = comm_nodes.size();  //number of nodes in graph2 - size of comm_nodes = 5
  g2.degrees.resize(comm_nodes.size());

  int comm_deg = comm_nodes.size(); //comm_deg = 5 how many supernodes are forming
  //cout<<"Howmany supernode"<<comm_deg<<endl;
  for (int comm=0 ; comm<comm_deg ; comm++) { //for each supernode ---1st loop
    map<int,double> m;
    map<int,double>::iterator it;

    int comm_size = comm_nodes[comm].size(); //size of 1st supernode
    //cout<<"supernodesize"<<comm_size<<endl;
    for (int node=0 ; node<comm_size ; node++) { //for each members in a supernode  ---nested loop
      pair<vector<unsigned int>::iterator, vector<double>::iterator> p = g.neighbors(comm_nodes[comm][node]);
      int deg = g.nb_neighbors(comm_nodes[comm][node]);
      for (int i=0 ; i<deg ; i++) {
    int neigh        = *(p.first+i);
    int neigh_comm   = renumber[n2c[neigh]];
    double neigh_weight = (g.weights.size()==0)?1.:*(p.second+i);

    it = m.find(neigh_comm);
    if (it==m.end())
      m.insert(make_pair(neigh_comm, neigh_weight));
    else
      it->second+=neigh_weight;
      }
    }
    g2.degrees[comm]=(comm==0)?m.size():g2.degrees[comm-1]+m.size(); //gives the cum degrees of each super node
    g2.nb_links+=m.size();
    //cout<<"degree of each supernode"<<g2.degrees[comm]<<endl;

    
    for (it = m.begin() ; it!=m.end() ; it++) {
      g2.total_weight  += it->second;
      g2.links.push_back(it->first);
      g2.weights.push_back(it->second); //weights have been modified
    }
    // cout<<"m.size"<<m.size()<<endl;
    // for(int node= 0; node<g2.links.size(); node++)
    //     {cout<<"node"<<node<<"g2.links"<<g2.links[node]; 
    //      cout<<" "<<endl;}
    // for(int i= 0; i<g2.weights.size(); i++)
    //     cout<<"node"<<i<<"g2.weights"<<g2.weights[i]<<endl; //g2.weights [node] gives weights of each member of a supernode
    // cout<<"g2.total_weight"<<g2.total_weight<<endl; //g2.total_weight gives the total weight of a supernode
    // cout<<"end of one supernode"<<endl;
  }

  return g2;
}


bool
Community::one_level() {
  bool improvement=false ;
  int nb_moves;
  int nb_pass_done = 0;
  double new_mod   = modularity();
  double cur_mod   = new_mod;

  vector<int> random_order(size);
  for (int i=0 ; i<size ; i++)
    random_order[i]=i;
  for (int i=0 ; i<size-1 ; i++) {
    int rand_pos = rand()%(size-i)+i;
    int tmp      = random_order[i];
    random_order[i] = random_order[rand_pos];
    random_order[rand_pos] = tmp;
  }
  // repeat while 
  //   there is an improvement of modularity
  //   or there is an improvement of modularity greater than a given epsilon 
  //   or a predefined number of pass have been done
  //cerr<<"size"<<size<<endl;
  do {
    cur_mod = new_mod;
    nb_moves = 0;
    nb_pass_done++;
    // for each node: remove the node from its community and insert it in the best community
    for (int node_tmp=0 ; node_tmp<size ; node_tmp++) {
//      int node = node_tmp;
      int node = random_order[node_tmp];
      int node_comm     = n2c[node];
      double w_degree = g.weighted_degree(node);


      // computation of all neighboring communities of current node
      neigh_comm(node);
      // remove node from its current community
      remove(node, node_comm, neigh_weight[node_comm]);
      //cerr<<"hey2"<<endl;

      // compute the nearest community for node
      // default choice for future insertion is the former community
      int best_comm        = node_comm;
      double best_nblinks  = 0.;
      double best_increase = 0.;
      for (unsigned int i=0 ; i<neigh_last ; i++) {
        double increase = modularity_gain(node, neigh_pos[i], neigh_weight[neigh_pos[i]], w_degree);
        //cerr<<"hey3"<<endl;
        if (increase>best_increase) {
          best_comm     = neigh_pos[i];
          best_nblinks  = neigh_weight[neigh_pos[i]];
          best_increase = increase;
        }
        
      }
      // insert node in the nearest community
      insert(node, best_comm, best_nblinks);
      if (best_comm!=node_comm)
        {nb_moves++;}
          //cerr<<"hey4"<<endl;}
    }
    // cerr<<"hey5"<<endl;
    // cerr<<"nb_moves"<<nb_moves<<endl;
    double total_tot=0;
    double total_in=0;
    //cerr<<"tot.size1"<<tot.size()<<endl;
    for (unsigned int i=0 ; i<tot.size() ;i++) {
      total_tot+=tot[i];
      total_in+=ini[i];
    }
    //cerr<<"tot.size2"<<tot.size()<<endl;
   
    new_mod = modularity();
    //cerr<<"modularity"<<new_mod<<endl;
    if (nb_moves>0)
      improvement=true;
    
  } while (nb_moves>0 && new_mod-cur_mod>.00001);

  return improvement;
}

