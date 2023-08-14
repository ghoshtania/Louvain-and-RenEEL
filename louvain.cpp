#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string>
#include <iostream> 
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "graph.h"
#include "community.h"
#include "louvain.h"

int display_level = -2;
using namespace std;
//int reservehistory_of_node(int node);

struct part Louvain(Graph graph,struct inserting I)
{
    struct part ans;
    //struct inserting Ia;
    //gcopy(I,&Ia);
    cout<<"3I.N "<<I.N<<"I.n "<<I.n<<"I.com "<<I.com<<endl;
    
    Community c(graph, -1);

    Graph g;
    bool improvement = true;
    double mod=c.modularity(), new_mod;
    int level=0;

    do  {
      //   {cerr << "level " << level << ":\n";
      //   //display_time("  start computation");
      //   cerr << "  network size: "
      //    << c.g.nb_nodes << " nodes, "
      //    << c.g.nb_links << " links, "
      //    << c.g.total_weight << " weight." << endl;
      // }
      //cerr<<"Hey"<<endl;
      improvement = c.one_level();
      
      new_mod = c.modularity();

      // if (++level==display_level)
      //   {cerr<<"true"<<endl;
      //   g.display();}
     // display_level = -1;
     //  if (display_level==-1)
      c.display_partition(); //show the partition
      
      g = c.partition2graph_binary();
      //g.display_for_check();
    
      c = Community(g, -1);
      

      
      //cerr << "  modularity increased from " << mod << " to " << new_mod << endl;

      mod=new_mod;
      // if (verbose)
      //   display_time("  end computation");

      // if (filename_part!=NULL && level==1) // do at least one more computation if partition is provided
      //   improvement=true;
    } while(improvement);
  vector<vector<int> >levels;
  ifstream finput;
  finput.open("Com.txt",fstream::in); //open file for reading

  int l=-1;
  while (!finput.eof()) {
    int node, nodecomm;
    finput >> node >> nodecomm;
    //cout<<"node"<<node<<" "<<"nodecomm"<<nodecomm<<endl;
    //cout<<"end"<<endl;
    if (finput) {
      if (node==0) {
         l++;
         levels.resize(l+1);
      }
      levels[l].push_back(nodecomm);
    }
  }
  //cout<<"nodes "<<levels[0].size()<<endl;
  vector<int> n2c(levels[0].size());
  vector<int>partition(levels[0].size()-1);
  vector<int>check;
  for (unsigned int i=0 ; i<levels[0].size() ; i++)
    {
      n2c[i]=i;
    }
    int highestlevel = levels.size();
    cout<<"Highestlevel"<<highestlevel<<endl;
    ans.pa=(int*)malloc(sizeof(int)*I.N);  //here levels[0].size() is 35 (where graph has 34 nodes total)
    ans.Q = new_mod;
    cout<<"Final Mod"<<ans.Q<<endl;
    ans.com =c.g.nb_nodes-1;
    cout<<ans.com<< "number of communities"<<endl;
    for (l=0 ; l<highestlevel ; l++)
      for (unsigned int node=0 ; node<levels[0].size() ; node++)  //n2c looks like - Node   Com
      {                                                           //                  0     0
        n2c[node] = levels[l][n2c[node]];                         //                  1     1
        //cout<<"node "<<node<< "n2c "<<n2c[node]<<endl;
                                                                 //                  ..    ..
      }                                                           //                  34    3 (total 3 communities)
    for (unsigned int node=0 ; node<levels[0].size() ; node++)
      partition[node] = n2c[node+1];

    // cout<<"Node"<<' '<<"Partition in Louvain"<<endl; 
    // for(int i = 0; i<34; i++)                                     //ans.pa looks like - Node   Com
    //     cout<<i<<' '<<ans.pa[i]<<endl;  
    int s =0;                                                        //                    0     1
    for (unsigned int node=0 ; node<I.com ; node++)
      {
        //cout<<"node "<<node<<"partition in vector "<<partition[node]<<endl;
        // ans.pa[node] = partition[node];
        // cout<<"node "<<node<<"ans.pa "<<ans.pa[node]<<endl;
        s+=I.ml[node][0];
        check.resize(s);
        //cout<<"check.size"<<check.size()<<endl;
      }
    //cout<<"check.finalsize"<<check.size()<<endl;
    for (int i=0; i<35; i++) check[i] = 0;
    for (int node=0; node<I.com; node++)
      { //if (!check[i]);
          for (int j=1; j<I.ml[node][0]+1; j++) {//cout<<"node "<<node<<"ml element"<<Ia.ml[node][j]<<endl;
          check[I.ml[node][j]-1] = partition[node];
          ans.pa[I.ml[node][j]-1] = partition[node];
        }
      }
    // for (int i =0; i<34; i++) cout<<"node "<<i<<"checkpartition "<<check[i]<<endl;  
    // for (int i =0; i<34; i++) cout<<"node "<<i<<"ans.partition "<<ans.pa[i]<<endl;
    //cout<<"Final size of check"<<check.size()<<endl;
    //freeG(Ia);
    return (ans);                                                 //                    ..    ..
}                                                                 //                   33     3 (total 3 communities)







// int reservehistory_of_node(int node)
// {
//    vector<vector<int> >levels;

//   ifstream finput;
//   finput.open("Com.txt",fstream::in); //open file for reading

//   int l=-1;
//   while (!finput.eof()) {
//     int node, nodecomm;
//     finput >> node >> nodecomm;
//     cout<<"node"<<node<<" "<<"nodecomm"<<nodecomm<<endl;
//     cout<<"end"<<endl;
//     if (finput) {
//       if (node==0) {
//          l++;
//          levels.resize(l+1);
//       }
//       levels[l].push_back(nodecomm);
//     }
//   }
//  cout<<"nodes "<<levels[0].size()<<endl;
//  // cout<<"level tania :" << endl;
//  //  for (l=0 ; l<levels.size() ; l++){
    
//  //    for (unsigned int node=0 ; node<levels[0].size() ; node++){
//  //            printf("%d ",levels[l][node]);
//  //      }
//  //      cout << endl;
//  //  }




//   vector<int> n2c(levels[0].size());

//  for (unsigned int i=0 ; i<levels[0].size() ; i++)
//     {
//       n2c[i]=i;
//     }
//     int highestlevel = levels.size();
//     cout<<"Highestlevel"<<highestlevel<<endl;
//     for (l=0 ; l<highestlevel ; l++)
//       for (unsigned int node=0 ; node<levels[0].size() ; node++)
//   n2c[node] = levels[l][n2c[node]];
//     cout << node << " " << n2c[node] << endl;
  
//      return (n2c[node]);
// //     for (unsigned int node=0 ; node<levels[0].size() ; node++)
// }

