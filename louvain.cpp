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
int reservehistory_of_node(int node);

struct part Louvain(Graph graph)
{
    struct part ans;
    
    Community c(graph, -1);
    // cout<<"Hi2"<<endl;
    Graph g;
    bool improvement = true;
    double mod=c.modularity(), new_mod;
    int level=0;
    cerr<<"old mod"<<mod<<endl;

    do  {
        {cerr << "level " << level << ":\n";
        //display_time("  start computation");
        cerr << "  network size: "
         << c.g.nb_nodes << " nodes, "
         << c.g.nb_links << " links, "
         << c.g.total_weight << " weight." << endl;
      }
      cerr<<"Hey"<<endl;
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
      

      
      cerr << "  modularity increased from " << mod << " to " << new_mod << endl;

      mod=new_mod;
      // if (verbose)
      //   display_time("  end computation");

      // if (filename_part!=NULL && level==1) // do at least one more computation if partition is provided
      //   improvement=true;
    } while(improvement);
    //cout<<"nb_nodes"<<c.g.nb_nodes;
    ans.pa=(int*)malloc(sizeof(int)*35);
    ans.Q = new_mod;
    cout<<"Final Mod"<<ans.Q<<endl;
    ans.com =c.g.nb_nodes;
    cout<<ans.com<< "number of communities"<<endl;
    //cout<<"Partition"<<c.partitioning(5);
    for(int i = 0; i<35; i++)
    {
        ans.pa[i] = reservehistory_of_node(i);
    }
    // for(int i = 0; i<35; i++)
    //     cout<<"partition"<<ans.pa[i]<<endl;
    return (ans);
}

int reservehistory_of_node(int node)
{
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

  vector<int> n2c(levels[0].size());
 for (unsigned int i=0 ; i<levels[0].size() ; i++)
    {
      n2c[i]=i;
    }
    int highestlevel = levels.size();
    //cout<<"Highestlevel"<<highestlevel<<endl;
    for (l=0 ; l<highestlevel ; l++)
      for (unsigned int node=0 ; node<levels[0].size() ; node++)
  n2c[node] = levels[l][n2c[node]];
  //remove("Com.txt");
     return (n2c[node]);
//     for (unsigned int node=0 ; node<levels[0].size() ; node++)
//       cout << node << " " << n2c[node] << endl;
}

