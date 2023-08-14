//
//  main_visualization.cpp
//  
//
//  Created by Tania Ghosh on 8/31/22.
#include<cstdio>
#include<cstdlib>
#include <climits>
#include <ctime>
//#include"help.h"
#include"graph.h"
#include "louvain.h"

#define Tsmall 0.0000001

int *cols;//keeps track
double cumulative_degree(unsigned int i,struct insert I);
unsigned int Links(int i,struct insert I);
double Weight(unsigned int i,struct insert I);
Graph convert_Insert_to_class_Graph(struct insert I);
void inputGraph(struct insert *I);
// int cleansort(struct part *ensem,int kmax,int N); //help to sort the partitions from lowest to highest value of Q
// int getscore(int *score,int *edgecc, int N,int know,int *ref);
// void renorm(struct graph G, int *ref);
// int mergroup(struct graph *G, int x, int y);
// int comps(int *p, int *q, int N); //return false if elements of array p and array q are not same
// void remolast(struct part *ensemble, int ens, int *edgecc, int N) //to remove one set (last) of partition from the ensemble
// {
//     long i;
//     long j;
//     for (i=0;i<N;i++) //initialize an array cols of element 0
//         cols[i]=0;
    
//     for (i = 0;i<N-1;i++) //for each node
//         if (!cols[i])
//         {
//             cols[i]=1;
//             for (j = i+1;j<N;j++) //for each node except ith node
//             {
//                 if (ensemble[0].pa[i] == ensemble[0].pa[j]) //if i and jth node belongs to same community for ensemble 0
//                 edgecc[i*N + j]--; //reduces by 1
//                 if (edgecc[i*N+j]==ens-1) //when the nodes are present in all ensembles (ens-1), jth element of
//                 cols[j]=1;                //cols is 1.
//             }
        
//         }
//     //delete in ensemble
//     i = 0;
//     free(ensemble[0].pa);     //delete the worst partition - ensemble[0].pa
//     for (i = 0;i<ens - 1;i++)
//     {
//         ensemble[i] = ensemble[i + 1]; //replace ensemble 0 with ensemble 1 and so on
//     }
//     ensemble[ens - 1].pa = NULL;
// }

// void replaceone(struct part* ensemble,int know, struct part* inensem, int max, int *edgecc, int N)
// {              //to replace ensemble[0] with inensemble[max]
//     long i,j;
//     for (i=0;i<N;i++)    cols[i]=0;
//     for (i = 0;i<N-1;i++)
//         if (!cols[i])
//         {
//             cols[i]=1;
//             for (j = i+1;j<N;j++)
//             {
//                 if (inensem[max].pa[i] == inensem[max].pa[j]) //if i and jth node are in same partition in inensem[max],
//                 edgecc[i*N + j]++; //increase by 1
//                 if (ensemble[0].pa[i]==ensemble[0].pa[j]) //if i and jth node are in same partition in ensemble[0]
//                 edgecc[i*N+j]--;                          //decrease by 1
//                 if (edgecc[i*N+j]==know)    //when the nodes are present in all ensembles, jth element of cols is 1
//                 cols[j]=1;
//             }
//         }

//     j = -1;
//     for (i = 1;i<know;i++)
//         if (ensemble[i].Q>inensem[max].Q) //not to disturb the ordering of elements of ensemble
//         {
//             j = i;    break;
//         }

//     if (j<0) //if inensem[max] is the largest one in ensemble
//     {
//         free(ensemble[0].pa);//to delete the worst partition in ensemble
//         for (i=1;i<know;i++)
//             ensemble[i-1]=ensemble[i];//ensemble 0 = ensemble 1
//         ensemble[know-1]=inensem[max];//put the inensem[max] as the last element in ensemble

//     }
//     else  //if inensem[max] is not the largest - ith element in ensemble is larger than that
//     {
//         free(ensemble[0].pa);
//         for (i = 1;i < j;i++)  //place inensem[max] after jth element in ensemble
//             ensemble[i - 1] = ensemble[i];
//         ensemble[j-1]=inensem[max];


//     }
// }


// void addone(struct part *ensemble, int know, struct part *inensem, int max, int *edgecc, int N)
// {         //to add one element(inensem[max]) in ensemble
//     long i, j;
//     for (i=0;i<N;i++)        cols[i]=0;

//     for (i = 0;i<N-1;i++)
//         if (!cols[i])
//         {
//             cols[i]=1;
//             for (j = i+1;j<N;j++)
//             {
//                 if (inensem[max].pa[i] == inensem[max].pa[j]) //f i and jth node are in same partition in       inensem[max],
//                 edgecc[i*N + j]++; // increase by 1
//                 if (edgecc[i*N+j]==know+1)  //when the nodes are present in all ensembles, jth element of cols is 1
//                 cols[j]=1;
//             }
//         }

//     j = -1;
//     for (i = 1;i<know;i++) //not to disturb the ordering of ensemble
//         if (ensemble[i].Q>inensem[max].Q)
//         {
//             j = i;    break;
//         }

//     if (j<0)  //if inensem[max] is the largest one in ensemble
//     {
//         ensemble[know]=inensem[max];   //put inensem[max] as the last element in  ensemble
//     }
//     else
//     {
//         for (i = know - 1;i >= j;i--)
//         {
//             ensemble[i + 1] = ensemble[i];  //o/w place the inensem[max] after jth element in ensemble
//         }
//         ensemble[j]=inensem[max];

//     }

// }
// clock_t TIME;
// void puttime(FILE *fo,char *mess)
// {
//     fprintf(fo,mess);
//     fprintf(fo,": %lf\n",(double)(clock()-TIME)/CLOCKS_PER_SEC);
//     printf(mess);
//     printf(": %lf\n",(double)(clock()-TIME)/CLOCKS_PER_SEC);
//     TIME=clock();
// }



int main(int argc, char *argv[])
{
   
    char command[256];
    long i,j,k;
    int size;
    size=1;
    unsigned int seed, *seedlist;
    int krg,know,copy1,copy2;
    int kmax,kp;
    FILE *ex;

    // ex = malloc(2*sizeof(FILE*));
    ex = fopen("ex1.dot","w");
    // ex[1] = fopen("ex2.dot","w");
    krg=atoi(argv[1]);
    copy1=atoi(argv[2]);
    kmax=copy1*size;
    seed=(unsigned int)time(NULL);

    seedlist=(unsigned int*)malloc(sizeof(unsigned int)*size);
    // for (i=0;i<size;i++)
    // {
    //     rand1(&seed);
    //     seedlist[i]=seed+i;
    //     //printf("%d", seedlist[i]);
    // }
    struct insert I;
    inputGraph(&I);
    // cout << G.nb_nodes<<endl;
    Graph g = convert_Insert_to_class_Graph(I);
    printf("HEllo");
    g.display_for_check();
    struct part ensemble;
    //struct part *ensemble; //ensemble array is for original graph, where iterensem is for reduced graph
    //ensemble=(struct part*)malloc(sizeof(struct part)*kmax);
    // printf("working on initialization\n");
    // int *edgecc, *score;
    // score=(int*)malloc(sizeof(int)*G.N);
    // edgecc=(int*)malloc(sizeof(int)*G.N*G.N);
    ensemble=Louvain(g);
    ofstream outfile;
    outfile.open("Com.txt", ofstream::trunc);
    outfile.close();

    // {
        

    //     for (j=0;j<copy1;j++)
    //     {

    //             ensemble[j*size]=RG(G,krg,seedlist);  //get the 1st ensemble - struct named part
    //             printf("Mod value %lf \n",ensemble[j].Q);
    //             // printf("Mod value d%d \n",ensemble[j].Q);
    //             // printf("Mod value %lf \n",ensemble[j].Q/(2*G.n));
    //             // printf("2.G.n %ld\n", 2*G.n);
    //             // printf("com value %d \n", ensemble[j].com);
    //     //get set of partitions
    //     }
    // }
    for (int i =1; i<I.N+1; i++)
    printf("%d partition %d \n", i, ensemble.pa[i]);
    fprintf(ex,"graph D1 {\n");
    fprintf(ex,"\tnode [size = fixed, shape = circle,  width = 0.20,  style= filled, colorscheme= paired12];\n\n");
    int numcom=0;
    for (int i=1; i<I.N+1; i++) {
        fprintf(ex,"\t%d [color=%d];\n",i,ensemble.pa[i]);
        if (ensemble.pa[i]>numcom) numcom=ensemble.pa[i];
    }
    fprintf(ex,"\n");
    // for (i=0; i<G.N; i++)
    //     for (j=1; j<=G.nl[i][0]; j++) if (i<G.nl[i][j])
    // {
    //     if (ensemble[1].pa[i]==ensemble[1].pa[G.nl[i][j]]) fprintf(ex[0],"edge [weight=10, len=1];\n");
    //     else fprintf(ex[0],"edge [weight=1, len=5];\n");
    //     fprintf(ex[0],"%d -- %d;\n",i+1,G.nl[i][j]+1);
    // }
    fprintf(ex,"}\n");
    fclose(ex);
    sprintf(command,"neato -Tpng -oex1.png ex1.dot > /dev/null 2>&1");
    system(command);
 //    know=cleansort(ensemble,kmax,G.N);
 //    cols=(int*)malloc(sizeof(int)*G.N);
 //    i=0;
 //    while (i<G.N) {cols[i]=0; i++;}

 //    for (i=0;i<G.N-1;i++)
 //       if (!cols[i])
 //    {
 //        cols[i]=1;
 //                for (j=i+1;j<G.N;j++)
 //        {
 //            edgecc[i*G.N+j]=0;
 //            for (k=0;k<know;k++)
 //                    if (ensemble[k].pa[i]==ensemble[k].pa[j])
 //                    edgecc[i*G.N+j]++;
 //                if (edgecc[i*G.N+j]==know)
 //            cols[j]=1;
 //        }
 //    }

 //    int iteration=0;
 //    int *ref, bp;
 //    // ref: new group label for each old group (group=community)
 //    // score: current group label for each original node (not reduced node)
 //    ref=(int*)malloc(sizeof(int)*G.N);
 //    for (i = 0;i < G.N;i++)
 //        score[i] = i + 1;

 //    int sizeofRN;
 //     sizeofRN=getscore(score,edgecc,G.N,know,ref);
 //        //update G, edgecc
 //    renorm(G,ref);
 //    G.com=sizeofRN;
 //    printf("second set\n");
 //    printf("G.com %d G.n %d\n", G.com, G.n);
 //     for (int i=0; i<G.com; i++)
 //    {
        
 //        printf("degreelist %d\n", G.dl[i]);   //print the degree for each node

 //        printf("neigh_list\n");
 //        for (int j =0; j<G.nl[i][0]+1; j++)
 //            {
 //            printf("%d %d\n",j, G.nl[i][j]); }  //print the neighbor list for each node
            
 //        printf("weight_list\n");
 //        for (int k =0; k<G.nl[i][0]+1; k++)
 //            printf("%d %d\n",k, G.wl[i][k]);   //print the weight list for each node
            
 //        printf("member_list \n");
 //        for (int k =0; k<G.ml[i][0]+1; k++)
 //            printf("%d %d\n",k, G.ml[i][k]);   //print the member list for each node
 //            printf("END\n");

            
 //    }
 //    struct part iterensem;
 //    iterensem=RG(G,krg,seedlist);
 //    printf("new_numcom %d\n", G.com);
 //    for (int i =0; i<G.com; i++)
 //    printf("partition %d \n", iterensem.pa[i]);
    // printf("hey\n");
 //    fprintf(ex[1],"graph D2 {\n");
 //    fprintf(ex[1],"\tnode [shape = circle,  width = 1.5,  style= filled, colorscheme= paired12];\n\n");
 //    int new_numcom=0;
    // for (int i=0; i<G.com; i++) {
    //     fprintf(ex[1],"\t%d [color=%d] [fontsize=50];\n",i+1,iterensem.pa[i]+1);
    //     if (iterensem.pa[i]>new_numcom) new_numcom=iterensem.pa[i];
    // }
    // fprintf(ex[1],"\n");
    // for (i=0; i<G.com; i++)
    //     if (G.wl[i][0]==0)
    // {
    //     for (j=1; j<=G.nl[i][0]; j++) if (i<G.nl[i][j])
    // {
    //     if (iterensem.pa[i]==iterensem.pa[G.nl[i][j]]) fprintf(ex[1],"edge [weight=10, len=1, penwidth=%d];\n", G.wl[i][j]);
    //     else fprintf(ex[1],"edge [weight=1, len=10, penwidth=%d];\n", G.wl[i][j]);

    //     fprintf(ex[1],"%d -- %d;\n",i+1,G.nl[i][j]+1);
    // }
    // }
    // else
    // {
    //     fprintf(ex[1],"edge [weight=10, len=50, penwidth=%d];\n", G.wl[i][0]);
    //     fprintf(ex[1],"%d:e -- %d:w;\n",i+1,i+1);
    //     for (j=1; j<=G.nl[i][0]; j++) if (i<G.nl[i][j])
    // {
    //     if (iterensem.pa[i]==iterensem.pa[G.nl[i][j]]) fprintf(ex[1],"edge [weight=10, len=5, penwidth=%d];\n", G.wl[i][j]);
    //     else fprintf(ex[1],"edge [weight=1, len=10, penwidth=%d];\n", G.wl[i][j]);

    //     fprintf(ex[1],"%d -- %d;\n",i+1,G.nl[i][j]+1);
    // }

    // }

    // //fprintf(ex[1],"edge [weight=10, penwidth = 10.];\n");
    // //fprintf(ex[1],"%d -- %d;\n",1,1);
 //    fprintf(ex[1],"}\n");
    // fclose(ex[1]);
    
    // sprintf(command,"neato -Tpng -oex2.png ex2.dot > /dev/null 2>&1");
    // system(command);
    return 0;
}

//read Graph data from file
void inputGraph(struct insert *I)
{

    
    int i,n1,n2;
    FILE *fi;
    fi=fopen("info.txt","r");
    fscanf(fi,"%d%d",&I->N,&I->n);
    fclose(fi);

    I->com=I->N;
    I->ml=(int**)malloc(sizeof(int*)*I->N);
    I->nl=(int**)malloc(sizeof(int*)*I->N);
    I->wl=(int**)malloc(sizeof(int*)*I->N);
    I->dl=(int*)malloc(sizeof(int)*I->N); // keeps the number of degrees for each node

    fi=fopen("degree.txt","r");
    for (i=0;i<I->N;i++)
        fscanf(fi,"%d",I->dl+i);  //from degree.txt, no of degrees has been printed in I->dl for each node
    fclose(fi);

    //fi=fopen("check.txt","r");
    for (i=0;i<I->N;i++)
    {
        I->ml[i]=(int*)malloc(sizeof(int)*2);
        I->ml[i][0]=1;
        I->ml[i][1]=i+1;
        I->nl[i]=(int*)malloc(sizeof(int)*(I->dl[i]+1));
        I->nl[i][0]=0;
        I->wl[i]=(int*)malloc(sizeof(int)*(I->dl[i]+1));
        I->wl[i][0]=0;
        //fprintf(fi,"%d",I->ml[i][0], I->ml[i][1]);
    }
    // for (int i =0; i<(sizeof(int)*2);i++)
    // {
    //     printf("%d\n",I->ml[4][i]);
    // }
    //fclose(fi);
    fi=fopen("clean.txt","r");
    for (i=0;i<I->n;i++)
    {
        fscanf(fi,"%d%d",&n1,&n2);
        n1--;    n2--;
        //printf("%d %d\n", n1,n2);

        I->nl[n1][0]++;    I->nl[n2][0]++;
        I->nl[n1][I->nl[n1][0]]=n2;
        I->wl[n1][I->nl[n1][0]]=1;
        I->nl[n2][I->nl[n2][0]]=n1;
        I->wl[n2][I->nl[n2][0]]=1;
       
    }
    //printf("I am %d %d %d %d\n", I->nl[1][0],I->nl[2][0],I->nl[3][0], I->nl[3][0]);
    //printf(sizeof(nl[1]));
    /*for (int i=0; i<I->dl[1]+1; i++)
    {
        printf("I am %d \n", I->nl[1][i]);
    }*/
    fclose(fi);

}


Graph convert_Insert_to_class_Graph(struct insert I)
{
   Graph g(I.N,I.n);
   cout<<g.nb_nodes<<endl;
   g.degrees[0]=0;
   //cout<<("check2")<<endl;
   {
    for (int i=0; i<g.nb_nodes; i++)
      {
         g.degrees[i+1] = cumulative_degree(i,I);
       }
   g.nb_links = g.degrees[g.nb_nodes-1];
   g.links.resize(g.nb_links);
   cout<< "nb_links"<<g.nb_links<<endl;
   cout<< "nb_nodes"<<g.nb_nodes<<endl;
   for (int i=0; i<g.nb_links; i++)
      {
        g.links[i] = Links(i,I);
        //cout<<g.links[i]<<endl;
        //cout<<"Links"<<Links(i,I)<<endl;
      }
    g.weights.resize(g.nb_links);
    //g.weights[0]=0;
    for (int i=0; i<g.nb_links; i++)
    {
      g.weights[i]=Weight(i,I);
      //cout<<"Weights"<<g.weights[i]<<endl;
    }
    // Compute total weight
    for (unsigned int i=0 ; i<g.nb_nodes ; i++)
        {
            g.total_weight += (double)g.weighted_degree(i);
        }
    g.denominator = I.n;
   }
   return g;
}

double Weight(unsigned int i,struct insert I)
    {
        vector <unsigned int> weights_list;
        for (int i=0; i<I.N; i++)
        {
            if (I.wl[i][0]==0)
                for (int j=0; j<I.nl[i][0]; j++)
                {
                    weights_list.push_back ((I.wl[i][j+1]));
                }
            else
                {//cout<<"check5"<<endl;
                    weights_list.push_back(I.wl[i][0]*2);
                    for (int j=0; j<I.nl[i][0]; j++)
                    {
                        weights_list.push_back ((I.wl[i][j+1]));
                    }
                }
        }
        
    return weights_list[i];

    }

double cumulative_degree(unsigned int i,struct insert I)
{
  vector <unsigned int> cum_deg(I.N);
   cum_deg[0]=I.nl[0][0];
  for (int i=1; i<I.N; i++)
   {
    if (I.wl[i][0]==0)
        cum_deg[i]=cum_deg[i-1]+I.nl[i][0];
    else
        cum_deg[i]=cum_deg[i-1]+I.nl[i][0]+1;
        //printf("cum_deg");}
    }
    //cout<<"Hi5"<<endl;
    return cum_deg[i];
 
}

unsigned int Links(int i, struct insert I)
{
   vector<unsigned int> getlink;
  for (int i = 0; i<I.N; i++)
    {
    if (I.wl[i][0]==0)
    {
        //cout<<"check4"<<endl;
        for (int j=0; j<I.nl[i][0]; j++)
        {
            getlink.push_back ((I.nl[i][j+1])+1);
        }
  
    }
    else
        {//cout<<"check5"<<endl;
        getlink.push_back(i+1);
        for (int j=0; j<I.nl[i][0]; j++)
            {
                getlink.push_back ((I.nl[i][j+1])+1);
            }
        }
}
return getlink[i];
}
