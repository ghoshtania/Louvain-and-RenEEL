/*
*
* File main.cpp
* This program implements the Reduced Network Extremal Ensemble Learning (RenEEL) scheme for maximizing modularity.
* Instructions for usage are contained in the file readme.txt
*
* If you use it for your work, please cite the following paper:
* Reduced network extremal ensemble learning (RenEEL) scheme for community detection in complex networks,
* For comments/questions and reporting any bugs please contact Kevin E. Bassler (bassler@uh.edu)
*
*/

#include<cstdio>
#include<cstdlib>
#include <climits>
#include <ctime>
//#include"help.h"
#include"graph.h"
#include "louvain.h"

#define Tsmall 0.0000001


int *cols;

int findxiny(int x, int *y);
double cumulative_degree(unsigned int i,struct insert I);
unsigned int Links(int i,struct insert I);
double Weight(unsigned int i,struct insert I);
Graph convert_Insert_to_class_Graph(struct insert I);
void inputGraph(struct insert *I);
int cleansort(struct part *ensem,int kmax,int N);
int getscore(int *score,int *edgecc, int N,int know,int *ref);
void renorm(struct insert I, int *ref);
int mergroup(struct insert *I, int x, int y);
int comps(int *p, int *q, int N);
void remolast(struct part *ensemble, int ens, int *edgecc, int N)
{
    long i;
    long j;
    for (i=0;i<N;i++)
        cols[i]=0;
    
    for (i = 0;i<N-1;i++)
        if (!cols[i])
        {
            cols[i]=1;
            for (j = i+1;j<N;j++)
            {
                if (ensemble[0].pa[i] == ensemble[0].pa[j])
                edgecc[i*N + j]--;
                if (edgecc[i*N+j]==ens-1)
                cols[j]=1;
            }
        
        }
    //delete in ensemble
    i = 0;
    free(ensemble[0].pa);
    for (i = 0;i<ens - 1;i++)
    {
        ensemble[i] = ensemble[i + 1];
    }
    ensemble[ens - 1].pa = NULL;
}

void replaceone(struct part* ensemble,int know, struct part* inensem, int max, int *edgecc, int N)
{
    long i,j;
    for (i=0;i<N;i++)    cols[i]=0;
    for (i = 0;i<N-1;i++)
        if (!cols[i])
        {
            cols[i]=1;
            for (j = i+1;j<N;j++)
            {
                if (inensem[max].pa[i] == inensem[max].pa[j])
                edgecc[i*N + j]++;
                if (ensemble[0].pa[i]==ensemble[0].pa[j])
                edgecc[i*N+j]--;
                if (edgecc[i*N+j]==know)
                cols[j]=1;
            }
        }

    j = -1;
    for (i = 1;i<know;i++)
        if (ensemble[i].Q>inensem[max].Q)
        {
            j = i;    break;
        }

    if (j<0)
    {
        free(ensemble[0].pa);
        for (i=1;i<know;i++)
            ensemble[i-1]=ensemble[i];
        ensemble[know-1]=inensem[max];

    }
    else
    {
        free(ensemble[0].pa);
        for (i = 1;i < j;i++)
            ensemble[i - 1] = ensemble[i];
        ensemble[j-1]=inensem[max];


    }
}


void addone(struct part *ensemble, int know, struct part *inensem, int max, int *edgecc, int N)
{
    long i, j;
    for (i=0;i<N;i++)        cols[i]=0;

    for (i = 0;i<N-1;i++)
        if (!cols[i])
        {
            cols[i]=1;
            for (j = i+1;j<N;j++)
            {
                if (inensem[max].pa[i] == inensem[max].pa[j])
                edgecc[i*N + j]++;
                if (edgecc[i*N+j]==know+1)
                cols[j]=1;
            }
        }

    j = -1;
    for (i = 1;i<know;i++)
        if (ensemble[i].Q>inensem[max].Q)
        {
            j = i;    break;
        }

    if (j<0)
    {
        ensemble[know]=inensem[max];
    }
    else
    {
        for (i = know - 1;i >= j;i--)
        {
            ensemble[i + 1] = ensemble[i];
        }
        ensemble[j]=inensem[max];

    }

}
clock_t TIME;
void puttime(FILE *fo,char *mess)
{
    fprintf(fo,mess);
    fprintf(fo,": %lf\n",(double)(clock()-TIME)/CLOCKS_PER_SEC);
    printf(mess);
    printf(": %lf\n",(double)(clock()-TIME)/CLOCKS_PER_SEC);
    TIME=clock();
}

double rand1(unsigned int* seed)
{   //printf(RAND_MAX);
    return (rand_r(seed)/(RAND_MAX*1.0+1));

}

double comp(struct insert I)
{
    double ans=0;
    int i;
    for (i=0;i<I.com;i++)
    {

        ans+=2*I.wl[i][0]-I.dl[i]*1.0*I.dl[i]/(2*I.n);
    }
    return (ans);
}

int main(int argc, char *argv[])
{


    long i,j,k;
    int size;
    // FILE *fout,*frc;
    // fout=fopen("result.txt","w");
    // frc=fopen("records.txt","a");

    size=1;


    unsigned int seed, *seedlist;
    int know,copy1,copy2;
    int kmax,kp;

    //krg=atoi(argv[1]);
    copy1=atoi(argv[1]);
    kmax=copy1*size;
    copy2=atoi(argv[2]);
    kp=copy2*size;

    srand(time(NULL));
    //fprintf(fout,"seed=%d\n",seed);
    //fprintf(frc,"seed=%d ",seed);
    //fprintf(fout,"kmax=%d\nkp=%d\n",kmax,kp);
    //fprintf(frc,"kmax=%d kp=%d ",kmax,kp);
    printf("seed=%d\n",seed);
    printf("kmax=%d\nkp=%d\n",kmax,kp);
    printf("size",size);
    seedlist=(unsigned int*)malloc(sizeof(unsigned int)*size);
    for (i=0;i<size;i++)
    {
        rand1(&seed);
        seedlist[i]=seed+i;
    }
    // cout << "hello"<<endl;
    struct insert I;
    inputGraph(&I);
    // cout << G.nb_nodes<<endl;
    Graph g = convert_Insert_to_class_Graph(I); 
    //display_for_check();
    cout<<"here it is"<<endl;
    printf("1st graph");
    for (int i=0; i<I.com; i++)
    {
        printf("degreelist %d\n", I.dl[i]);   //print the degree for each node

        printf("neigh_list\n");
        for (int j =0; j<I.nl[i][0]+1; j++) 
            printf("%d %d\n",j, I.nl[i][j]+1);   //print the neighbor list for each node    
        // printf("weight_list\n");
        // for (int k =0; k<I.nl[i][0]+1; k++)
        //     printf("%d %d\n",k, I.wl[i][k]);   //print the weight list for each node
            
        // printf("member_list \n");
        // for (int k =0; k<I.ml[i][0]+1; k++)
        //     printf("%d %d\n",k, I.ml[i][k]);   //print the member list for each node
            printf("END\n");
    }
    //g.display_for_check();
    //Louvain(g);
 //    fprintf(fout,"links=%d nodes=%d\n",nb_links,nb_nodes);
 //    printf("links=%d nodes=%d\n",nb_links,nb_nodes);
 //    fprintf(frc,"links=%d nodes=%d ",nb_links,nb_nodes);
 //       puttime(fout,"reading time");
 // //input network

    struct part *ensemble,*iterensem;
    ensemble=(struct part*)malloc(sizeof(struct part)*kmax);
    iterensem=(struct part*)malloc(sizeof(struct part)*kp);

 printf("working on initialization\n");
 // fprintf(fout,"working on initialization\n");
   
    int *edgecc, *score;
    score=(int*)malloc(sizeof(int)*I.N);
    edgecc=(int*)malloc(sizeof(int)*I.N*I.N);

    
    {    for (j=0;j<copy1;j++)
        {

                ensemble[j*size]=Louvain(g);
                ofstream outfile;
                outfile.open("Com.txt", ofstream::trunc);
                outfile.close();
                //cout<<"Mod value"<<ensemble[j].Q;
            
        //get set of partitions
        }
    }
    // for (j=0;j<copy1;j++)
    // {
    //     cout<<"Mod Value"<<ensemble[j].Q<<endl;
    //     for (int i =0; i<35; i++)
    //         cout<<i<<" "<<ensemble[j].pa[i]<<endl;
    //     cout<<endl;
    // }

 //    //puttime(fout,"initialization time");
     know=cleansort(ensemble,kmax,I.N);
     cerr<<"cleansort"<<know<<endl;
    // puttime(fout,"sort time");
     //cout<<"know value\n"<<know;
    cols=(int*)malloc(sizeof(int)*I.N);
    i=0;
    while (i<I.N) {cols[i]=0; i++;}

    for (i=0;i<I.N-1;i++)
       if (!cols[i])
    {
        cols[i]=1;
                for (j=i+1;j<I.N;j++)
        {
            edgecc[i*I.N+j]=0;
            for (k=0;k<know;k++)
                    if (ensemble[k].pa[i]==ensemble[k].pa[j])
                    edgecc[i*I.N+j]++;
                if (edgecc[i*I.N+j]==know)
            cols[j]=1;
        }
    }
   

 // // //    //initialization
 // // // //    puttime(fout,"count edgecc time");

    int iteration=0;
    int *ref, bp;
    ref=(int*)malloc(sizeof(int)*I.N);
    for (i = 0;i < I.N;i++)
        score[i] = i + 1;


    cerr<<"2nd graph"<<endl;
      int sizeofRN;
        sizeofRN=getscore(score,edgecc,I.N,know,ref);
        //update I, edgecc

    renorm(I,ref);
        I.com=sizeofRN;
     
 
    printf("I.com %d I.N %d\n", I.com,I.N); 
    //FILE *fout1;
    //fout1=fopen("second_graph.txt","w"); //to write
    for (int i=0; i<I.com; i++)
    {
        printf("degreelist %d\n", I.dl[i]);   //print the degree for each node

        printf("neigh_list\n");
        for (int j =0; j<I.nl[i][0]+1; j++) 
            printf("%d %d\n",j, I.nl[i][j]+1);   //print the neighbor list for each node
            
        printf("weight_list\n");
        for (int k =0; k<I.nl[i][0]+1; k++)
            printf("%d %d\n",k, I.wl[i][k]);   //print the weight list for each node
            
        printf("member_list \n");
        for (int k =0; k<I.ml[i][0]+1; k++)
            printf("%d %d\n",k, I.ml[i][k]);   //print the member list for each node
            printf("END\n");

        // if (I.wl[i][0] == 0)
        //     {for (int j =1; j<I.nl[i][0]+1; j++) 
        //     fprintf(fout1,"%d %d %.2f\n",i+1, I.nl[i][j]+1, float(I.wl[i][j])/2.);  } //print the neighbor list for each node 
        // else
        //    { fprintf(fout1,"%d %d %d\n",i+1, i+1, I.wl[i][0]*2);
        //      for (int j =1; j<I.nl[i][0]+1; j++) 
        //         fprintf(fout1,"%d %d %.2f\n",i+1, I.nl[i][j]+1, float(I.wl[i][j])/2.);   //print the neighbor list for each node 
        //     } 

    }
   //fclose(fout1);
   Graph g2 = convert_Insert_to_class_Graph(I);
   //display_for_check();
    cout<<"here it is"<<endl;
        //partitions of newG
        for (j=0;j<copy2;j++)
    {

        {
            
                    iterensem[j*size]=Louvain(g2);
                    ofstream outfile;
                    outfile.open("Com.txt", ofstream::trunc);
                    outfile.close();

        }
           }
    for (j=0;j<copy1;j++) cout<<"Mod value 1st ensemble"<<ensemble[j].Q<<endl;
    for (j=0;j<copy2;j++) cout<<"Mod value 2nd ensemble"<<iterensem[j].Q<<endl; 


    //puttime(fout,"RG time");
        //pick the best
        bp = 0;    j = 1;
        for (i=1;i<kp;i++)
            if (abs(iterensem[i].Q - iterensem[bp].Q) < Tsmall)
            {
                j++;
                if (rand1(seedlist)*j < 1)
                {
                    free(iterensem[bp].pa);
                    bp = i;
                }
                else
                    free(iterensem[i].pa);
            }
            else if (iterensem[i].Q > iterensem[bp].Q)
            {
                j = 1;
                free(iterensem[bp].pa);
                bp=i;
            }
            else
            {
                free(iterensem[i].pa);
            }
    cout<<"bestpartition"<<iterensem[bp].Q<<endl;
        //update ensemble and edgecc
        i = -1;
        for (j = 0;j<know;j++)
            if (abs(iterensem[bp].Q -ensemble[j].Q )<Tsmall)
            {   cout<<"iterensem[bp].com"<<iterensem[bp].com<<endl;
                cout<<"ensemble[j].com"<<ensemble[j].com<<endl;
                if (iterensem[bp].com == ensemble[j].com)
                { for (int i=0; i<I.N; i++) cout<<"ensemble[j]"<<ensemble[j].pa[i]<<endl;
                    cout<<endl;
                  for (int i=0; i<I.N; i++) cout<<"iterensem[bp]"<<iterensem[bp].pa[i]<<endl;
                    if (comps(iterensem[bp].pa, ensemble[j].pa, I.N))
                    {   cout<<"yes"<<endl;
                        i = j;    break;
                    }
                }
            }

        if (i >= 0)
        {
            cout<<"yes1"<<endl;
            remolast(ensemble, know, edgecc, I.N);
            know--;
        }
        else
        {
            if (iterensem[bp].Q<ensemble[0].Q)
            {
                cout<<"yes2"<<endl;
                remolast(ensemble, know, edgecc, I.N); know--;
            }
            else
            {
                if (know<kmax)
                {

                    cout<<"yes3"<<endl;
                    addone(ensemble, know, iterensem, bp, edgecc, I.N);
                    know++;
                }
                else
                {

                    cout<<"yes4"<<endl;
                    cout<<"know and kmax"<<know<<endl<<kmax<<endl;
                    cout<<"Mod of bp"<<iterensem[bp].Q<<endl<<"worst ensemble mod"<<ensemble[0].Q<<endl;
                    replaceone(ensemble, know, iterensem, bp, edgecc, I.N);
                    
                }
            }
        }
   
    //T1=clock();
 //    while (know>1)
 //    { 
 //        iteration++;

 //        printf("%d:size=%d,ensem=%d,%lf,%lf\n",iteration,I.com,know,ensemble[0].Q,ensemble[know-1].Q);
 //       // fprintf(fout,"%d:size=%d,ensem=%d,%lf,%lf\n",iteration,G.com,know,ensemble[0].Q/(2*nb_links),ensemble[know-1].Q/(2*nb_links));

 //      int sizeofRN;
 //        sizeofRN=getscore(score,edgecc,I.N,know,ref);
 //        //update I, edgecc
 //        cout<<"Size of RN "<<sizeofRN<<endl;
 //    renorm(I,ref);
 //        //I.N=sizeofRN;
 //    printf("I.com %d I.N %d\n", I.com,I.N); 
 //        I.com=sizeofRN;
 //        printf("I.com %d I.N %d\n", I.com,I.N); 
 //    printf("2nd set\n"); 
 //     FILE *fout1;
 //    fout1=fopen("second_graph.txt","w"); //to write
 //    for (int i=0; i<I.com; i++)
 //    {
 //        printf("degreelist %d\n", I.dl[i]);   //print the degree for each node

 //        printf("neigh_list\n");
 //        for (int j =0; j<I.nl[i][0]+1; j++) 
 //            printf("%d %d\n",j, I.nl[i][j]+1);   //print the neighbor list for each node
            
 //        printf("weight_list\n");
 //        for (int k =0; k<I.nl[i][0]+1; k++)
 //            printf("%d %d\n",k, I.wl[i][k]);   //print the weight list for each node
            
 //        printf("member_list \n");
 //        for (int k =0; k<I.ml[i][0]+1; k++)
 //            printf("%d %d\n",k, I.ml[i][k]);   //print the member list for each node
 //            printf("END\n");
 //        if (I.wl[i][0] == 0)
 //            {for (int j =1; j<I.nl[i][0]+1; j++) 
 //            fprintf(fout1,"%d %d %.2f\n",i+1, I.nl[i][j]+1, float(I.wl[i][j])/2.);  } //print the neighbor list for each node 
 //        else
 //           { fprintf(fout1,"%d %d %d\n",i+1, i+1, I.wl[i][0]*2);
 //             for (int j =1; j<I.nl[i][0]+1; j++) 
 //                fprintf(fout1,"%d %d %.2f\n",i+1, I.nl[i][j]+1, float(I.wl[i][j])/2.);   //print the neighbor list for each node 
 //            } 

 //    }
 //   fclose(fout1);
 
 //    printf("I.com %d I.N %d\n", I.com,I.N); 
 //    Graph g2 = convert_Insert_to_class_Graph(I);
 //    //g2.display_for_check();



 // //    //puttime(fout,"reducing time");
    

 //        //partitions of newG
 //        for (j=0;j<copy2;j++)
 //    {

 //        {
            
 //                    iterensem[j*size]=Louvain(g2);
 //                    ofstream outfile;
 //                    outfile.open("Com.txt", ofstream::trunc);
 //                    outfile.close();

 //        }
 //           }
 //    for (j=0;j<copy1;j++) cout<<"Mod value 1st ensemble"<<ensemble[j].Q<<endl;
 //    for (j=0;j<copy2;j++) cout<<"Mod value 2nd ensemble"<<iterensem[j].Q<<endl; 


 //    //puttime(fout,"RG time");
 //        //pick the best
 //        bp = 0;    j = 1;
 //        for (i=1;i<kp;i++)
 //            if (abs(iterensem[i].Q - iterensem[bp].Q) < Tsmall)
 //            {
 //                j++;
 //                if (rand1(seedlist)*j < 1)
 //                {
 //                    free(iterensem[bp].pa);
 //                    bp = i;
 //                }
 //                else
 //                    free(iterensem[i].pa);
 //            }
 //            else if (iterensem[i].Q > iterensem[bp].Q)
 //            {
 //                j = 1;
 //                free(iterensem[bp].pa);
 //                bp=i;
 //            }
 //            else
 //            {
 //                free(iterensem[i].pa);
 //            }
 //    cout<<"bestpartition"<<iterensem[bp].Q<<endl;
 //        //update ensemble and edgecc
 //        i = -1;
 //        for (j = 0;j<know;j++)
 //            if (abs(iterensem[bp].Q -ensemble[j].Q )<Tsmall)
 //            {   cout<<"iterensem[bp].com"<<iterensem[bp].com<<endl;
 //                cout<<"ensemble[j].com"<<ensemble[j].com<<endl;
 //                if (iterensem[bp].com == ensemble[j].com)
 //                { for (int i=0; i<I.N; i++) cout<<"ensemble[j]"<<ensemble[j].pa[i]<<endl;
 //                    cout<<endl;
 //                  for (int i=0; i<I.N; i++) cout<<"iterensem[bp]"<<iterensem[bp].pa[i]<<endl;
 //                    if (comps(iterensem[bp].pa, ensemble[j].pa, I.N))
 //                    {   cout<<"yes"<<endl;
 //                        i = j;    break;
 //                    }
 //                }
 //            }

 //        if (i >= 0)
 //        {
 //            cout<<"yes1"<<endl;
 //            remolast(ensemble, know, edgecc, I.N);
 //            know--;
 //        }
 //        else
 //        {
 //            if (iterensem[bp].Q<ensemble[0].Q)
 //            {
 //                cout<<"yes2"<<endl;
 //                remolast(ensemble, know, edgecc, I.N); know--;
 //            }
 //            else
 //            {
 //                if (know<kmax)
 //                {

 //                    cout<<"yes3"<<endl;
 //                    addone(ensemble, know, iterensem, bp, edgecc, I.N);
 //                    know++;
 //                }
 //                else
 //                {

 //                    cout<<"yes4"<<endl;
 //                    cout<<"know and kmax"<<know<<endl<<kmax<<endl;
 //                    cout<<"Mod of bp"<<iterensem[bp].Q<<endl<<"worst ensemble mod"<<ensemble[0].Q<<endl;
 //                    replaceone(ensemble, know, iterensem, bp, edgecc, I.N);
                    
 //                }
 //            }
 //        }

 //    //puttime(fout,"update ensemble time");
 //    cout<<"I.N and I.com"<<endl<<I.N<<endl<<I.com<<endl;


 //        //free iterensem
 //   }
 
 //    // fprintf(fout,"rest time=%lf\n",(double)(clock()-T1)/CLOCKS_PER_SEC);
 //    // fprintf(fout,"real walltime=%d\n",time(NULL)-sec);
 //    // fprintf(frc,"real time=%d ",time(NULL)-sec);
 //    // printf("rest time=%lf\n",(double)(clock()-T1)/CLOCKS_PER_SEC);
 //    // printf("real walltime=%d\n",time(NULL)-sec);

 //    //get score of reduced network
 //       int sizeofRN;
 //        sizeofRN=getscore(score,edgecc,I.N,know,ref);
 //        //update G, edgecc
 //        renorm(I,ref);
 //        I.com=sizeofRN;
 // //   for (i=0;i<know;i++)
 //    //outpart(ensemble[0],nb_links, nb_nodes);
 //    double Qfinal;
 //    Qfinal=comp(I);                                ////need change.....
 //    printf("Qfinal=%lf\n",Qfinal/(2*I.N));
    //fprintf(fout,"Qfinal=%lf\n",Qfinal/(2*I.N));
    //fprintf(frc,"Qf=%lf\n",Qfinal/(2*I.N));

    //fclose(frc);
    //fclose(fout);
    return 0;
 }


int mergroup(struct insert *I, int x, int y)
{
    int i;
    int *newlist,*newl2;
    
    
    newlist=(int*)malloc(sizeof(int)*(I[0].ml[x][0]+I[0].ml[y][0]+1));
  
    for (i=1;i<=I[0].ml[x][0];i++)    newlist[i]=I[0].ml[x][i];
    for (i=1;i<=I[0].ml[y][0];i++)    newlist[i+I[0].ml[x][0]]=I[0].ml[y][i];

    newlist[0]=I[0].ml[x][0]+I[0].ml[y][0];

    free(I[0].ml[x]);
    free(I[0].ml[y]);

    I->ml[x]=newlist;
    I->ml[y]=NULL;
    
    //printf("update1\n");
    newlist=(int*)malloc(sizeof(int)*(I[0].nl[x][0]+I[0].nl[y][0]+1));
    newl2  =(int*)malloc(sizeof(int)*(I[0].nl[x][0]+I[0].nl[y][0]+1));
    newl2[0]=I->wl[x][0]+I->wl[y][0];
    newlist[0]=0;
    for (i=1;i<=I[0].nl[x][0];i++)
    {
        if (I[0].nl[x][i]==y){    newl2[0]+=I[0].wl[x][i];    }
        else{
            newlist[0]++;
            newlist[newlist[0]]=I[0].nl[x][i];
            newl2[newlist[0]]=I[0].wl[x][i];
        }
    }

    for (i=1;i<=I[0].nl[y][0];i++)
        if (I[0].nl[y][i]!=x)
        {
            int sg,k1,k2;
            sg=I[0].nl[y][i];
            k1=findxiny(x,I[0].nl[sg]);
            k2=findxiny(y,I[0].nl[sg]);
            if (k1>0)
            {
                I[0].wl[sg][k1]+=I[0].wl[sg][k2];
                int test;
                test=findxiny(sg,newlist);
                if (test>0)
                    newl2[test]=I[0].wl[sg][k1];
                else
                {
                    printf("Error2!\n");
                    exit(0);
                }
                
                if (k2!=I[0].nl[sg][0])
                {
                    I[0].nl[sg][k2]=I[0].nl[sg][I[0].nl[sg][0]];
                    I[0].wl[sg][k2]=I[0].wl[sg][I[0].nl[sg][0]];
                }
                I[0].nl[sg][0]--;
                
            }
            else
            {
                if (k2>0)
                    I[0].nl[sg][k2]=x;
                else
                {
                    printf("Error3!\n");    exit(0);
                }
                newlist[0]++;
                newlist[newlist[0]]=sg;
                newl2[newlist[0]]=I[0].wl[sg][k2];
            }
        }
    free(I[0].nl[x]);    free(I[0].wl[x]);
    I->nl[x]=newlist;
    I->wl[x]=newl2;
    free(I[0].nl[y]);    free(I[0].wl[y]);
    I->nl[y]=NULL;
    I->wl[y]=NULL;

    
    I[0].dl[x]+=I[0].dl[y];

    return 0;
}

void movegroup(struct insert *I, int f, int t)
{
    int i,j;
    I->ml[t]=I->ml[f];
    I->ml[f] = NULL;
   
    for (i=1;i<=I->nl[f][0];i++)
    {
        j=findxiny(f,I->nl[I->nl[f][i]]);
        I->nl[I->nl[f][i]][j]=t;
    }
    
    I->nl[t]=I->nl[f];
    I->wl[t]=I->wl[f];
    I->dl[t]=I->dl[f];
    I->nl[f] = NULL;
    I->wl[f] = NULL;
}


void renorm(struct insert I,int *ref)
{ 
    int i,j;
    //renormlize G
    for (i=0;i<I.com;i++)
    {  printf("I.dl in renorm subroutine %d\n",I.dl[i]);
       printf("ref[i]  %d\n",ref[i]);
        if (ref[i]!=i)
        {

            if (I.ml[ref[i]] == NULL)
                movegroup(&I,i,ref[i]);
            else
            mergroup(&I,ref[i],i);
        }
    }
    
}


int comps(int *p, int *q, int N)  //compare-they are same, return 1, ow return 0
{
    int i;
    for (i = 0;i < N;i++)
        if (p[i] != q[i])
            return 0;
    
    return 1;
}

int getscore(int *score,int *edgecc, int N,int know,int *ref)
{
    cout<<"N "<<N<<endl;
    long i,j,k;
    int groupsnumber;
    int *newsc;
    newsc = (int*)malloc(sizeof(int)*N);
    for (i=0;i<N;i++)
        newsc[i]=0;
    groupsnumber=0;

    for (i=0;i<N;i++)
    {
        if (!newsc[i])
        {
            ref[score[i]-1]=groupsnumber;
            groupsnumber++;
            newsc[i]=groupsnumber;
        
//below may include in above section
            for (j=i+1;j<N;j++)
            {
 
                if (edgecc[i*N+j]==know)
                {
            ref[score[j] - 1] = newsc[i]-1;
                newsc[j]=newsc[i];
                
                }
            }
    }
    }
    for (i = 0;i < N;i++)
        {score[i] = newsc[i];
            cout<<"score[i] "<<score[i]<<endl;}

    free(newsc);
    return (groupsnumber);
}

int cleansort(struct part *ensem, int kmax,int N)
{   //cout<<"HI"<<endl;
    int i,j;
    struct part temp;
    for (i=0; i<kmax; i++)
        cout<<ensem[i].Q<<endl;
    for (i=0;i<kmax-1;i++)
        for (j=i+1;j<kmax;j++)
        {
            if (ensem[i].Q>ensem[j].Q)
            {
                temp=ensem[i];  ensem[i]=ensem[j]; ensem[j]=temp;
            }
        }

    int *del,len;
    del=(int*)malloc(sizeof(int)*kmax);
    for (i=0;i<kmax;i++)
        del[i]=0;
    for (i=0;i<kmax-1;i++)
        if (del[i]==0)
        {
            j=i+1;
            //cout<<"mod value of i and j"<<ensem[i].Q<<" "<<ensem[j].Q;
            while (j<kmax&&(ensem[j].Q-ensem[i].Q)<0.1)
            {
                if (ensem[j].com==ensem[i].com)
                    if (comps(ensem[i].pa,ensem[j].pa,N))
                        del[j]=1;
                j++;
                //cout<<"true"<<endl;
            }
        }
        //cout<<"value of j"<<j;
    for (i=0;i<kmax-1;i++)
        for (j=i+1;j<kmax;j++)
        {
            if (del[i]>del[j])
            {
                del[i]=0;   del[j]=1;
                temp=ensem[i];  ensem[i]=ensem[j]; ensem[j]=temp;
            }
            else if (del[i]+del[j]==0)
            {
                if (ensem[i].Q>ensem[j].Q)
                {
                    temp=ensem[i];  ensem[i]=ensem[j]; ensem[j]=temp;
                }
            }
        }
    len=0;
    while (len<kmax&&del[len]==0) len++;
    i=len;
    while (i<kmax)
    {
        free(ensem[i].pa);
        i++;
    }
    free(del);
    //cout<<"check"<<endl;
    
    return len;
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
   Graph g(I.com,I.n);
   g.nb_com = I.com;
   cout<<"nodes here"<<g.nb_nodes<<endl;
   cout<<"Com here"<<g.nb_com<<endl;
   g.degrees.resize(g.nb_com+1);
   g.degrees[0]=0;
   
   {
    for (int i=0; i<g.nb_com; i++)
      {  
         g.degrees[i+1] = cumulative_degree(i,I);
       }
    
   g.nb_links = g.degrees[g.nb_com];
   cout<<"check"<<endl;
   g.links.resize(g.nb_links);

   int row;
   int col;
   row = g.nb_com;
   g.member.resize(row+1);
   
   // g.member.resize(row+1, vector<int> (col+1));
   cout<<"check2"<<endl; ////?????
   cout<<"row "<<row<<endl;
   
   for (int i=1; i<=row; i++)
   {
     col = I.ml[i-1][0];
     cout << "col"<<col<<endl;
     g.member[i].resize(col+1);
     for (int j=1; j<=col; j++){
        g.member[i][j]=I.ml[i-1][j];
     }
   }
   cout<<"check1"<<endl; ////?????
    for (int i=1; i<=row; i++)
    {

      for (int j=1; j<g.member[i].size(); j++) 
        {
            cout<<"member"<<g.member[i][j]<<" ";
        }
        cout<<endl;
    }
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

double cumulative_degree(unsigned int i,struct insert I)
{  
  vector <unsigned int> cum_deg(i+1);
   if (I.wl[0][0] == 0) cum_deg[0]=I.nl[0][0];
   else cum_deg[0]=I.nl[0][0]+1;
  for (int m=1; m<i+1; m++) 
    {
    if (I.wl[m][0] == 0)
        cum_deg[m]=cum_deg[m-1]+I.nl[m][0];

    else
       cum_deg[m]=cum_deg[m-1]+I.nl[m][0]+1;
    }
     //cout<<"cum_deg"<<cum_deg[i]<<' '<<i<<endl;
  
  return cum_deg[i];
 
}

double Weight(unsigned int i,struct insert I)
    {
        vector <unsigned int> weights_list;
        for (int i=0; i<I.com; i++)
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
    //cout<<"weights_list"<<weights_list[i]<<endl; 
    return weights_list[i];

    }
unsigned int Links(int i, struct insert I)
{ 
   vector<unsigned int> getlink;
  for (int i = 0; i<I.com; i++)
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
//cout<<"links_list"<<getlink[i]<<endl;        
return getlink[i];
}

int findxiny(int x, int *y)
{
    int j;
    for (j=1;j<=y[0];j++)
    
    if (y[j]==x) return (j);
    return (-1);
}