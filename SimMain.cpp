#include<cstdio>
#include<cstdlib>
#include <climits>
#include <ctime>
#include"SimCommunity.h"

void inputGraph(struct graph *G);

int main(int argc, char *argv[]){

    int krg,know,copy1,copy2;
    int kmax,kp;
    long i,j,k;
    int size;
    FILE *fout,*frc, *fout1;
    fout=fopen("result.txt","w");
    frc=fopen("records.txt","a");

    size=1;
    krg=atoi(argv[1]);
    copy1=atoi(argv[2]);
    kmax=copy1*size;
    copy2=atoi(argv[3]);
    kp=copy2*size;
    fprintf(fout,"krg=%d\nkmax=%d\nkp=%d\n",krg,kmax,kp);
    fprintf(frc,"krg=%d kmax=%d kp=%d ",krg,kmax,kp);
    printf("krg=%d\nkmax=%d\nkp=%d\n",krg,kmax,kp);
    struct graph G;
    inputGraph(&G); //insert the graph information in dl, ml, nl,wl
    Community(struct G);
    fclose(frc);
    fclose(fout);
    return 0;
}


void inputGraph(struct graph *G)
{

    
    int i,n1,n2;
    FILE *fi;
    fi=fopen("info.txt","r");
    fscanf(fi,"%d%d",&G->N,&G->n);
    fclose(fi);

    G->com=G->N;
    G->ml=(int**)malloc(sizeof(int*)*G->N);
    G->nl=(int**)malloc(sizeof(int*)*G->N);
    G->wl=(int**)malloc(sizeof(int*)*G->N);
    G->dl=(int*)malloc(sizeof(int)*G->N); // keeps the number of degrees for each node

    fi=fopen("degree.txt","r");
    for (i=0;i<G->N;i++)
        fscanf(fi,"%d",G->dl+i);  //from degree.txt, no of degrees has been printed in G->dl for each node
    fclose(fi);

    //fi=fopen("check.txt","r");
    for (i=0;i<G->N;i++)
    {
        G->ml[i]=(int*)malloc(sizeof(int)*2);
        G->ml[i][0]=1;
        G->ml[i][1]=i+1;
        G->nl[i]=(int*)malloc(sizeof(int)*(G->dl[i]+1));
        G->nl[i][0]=0;
        G->wl[i]=(int*)malloc(sizeof(int)*(G->dl[i]+1));
        G->wl[i][0]=0;
        //fprintf(fi,"%d",G->ml[i][0], G->ml[i][1]);
    }
    // for (int i =0; i<(sizeof(int)*2);i++)
    // {
    //     printf("%d\n",G->ml[4][i]);
    // }
    //fclose(fi);
    fi=fopen("clean.txt","r");
    for (i=0;i<G->n;i++)
    {
        fscanf(fi,"%d%d",&n1,&n2);
        n1--;    n2--;
        //printf("%d %d\n", n1,n2);

        G->nl[n1][0]++;    G->nl[n2][0]++;
        G->nl[n1][G->nl[n1][0]]=n2;
        G->wl[n1][G->nl[n1][0]]=1;
        G->nl[n2][G->nl[n2][0]]=n1;
        G->wl[n2][G->nl[n2][0]]=1;
       
    }
    //printf("I am %d %d %d %d\n", G->nl[1][0],G->nl[2][0],G->nl[3][0], G->nl[3][0]);
    //printf(sizeof(nl[1]));
    
    fclose(fi);

}
