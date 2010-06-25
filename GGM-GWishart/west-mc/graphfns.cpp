#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <iomanip>
#include <time.h>
#include "graph.h"
#include "various.h"

int rand_int(int low, int high){
double cutoff;
int i;

cutoff=drand48();
i=low-1;

while(cutoff> ((double) i-(low-1))/(double (high-(low-1)))) i++;

return i;
}

void Add(LPGraph graph, int*  edge){
  graph->Edge[edge[0]][edge[1]]=1;
  graph->Edge[edge[1]][edge[0]]=1;
}

void remove(LPGraph graph, int*  edge){
  graph->Edge[edge[0]][edge[1]]=0;
  graph->Edge[edge[1]][edge[0]]=0;
}


int * rand_edge_add(LPGraph graph, int current_edges,  int total_edges){
  int *edge;
  int edge_index,edge_counter,i,j;
  extern mwSize NumberOfGenes;
  edge=new int[2];
  edge_index=rand_int(1, total_edges-current_edges);
  //cout<<"edge index "<<edge_index<<endl;
  edge_counter=0;
  i=0;
  j=1;
  while(edge_counter < edge_index){
    if(graph->Edge[i][j]==0){
      edge_counter++;
      //cout<<i<<"  "<<j<<"  "<<edge_counter<<endl;
    }
    if(edge_counter < edge_index && j== NumberOfGenes-1){
      i=i+1;
      j=i+1;
    }
    else if(edge_counter < edge_index) j++;
  }

  edge[0]=i;
  edge[1]=j;
  Add(graph, edge);

  return edge;
}


int * rand_edge_delete(LPGraph graph, int current_edges){
  int *edge;
  int edge_index,edge_counter,i,j;
  extern mwSize NumberOfGenes;
  edge=new int[2];
  edge_index=rand_int(1, current_edges);

  edge_counter=0;
  i=0;
  j=1;
  while(edge_counter < edge_index){
    if(graph->Edge[i][j]==1) edge_counter++;
    if(edge_counter < edge_index && j== NumberOfGenes-1){
      i=i+1;
      j=i+1;
    }
    else if(edge_counter < edge_index) j++;
  }

  edge[0]=i;
  edge[1]=j;
  remove(graph, edge);

  return edge;
}
