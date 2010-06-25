/*
 g++ newgraph.cpp -o newgraph.exe
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <iostream>
#include <iomanip>
#include "matrix.h"
#include "graph.h"

#define max(a,b) ((a>b)?a:b)
#define MAXDIMENS 50

///////////////////////////////////////////////
extern void CheckPointer(void* pointer);
extern int  choose(int m, int n);
///////////////////////////////////////////////

//for qsort
int numeric (const void *p1, const void *p2)
{
   return(*((int*)p1) - *((int*)p2));
}

//class Graph::Begins
Graph::Graph()
{
   nVertices                 = 0;
   Edge                      = NULL;
   Labels                    = NULL;
   nLabels                   = 0;
   Cliques                   = NULL;
   CliquesDimens             = NULL;
   nCliques                  = 0;
   StarComp                  = NULL;
   ConnectedComponents       = NULL;
   ConnectedComponentsDimens = NULL;
   nConnectedComponents      = 0;
   TreeEdgeA                 = NULL;
   TreeEdgeB                 = NULL;
   nTreeEdges                = 0;
   nMss                      = 0;
   Mss                       = NULL;
   MssDimens                 = NULL;
   ordering                  = NULL;
   Separators                = NULL;
   SeparatorsDimens          = NULL;
   nSeparators               = 0;
   localord                  = NULL;
   return;
}

Graph::Graph(LPGraph InitialGraph)
{
   nVertices                 = 0;
   Edge                      = NULL;
   Labels                    = NULL;
   nLabels                   = 0;
   Cliques                   = NULL;
   CliquesDimens             = NULL;
   nCliques                  = 0;
   ConnectedComponents       = NULL;
   ConnectedComponentsDimens = NULL;
   nConnectedComponents      = 0;
   StarComp                  = NULL;
   TreeEdgeA                 = NULL;
   TreeEdgeB                 = NULL;
   nTreeEdges                = 0;
   nMss                      = 0;
   Mss                       = NULL;
   MssDimens                 = NULL;
   ordering                  = NULL;
   Separators                = NULL;
   SeparatorsDimens          = NULL;
   nSeparators               = 0;
   localord                  = NULL;
	
   ///////////////////////////////////////
   int i,j;
   InitGraph(InitialGraph->nVertices);	
   for(i=0;i<nVertices;i++)
   {	
      for(j=0;j<nVertices;j++)
      {
         Edge[i][j] = InitialGraph->Edge[i][j];
      }
   }
   nMss      = InitialGraph->nMss;
   MssDimens = new int[nMss];
   CheckPointer(MssDimens);
   memset(MssDimens,0,nMss*sizeof(int));
   for(i=0;i<nMss;i++)
   {
      MssDimens[i] = InitialGraph->MssDimens[i];		
   }	
   Mss  = new int*[nMss];
   CheckPointer(Mss);
   memset(Mss,0,nMss*sizeof(int*));
   for(i=0;i<nMss;i++)
   {
      Mss[i] = new int[MssDimens[i]];
      CheckPointer(Mss[i]);
      memset(Mss[i],0,MssDimens[i]*sizeof(int));
      for(j=0;j<MssDimens[i];j++)
      {
         Mss[i][j] = InitialGraph->Mss[i][j];
      }	
   }	
   return;
}	

Graph::~Graph()
{
   int i;
	
   for(i=0; i<nVertices; i++)
   {
      delete[] Edge[i];
      Edge[i] = NULL; 
   }
   delete[] Edge; Edge = NULL;
   delete[] Labels; Labels = NULL;
	
   for(i=0; i<nVertices; i++)
   {
      delete[] Cliques[i];
      Cliques[i] = NULL;
   }
   delete[] Cliques; Cliques = NULL;
   delete[] CliquesDimens; CliquesDimens = NULL;
  
   if((nConnectedComponents>0)&&(ConnectedComponents!=NULL))
   {
      for(i=0;i<nConnectedComponents;i++)
      {  
         delete[] ConnectedComponents[i];
         ConnectedComponents[i] = NULL;
      }
      delete[] ConnectedComponents;
      ConnectedComponents = NULL;
      delete[] ConnectedComponentsDimens;
      ConnectedComponentsDimens = NULL;
   }
   if(NULL!=StarComp)
   {
      for(i=0;i<nTreeEdges;i++)
      {
         delete[] StarComp[i];
         StarComp[i] = NULL;
      }
      delete[] StarComp; StarComp = NULL;
   }
   delete[] TreeEdgeA; TreeEdgeA = NULL;
   delete[] TreeEdgeB; TreeEdgeB = NULL;	
   for(i=0;i<nMss;i++)
   {
      delete[] Mss[i];
      Mss[i] = NULL;
   }	
   delete[] Mss; Mss = NULL;
   delete[] MssDimens; MssDimens = NULL;	
   delete[] ordering; ordering = NULL;
	
   for(i=0; i<nVertices; i++)
   {
      delete[] Separators[i];
      Separators[i] = NULL;
   }
   delete[] Separators; Separators = NULL;
   delete[] SeparatorsDimens; SeparatorsDimens = NULL;
   delete[] localord; localord = NULL;
   return;
}

//REMEMBER!!!
//THIS FUNCTION DOES NOT INITIALIZE THE MSS COMPONENTS
//Should Call Function InitGraphFromMss() to properly initialize the graph
void Graph::InitGraph(int n)
{
   int i, j;
	
   nVertices = n;
   //alloc the matrix of vertices
   Edge = new int*[nVertices];
   CheckPointer(Edge);
   memset(Edge,0,nVertices*sizeof(int*));
   for(i=0; i<n; i++)
   {
      Edge[i] = new int[nVertices];
      CheckPointer(Edge[i]);
      memset(Edge[i],0,nVertices*sizeof(int));
   }
   nLabels = 0;
   Labels = new int[nVertices];
   CheckPointer(Labels);
   memset(Labels,0,nVertices*sizeof(int));
   nCliques = 0;
   Cliques = new int*[nVertices];
   CheckPointer(Cliques);
   memset(Cliques,0,nVertices*sizeof(int*)); 
   for(i=0;i<n;i++)
   {
      Cliques[i] = new int[nVertices];
      CheckPointer(Cliques[i]);
      memset(Cliques[i],0,nVertices*sizeof(int));
   }	
   CliquesDimens = new int[nVertices];
   CheckPointer(CliquesDimens);
   memset(CliquesDimens,0,nVertices*sizeof(int));
   nTreeEdges = 0;
   TreeEdgeA = new int[nVertices];
   CheckPointer(TreeEdgeA);
   memset(TreeEdgeA,0,nVertices*sizeof(int));
   TreeEdgeB = new int[nVertices];
   CheckPointer(TreeEdgeB);
   memset(TreeEdgeB,0,nVertices*sizeof(int));
   ordering = new int[nVertices];
   CheckPointer(ordering);
   memset(ordering,0,nVertices*sizeof(int));
   Separators = new int*[nVertices];
   CheckPointer(Separators);
   memset(Separators,0,nVertices*sizeof(int*));
   for(i=0; i<n; i++)
   {
      Separators[i] = new int[nVertices];
      CheckPointer(Separators[i]);
      memset(Separators[i],0,nVertices*sizeof(int));
   }
   SeparatorsDimens = new int[nVertices];
   CheckPointer(SeparatorsDimens);
   memset(SeparatorsDimens,0,nVertices*sizeof(int));
   localord = new int[nVertices];
   CheckPointer(localord);
   memset(localord,0,nVertices*sizeof(int));
   return;
}	

int Graph::ReadGraph(char* sFileName)
{
	int n = 0;
	int i, j, k;
	
	FILE* in = fopen(sFileName, "r");
	if(NULL == in)
	{
		printf("Cannot open input file :: %s\n", sFileName);
		return 0;
	}
	
	fscanf(in, "%d", &n);
	if(n <= 0) return 0;
	InitGraph(n);
	fscanf(in, "%d", &n);
	if(n <= 0) return 0; //number of edges in the graph
	
	for(k=0; k<n; k++)
	{
		if(2 != fscanf(in, "%d %d", &i, &j))
		{
			printf("Error reading input file :: %s\n", sFileName);
			return 0;
		}		
		Edge[i-1][j-1] = 1;
		Edge[j-1][i-1] = 1;
	}
	
	fclose(in);
	return 1;
}

int Graph::ReadMss(char* sFileName)
{
	FILE* in;	
	int s;
	int i,j;
	
	if(NULL==(in=fopen(sFileName,"r")))
	{
		printf("Could not open file %s!!\n", sFileName);
		return 0;
	}	
	if(!fscanf(in,"%d",&s))
	{
		printf("File %s has an unknown format!\n", sFileName);
		fclose(in);
		return 0;
	}
	nMss = s;	
	MssDimens = new int[nMss];
	CheckPointer(MssDimens);	
	memset(MssDimens,0,nMss*sizeof(int));	
	Mss = new int*[nMss];
	CheckPointer(Mss);
	memset(Mss,0,nMss*sizeof(int*));	
        //read the mss one  by one
	for(i=0; i<nMss; i++)
	{
		if(!fscanf(in,"%d",&s))
		{
			printf("File %s has an unknown format!\n", sFileName);
			fclose(in);
			//clean the memory we allocated already
			for(j=0;j<i;j++)
			{
				delete[] Mss[i];Mss[i] = NULL;
			}
			delete[] Mss;Mss = NULL;
			delete[] MssDimens; MssDimens = NULL;
			nMss = 0;
			return 0;
		}
		MssDimens[i] = s;
		Mss[i] = new int[MssDimens[i]];
		CheckPointer(Mss[i]);
		memset(Mss[i],0,MssDimens[i]*sizeof(int));
		for(j=0;j<MssDimens[i];j++)
		{
			if(!fscanf(in,"%d",&s))
			{
				printf("File %s has an unknown format!\n",sFileName);
				fclose(in);
			     //clean the memory we allocated already
				for(j=0;j<i;j++)
				{
					delete[] Mss[i];Mss[i] = NULL;
				}
				delete[] Mss;Mss = NULL;
				delete[] MssDimens; MssDimens = NULL;
				nMss = 0;				
				return 0;
			}
			Mss[i][j] = s-1;//index begins at zero
		}		
          //the mss need to be sorted
		qsort((void*)Mss[i],MssDimens[i], sizeof(int), numeric);		
	}	
	fclose(in);
	return 1;
}	

int Graph::EdgeListToMss(int* edgelist, int s)
{
	int i,j,a,b;
	
	nMss = s;	
	MssDimens = new int[nMss];
	CheckPointer(MssDimens);	
	memset(MssDimens,0,nMss*sizeof(int));	
	Mss = new int*[nMss];
	CheckPointer(Mss);
	memset(Mss,0,nMss*sizeof(int*));	
        //read the mss one  by one
	for(i=0; i<nMss; i++)
	{
                a = edgelist[i];
                b = edgelist[nMss + i];
                if(a==b){
                  MssDimens[i] = 1;
  		  Mss[i] = new int[MssDimens[i]];
		  CheckPointer(Mss[i]);
		  memset(Mss[i],0,MssDimens[i]*sizeof(int));
                  Mss[i][0] = a-1;
                  //std::cout << "Mss" << i << "[0] = " << a-1 << std::endl;
                }
                else{
                  MssDimens[i] = 2;
  		  Mss[i] = new int[MssDimens[i]];
		  CheckPointer(Mss[i]);
		  memset(Mss[i],0,MssDimens[i]*sizeof(int));
                  Mss[i][0] = a-1;
                  Mss[i][1] = b-1;
		  qsort((void*)Mss[i],MssDimens[i], sizeof(int), numeric);		
                  //std::cout << "Mss" << i << "[0] = " << a-1 <<std::endl;
                  //std::cout << "Mss" << i << "[1] = " << b-1 <<std::endl;
                }	
	}	
	return 1;
}


void Graph::InitGraphFromMss()
{
	int i,j,k;
	int n = 0;
	//first determine the number of vertices
	//remember that the Mss are now sorted in ascending order
	for(i=0;i<nMss;i++)
	{
		n = max(n,Mss[i][MssDimens[i]-1]);
	}
	InitGraph(n+1); //since the index started at 0
	//now determine the adjacency matrix based on the Mss
	for(i=0;i<nMss;i++)
	{
		for(j=0;j<MssDimens[i];j++)
		{
			for(k=j+1;k<MssDimens[i];k++)
			{
				Edge[Mss[i][j]][Mss[i][k]] = 1;
				Edge[Mss[i][k]][Mss[i][j]] = 1;
			}	
		}	
	}
	return;
}	

void Graph::WriteInfo(FILE* out)
{
   int i, j;
	
   fprintf(out,"Minimal Sufficient Statistics\n");
   for(i=0;i<nMss;i++)
   {
      fprintf(out,"%d ::",MssDimens[i]);
      for(j=0;j<MssDimens[i];j++)
      {
         fprintf(out," %d",Mss[i][j]+1);
      }
      fprintf(out,"\n");
   }
   fprintf(out,"\n");
   for(i=0; i<nCliques; i++)
   {
      fprintf(out, "Clique %d :: ", i+1);
      for(j=0; j<CliquesDimens[i]; j++)
      {
         fprintf(out, " %d", Cliques[i][j]+1);
      }
      fprintf(out, "\n");
   }
   for(i=0;i<nConnectedComponents;i++)
   {
      fprintf(out,"Connected component %d :: ", i+1);
      for(j=0;j<ConnectedComponentsDimens[i];j++)
      {
         fprintf(out," %d",ConnectedComponents[i][j]+1);
      }
      fprintf(out, "\n");
   }
   for(i=0; i<nTreeEdges; i++)
   {
      fprintf(out, "Edge %d :: (%d, %d)\n", i+1, TreeEdgeA[i]+1, TreeEdgeB[i]+1);
   }
	
   for(i=0; i<nSeparators; i++)
   {
      fprintf(out, "Separator %d :: ", i+1);
      for(j=0; j<SeparatorsDimens[i]; j++)
      {
         fprintf(out, " %d", Separators[i][j]+1);
      }
      fprintf(out, "\n");
   }
   if(NULL!=StarComp)
   {
      for(i=0;i<nTreeEdges;i++)
      {
	 fprintf(out,"StarComp %d :: ",i);
         for(j=0;j<nVertices;j++)
	    fprintf(out," %d",StarComp[i][j]);
         fprintf(out,"\n");
      }
   }
   return;
}

void Graph::GenerateCliques(int label)
{
	int i, j, k, p, r;
	int n = nVertices;
	int* clique = new int[nVertices];
	CheckPointer(clique);
	int* LRound = new int[nVertices];
	CheckPointer(LRound);
	
	//clean memory
	memset(localord,0,nVertices*sizeof(int));
	memset(clique,0,nVertices*sizeof(int));
	memset(LRound,0,nVertices*sizeof(int));
	for(i=0;i<n;i++)
	{
		memset(Cliques[i],0,nVertices*sizeof(int));
	}
	memset(CliquesDimens,0,nVertices*sizeof(int));
	
	int v, vk;
	int PrevCard = 0;
	int NewCard;
	int s = -1;
	
	for(i = n-1; i>=0; i--)
	{
		NewCard = -1;
          //choose a vertex v...
		for(j=0; j<n; j++)
		{
      	     //test vertex j
			if((Labels[j] == label) && (LRound[j] == 0))
			{
				int maxj = 0;
				for(r=0; r<n; r++)
				{
					if(Labels[r] == label)
					{
						if(Edge[j][r] && LRound[r])
						{
							maxj++;
						}
					}
				}
				
				if(maxj > NewCard)
				{
					v = j;
					NewCard = maxj;
				}
			}
		}
		
		if(NewCard == -1)
		{
			break;
		}
		
		localord[v] = i;
		if(NewCard <= PrevCard)
		{
			//begin new clique
			s++;
			for(r=0; r<n; r++)
			{
				if(Labels[r] == label)
				{
					if(Edge[v][r] && LRound[r])
					{
						Cliques[s][CliquesDimens[s]] = r;
						CliquesDimens[s]++;
					}
				}
			}
			if(NewCard != 0)
			{
         	          //get edge to parent
				vk = Cliques[s][0];
				k  = localord[vk];
				for(r=1; r<CliquesDimens[s]; r++)
				{
					if(localord[Cliques[s][r]] < k)
					{
						vk = Cliques[s][r];
						k  = localord[vk];
					}
				}				
				p = clique[vk];
				TreeEdgeA[nTreeEdges] = s;
				TreeEdgeB[nTreeEdges] = p;
				nTreeEdges++;
			}
		}		
		clique[v] = s;
		Cliques[s][CliquesDimens[s]] = v;
		CliquesDimens[s]++;
		LRound[v] = 1;
		PrevCard = NewCard;
	}
	
	nCliques = s+1;
	
	delete[] clique;
	delete[] LRound;
	return;
}

int Graph::CheckCliques()
{
	int i, j, k;
	
	for(i=0; i<nCliques; i++)
	{
		for(j = 0; j<CliquesDimens[i]-1; j++)
		{
			for(k = j+1; k<CliquesDimens[i]; k++)
			{
				if(Edge[Cliques[i][j]][Cliques[i][k]] == 0)
				{
					return(-i-1);
				}
			}
		}
		qsort((void*)Cliques[i], CliquesDimens[i], sizeof(int), numeric);
	}
	
	return 1;
}

int Graph::IsClique(int* vect,int nvect)
{
   int i,j;
   int okay = 1;
   for(i=0;i<nvect;i++)
   {
      for(j=i+1;j<nvect;j++)
      {
         if(Edge[vect[i]][vect[j]]==0)
	 {
	    okay = 0; break;				
	 }
      }
      if(!okay) break;
   }
   return okay;
}	

int Graph::IsSubsetMss(int* vect,int nvect)
{
   int i,j,l;
   int okay;
   for(i=0;i<nMss;i++)
   {
      okay = 1;
      if(nvect>MssDimens[i])
      {
         okay=0; continue;
      }
      for(j=0;j<nvect;j++)
      {
	 int found=0;
	 for(l=0;l<MssDimens[i];l++)
	 {
	    if(vect[j]==Mss[i][l])
	    {
	       found = 1;
               break;
            } 
         }
         if(!found)
	 {
	    okay = 0;
            break;
         }  
      }
      if(okay) break;
   }
   return okay;
}

void Graph::GenerateSeparators()
{
	int i;
	int j, k;
	int FirstClique, SecondClique;
	int v;
	
	for(i=0; i<nTreeEdges; i++)
	{
		FirstClique = TreeEdgeA[i];
		SecondClique = TreeEdgeB[i];
		
		for(j=0; j<CliquesDimens[FirstClique]; j++)
		{
			v = Cliques[FirstClique][j];
			for(k=0; k<CliquesDimens[SecondClique]; k++)
			{
				if(v == Cliques[SecondClique][k])
				{
					Separators[i][SeparatorsDimens[i]] = v;
					SeparatorsDimens[i]++;
					break;
				}
			}
		}
		qsort((void*)Separators[i], SeparatorsDimens[i], sizeof(int), numeric);
	}	
	return;
}

void Graph::AttachLabel(int v, int label)
{
	int i;
	
	//only if v has not been labeled yet
	if(Labels[v] == 0)
	{
		Labels[v] = label;
		for(i=0; i<nVertices; i++)
		{
			if(Edge[v][i] == 1)
			{
				AttachLabel(i, label);
			}
		}
	}	
	return;
}

void Graph::GenerateLabels()
{
	int i;
	int NotFinished = 1;
	int label = 0;
	int v;
	
	memset(Labels,0,nVertices*sizeof(int));
	nLabels = 0;
	while(NotFinished)
	{
		v = -1;
		for(i=0; i<nVertices; i++)
		{
			if(Labels[i] == 0)
			{
				v = i;
				break;
			}
		}
		
		if(v == -1)
		{
			NotFinished = 0;
		}
		else
		{
			label++;
			AttachLabel(v, label);
		}
	}	
	nLabels = label;
	return;
}

int Graph::GenerateAllCliques()
{
	int i, j;
	int n = nVertices;
	int label;
	int nAssigned = 0;
	
     //Alloc Memory :: Begin
	int nAllCliques  = 0;
	int** AllCliques = new int*[n];
	CheckPointer(AllCliques);
	memset(AllCliques,0,n*sizeof(int*));
	for(i=0;i<n;i++)
	{
		AllCliques[i] = new int[n];
		CheckPointer(AllCliques[i]);
		memset(AllCliques[i],0,n*sizeof(int));
	}
	
	int* AllCliquesDimens = new int[n];
	CheckPointer(AllCliquesDimens);
	memset(AllCliquesDimens,0,n*sizeof(int));
	
	int nAllTreeEdges = 0;
	int* AllTreeEdgeA = new int[n];
	CheckPointer(AllTreeEdgeA);
	memset(AllTreeEdgeA,0,n*sizeof(int));
	int* AllTreeEdgeB = new int[n];
	CheckPointer(AllTreeEdgeB);
	memset(AllTreeEdgeB,0,n*sizeof(int));
	
	int** AllSeparators = new int*[n];
	CheckPointer(AllSeparators);
	memset(AllSeparators,0,n*sizeof(int*));
	for(i=0;i<n;i++)
	{
		AllSeparators[i] = new int[n];
		CheckPointer(AllSeparators[i]);
		memset(AllSeparators[i],0,n*sizeof(int));
	}
	int nAllSeparators = 0;	
	int* AllSeparatorsDimens = new int[n];
	CheckPointer(AllSeparatorsDimens);
	memset(AllSeparatorsDimens,0,n*sizeof(int));
	//Alloc Memory :: End
	
	//clean memory
	nCliques = 0;	
	for(i=0;i<n;i++)
	{		
		memset(Cliques[i],0,n*sizeof(int));
	}		
	memset(CliquesDimens,0,n*sizeof(int));
	nTreeEdges = 0;	
	memset(TreeEdgeA,0,n*sizeof(int));	
	memset(TreeEdgeB,0,n*sizeof(int));	
	memset(ordering,0,n*sizeof(int));	
	for(i=0; i<n; i++)
	{
		memset(Separators[i],0,n*sizeof(int));
	}	
	memset(SeparatorsDimens,0,n*sizeof(int));	
	
        //find all the connected components	
	GenerateLabels();
	
	for(label = 1; label<=nLabels; label++)
	{		
		GenerateCliques(label);
		if(CheckCliques() < 0)
		{
		  for(i=0; i<n; i++)
		    {
		      delete[] AllCliques[i];
		      AllCliques[i] = NULL;
		    }
		  delete[] AllCliques; AllCliques = NULL;
		  delete[] AllCliquesDimens; AllCliquesDimens = NULL;
		  delete[] AllTreeEdgeA; AllTreeEdgeA = NULL;
		  delete[] AllTreeEdgeB; AllTreeEdgeB = NULL;
		  for(i=0; i<n; i++)
		    {
		      delete[] AllSeparators[i];
		      AllSeparators[i] = NULL;
		    }
		  delete[] AllSeparators; AllSeparators = NULL;
		  delete[] AllSeparatorsDimens; AllSeparatorsDimens = NULL;
		  
		  return 0; //this is not a decomposable model
		}
		GenerateSeparators();		
          //store the newly generated cliques
		for(i=0; i<nTreeEdges; i++)
		{
			AllTreeEdgeA[nAllTreeEdges] = nAllCliques + TreeEdgeA[i];
			AllTreeEdgeB[nAllTreeEdges] = nAllCliques + TreeEdgeB[i];
			TreeEdgeA[i] = 0;
			TreeEdgeB[i] = 0;
			nAllTreeEdges++;
		}		
		for(i=0; i<nCliques; i++)
		{
			for(j=0; j<CliquesDimens[i]; j++)
			{
				AllCliques[nAllCliques][j] = Cliques[i][j];
				Cliques[i][j] = 0;
			}
			AllCliquesDimens[nAllCliques] = CliquesDimens[i];
			CliquesDimens[i] = 0;
			nAllCliques++;
		}
		nCliques = 0;
		
		for(i=0; i<nTreeEdges; i++)
		{
			for(j=0; j<SeparatorsDimens[i]; j++)
			{
				AllSeparators[nAllSeparators][j] = Separators[i][j];
				Separators[i][j] = 0;
			}
			AllSeparatorsDimens[nAllSeparators] = SeparatorsDimens[i];
			SeparatorsDimens[i] = 0;
			nAllSeparators++;
		}
		/*
		//add an extra (null) separator between two connected components
		if(label<nLabels)
		{	
			AllSeparatorsDimens[nAllSeparators] = 0;
			nAllSeparators++;
		}
                */	
          //clean memory
		nSeparators = 0;		
		nTreeEdges = 0;
		
		//printf("Perfect ordering :: ");
		int partialAssigned = 0;
		for(i=0; i<n; i++)
		{
			//printf("%d ",localord[i]);
			if(Labels[i] == label)
			{	
				ordering[i] = localord[i] - nAssigned;
				partialAssigned++;
			}
		}
						//printf("\n");
		nAssigned += partialAssigned;
	}	
	for(i=0; i<nAllCliques; i++)
	{
		for(j=0; j<AllCliquesDimens[i]; j++)
		{
			Cliques[nCliques][j] = AllCliques[i][j];
		}
		CliquesDimens[nCliques] = AllCliquesDimens[i];
		nCliques++;
	}
	
	for(i=0; i<nAllTreeEdges; i++)
	{
		TreeEdgeA[nTreeEdges] = AllTreeEdgeA[i];
		TreeEdgeB[nTreeEdges] = AllTreeEdgeB[i];
		nTreeEdges++;
	}
	
	for(i=0; i<nAllSeparators; i++)
	{
		for(j=0; j<AllSeparatorsDimens[i]; j++)
		{
			Separators[nSeparators][j] = AllSeparators[i][j];
		}
		SeparatorsDimens[nSeparators] = AllSeparatorsDimens[i];
		nSeparators++;
	}
        //free memory
	for(i=0; i<n; i++)
	{
		delete[] AllCliques[i];
		AllCliques[i] = NULL;
	}
	delete[] AllCliques; AllCliques = NULL;
	delete[] AllCliquesDimens; AllCliquesDimens = NULL;
	delete[] AllTreeEdgeA; AllTreeEdgeA = NULL;
	delete[] AllTreeEdgeB; AllTreeEdgeB = NULL;
	for(i=0; i<n; i++)
	{
		delete[] AllSeparators[i];
		AllSeparators[i] = NULL;
	}
	delete[] AllSeparators; AllSeparators = NULL;
	delete[] AllSeparatorsDimens; AllSeparatorsDimens = NULL;
	return 1;
}

int Graph::SearchVertex()
{
	int x, u, v;
	int okay;
	int* sxAdj = new int[nVertices];
	CheckPointer(sxAdj);
	memset(sxAdj,0,nVertices*sizeof(int));
	
	for(x=0;x<nVertices;x++)
	{
		memmove(sxAdj,Edge[x],nVertices*sizeof(int));	
		sxAdj[x] = 1;
		okay = 1;
		for(u=0;u<nVertices;u++)
		{
			if((u!=x)&&(Edge[x][u]==1))
			{
				sxAdj[u] = 0; //we take u out
				for(v=u+1;v<nVertices;v++)
				{
					if((v!=x)&&(Edge[x][v]==1)&&(Edge[u][v]==0))
					{
						sxAdj[v] = 0;//we take v out
						SectionGraph sgraph(this,sxAdj);
						okay = sgraph.IsChain(u,v);
						sxAdj[v] = 1;//now put v back in the adjacency list of x
					}
					if(!okay) break;
				}
				sxAdj[u] = 1; //we put u back
			}
			if(!okay) break;
		}
		if(okay) break;
	}
	delete[] sxAdj;
	if(x==nVertices) x = -1;
	return x;
}

int Graph::IsDecomposable()
{
   return GenerateAllCliques();
}

void Graph::InitConnectedComponents()
{
   int i,label;
   if(ConnectedComponents !=NULL){
     for(i=0; i< nConnectedComponents; i++){
       if(ConnectedComponents[i]!=NULL) delete[] ConnectedComponents[i];}
     delete[] ConnectedComponents;}
   if(ConnectedComponentsDimens!=NULL) delete[] ConnectedComponentsDimens;

   nConnectedComponents=nLabels;
   
   ConnectedComponents = new int*[nConnectedComponents];
   CheckPointer(ConnectedComponents);

   ConnectedComponentsDimens = new int[nConnectedComponents];
   CheckPointer(ConnectedComponentsDimens);
   for(label=1;label<=nLabels;label++)
   {
      //count the number of vertices being labeled with label
      int count=0;
      for(i=0;i<nVertices;i++)
      {
         if(Labels[i]==label) count++;  
      }
      //printf("label = %d :: count = %d\n",label,count);
      ConnectedComponentsDimens[label-1]=count;
      ConnectedComponents[label-1] = new int[count];
      CheckPointer(ConnectedComponents[label-1]);
      count=0;
      for(i=0;i<nVertices;i++)
      {
         if(Labels[i]==label)
	 {
	    ConnectedComponents[label-1][count]=i;
            count++;
         }  
      }         
   }
   return;
}

void Graph::GetMPSubgraphs()
{
   int i,j,s;
	
   if(IsDecomposable()) return;//easy task if the graph is decomposable
   //if not, generate the minimal fill-in graph
   LPGraph gfill = MakeFillInGraph(this);
   if(!gfill->IsDecomposable())
   {
      printf("The fill-in graph is not decomposable!\n Something is wrong.\n");
      exit(1);
   }	
   //gfill->WriteInfo(stdout);
   //we clean the memory a bit, just to be on the safe side
   nCliques = nSeparators = 0;	
   for(i=0;i<nVertices;i++)
   {		
      memset(Cliques[i],0,nVertices*sizeof(int));
      memset(Separators[i],0,nVertices*sizeof(int));
   }		
   memset(CliquesDimens,0,nVertices*sizeof(int));
   memset(SeparatorsDimens,0,nVertices*sizeof(int));
   nTreeEdges = 0;	
   memset(TreeEdgeA,0,nVertices*sizeof(int));	
   memset(TreeEdgeB,0,nVertices*sizeof(int));	
   memset(ordering,0,nVertices*sizeof(int));					
   //////////////////////////////////////////////////////////
   //done cleaning memory                                  //
   //////////////////////////////////////////////////////////	
   int* UsedEdge = new int[gfill->nTreeEdges]; CheckPointer(UsedEdge);
   //mark the edges as "not used"
   memset(UsedEdge,0,gfill->nTreeEdges*sizeof(int));   
   int* MarkC = new int[gfill->nCliques]; CheckPointer(MarkC); 
   memset(MarkC,0,gfill->nCliques*sizeof(int));
   int* MarkS = new int[gfill->nSeparators]; CheckPointer(MarkS);
   memset(MarkS,0,gfill->nSeparators*sizeof(int)); 
   
   //printf("nTreeEdges = %d\n",gfill->nTreeEdges);
   ////////////////////////////////////////////////////////////   
   while(1)
   {
      //identify a terminal clique Cj
      int edg;
      int Ci;
      int Cj;
      for(edg=0;edg<gfill->nTreeEdges;edg++)
      {
	 //if we already used that edge, go to the next one
	 if(UsedEdge[edg]) continue;
         Ci = gfill->TreeEdgeB[edg];
         Cj = gfill->TreeEdgeA[edg];
         //printf("Cj = %d\n",Cj+1);         
         int foundterminal = 1;
         for(i=0;i<gfill->nTreeEdges;i++)
	 {
	    if(UsedEdge[i]) continue;
            if(gfill->TreeEdgeB[i]==Cj)
	    {
	       foundterminal=0;
               break;
            } 
         }
         if(foundterminal) break;
      }
      if(edg==gfill->nTreeEdges) break;
      //printf("Cj = %d :: Ci = %d\n",Cj+1,Ci+1);
      //mark the edge as used
      UsedEdge[edg]=1;
      //Step 4
      if(IsClique(gfill->Separators[edg],
                  gfill->SeparatorsDimens[edg]))
      {
	 //printf("%d is clique\n",edg+1);
         MarkC[Cj]  = 1;
         MarkS[edg] = 1;  
      }
      else
      {
	 MarkC[Cj] = -1;
         //combine the delta sets associated with Ci and Cj
         int  len1 = gfill->CliquesDimens[Ci]+
                     gfill->CliquesDimens[Cj];
         int  len2 = 0;
         int* buffer1 = new int[len1]; CheckPointer(buffer1);		
         int* buffer2 = new int[len1]; CheckPointer(buffer2);
         len1 = 0;
         for(i=0;i<gfill->CliquesDimens[Ci];i++)
	 {
	    buffer1[len1] = gfill->Cliques[Ci][i];
            len1++;   
         }
         for(i=0;i<gfill->CliquesDimens[Cj];i++)
	 {
	    buffer1[len1] = gfill->Cliques[Cj][i];
            len1++;
	 }
         qsort((void*)buffer1,len1,sizeof(int),numeric);
         buffer2[len2]=buffer1[0];
         for(i=0;i<len1;i++)
         {
            if(buffer2[len2]<buffer1[i])
	    {
	       len2++;
               buffer2[len2]=buffer1[i];
	    }	
	 }
         len2++;
         for(i=0;i<len2;i++)
	 {
	    gfill->Cliques[Ci][i] = buffer2[i];
         } 
         gfill->CliquesDimens[Ci] = len2;    
         ///////////////
         delete[] buffer1;
         delete[] buffer2;
      }
   }
   for(i=0;i<gfill->nCliques;i++)
   {
      if(MarkC[i]==-1) continue;
      for(j=0;j<gfill->CliquesDimens[i];j++)
      {
	 Cliques[nCliques][j] = gfill->Cliques[i][j];
      }
      CliquesDimens[nCliques] = gfill->CliquesDimens[i];
      nCliques++;
   }
   for(i=0;i<gfill->nSeparators;i++)
   {
      if(MarkS[i]==0) continue;
      for(j=0;j<gfill->SeparatorsDimens[i];j++)
      {
	 Separators[nSeparators][j] = gfill->Separators[i][j];
      }
      SeparatorsDimens[nSeparators] = gfill->SeparatorsDimens[i];
      nSeparators++;
   }
   //////////////////////////////////////////////////
   delete[] MarkS;
   delete[] MarkC;
   delete[] UsedEdge;
   delete gfill;	
   return;
}


void Graph::FindCliqueTree()
{
   int i,j,k;
   LPGraph tempGraph = new Graph; CheckPointer(tempGraph);
   tempGraph->InitGraph(nVertices);
   for(i=0;i<nCliques;i++)
   {
      for(j=0;j<CliquesDimens[i]-1;j++)
      {
         for(k=j+1;k<CliquesDimens[i];k++)
	 {
	    tempGraph->Edge[Cliques[i][j]][Cliques[i][k]] = 1;
            tempGraph->Edge[Cliques[i][k]][Cliques[i][j]] = 1;
         }
      }
   }
   if(!tempGraph->IsDecomposable())
   {
      printf("Something is very wrong in FindCliqueTree. Program exits...\n");
      exit(1);
   }
   ///////////////////////////////////////////////////////////////
   //the ordering of the cliques/mp-subgraphs might have changed//
   //so we need to copy the cliques from tempGraph              //
   //we do the same with the clique tree that was generated     //
   ///////////////////////////////////////////////////////////////
   nCliques = tempGraph->nCliques; //this should be redundant
   for(i=0;i<nCliques;i++)
   {  
      CliquesDimens[i] = tempGraph->CliquesDimens[i];
      for(j=0;j<CliquesDimens[i];j++)
      {
         Cliques[i][j] = tempGraph->Cliques[i][j];
      }
   }
   nSeparators = tempGraph->nSeparators; //this should also be redundant
   for(i=0;i<nSeparators;i++)
   {  
      SeparatorsDimens[i] = tempGraph->SeparatorsDimens[i];
      for(j=0;j<SeparatorsDimens[i];j++)
      {
         Separators[i][j] = tempGraph->Separators[i][j];
      }
   }
   nTreeEdges = tempGraph->nTreeEdges;
   for(i=0;i<nTreeEdges;i++)
   {
      TreeEdgeA[i] = tempGraph->TreeEdgeA[i];
      TreeEdgeB[i] = tempGraph->TreeEdgeB[i];
   }
   ////////////////////////////////////////////////////
   //initialize the connected components of the graph//
   ////////////////////////////////////////////////////
   tempGraph->InitConnectedComponents();
   int* LabelEdges = new int[nTreeEdges]; CheckPointer(LabelEdges);
   for(i=0;i<nTreeEdges;i++)
   {
      LabelEdges[i] = tempGraph->Labels[Cliques[TreeEdgeB[i]][0]];
      printf("Edge %d :: label %d\n",
             i+1,LabelEdges[i]); 
   }
   int* RootConComp = new int[tempGraph->nLabels]; CheckPointer(RootConComp);
   int* LeafConComp = new int[tempGraph->nLabels]; CheckPointer(LeafConComp);
   //printf("Concomp = %d\n",tempGraph->nLabels);
   for(i=0;i<tempGraph->nLabels;i++)
   {
      ////////////////////////////////////////////
      //identify the root of the clique tree for//
      //the i-th connected component            //
      ////////////////////////////////////////////
      int edg;
      int Ci;
      int Cj;
      for(edg=0;edg<nTreeEdges;edg++)
      {
	 if(LabelEdges[edg]!=i+1) continue;
         Ci = TreeEdgeB[edg];
         Cj = TreeEdgeA[edg];
         //printf("Cj = %d\n",Cj+1);         
         int foundroot = 1;
         for(j=0;j<nTreeEdges;j++)
	 {
	    if(LabelEdges[j]!=i+1) continue;
            if(TreeEdgeA[j]==Ci)
	    {
	       foundroot=0;
               break;
            } 
         }
         if(foundroot)
         {
            RootConComp[i] = Ci;
            break;
         }
      }
      for(edg=0;edg<nTreeEdges;edg++)
      {
	 if(LabelEdges[edg]!=i+1) continue;
         Ci = TreeEdgeB[edg];
         Cj = TreeEdgeA[edg];
         //printf("Cj = %d\n",Cj+1);         
         int foundterminal = 1;
         for(j=0;j<nTreeEdges;j++)
	 {
	    if(LabelEdges[j]!=i+1) continue;
            if(TreeEdgeB[j]==Cj)
	    {
	       foundterminal=0;
               break;
            } 
         }
         if(foundterminal)
         {
            LeafConComp[i] = Cj;
            break;
         }
      }
      printf("Component %d :: root = %d :: leaf = %d\n",
             i,RootConComp[i],LeafConComp[i]);  
   }
   ////////////////////////////////////////////////////////
   //add the edges that connect the trees associated with//
   //the connected components of the initial graph       //
   //instead of having a forest (a set of trees) we      //
   //will obtain one big tree connecting all cliques     //
   ////////////////////////////////////////////////////////
   for(i=0;i<tempGraph->nLabels-1;i++)
   {
      TreeEdgeA[nTreeEdges] = RootConComp[i];
      TreeEdgeB[nTreeEdges] = LeafConComp[i+1];
      nTreeEdges++;
      SeparatorsDimens[nSeparators] = 0;
      nSeparators++;
   }   
   /////////////////////////////
   //initialize the tree graph//
   /////////////////////////////
   LPGraph treeGraph = new Graph; CheckPointer(treeGraph);
   treeGraph->InitGraph(nCliques);
   for(i=0;i<nTreeEdges;i++)
   {
      treeGraph->Edge[TreeEdgeA[i]][TreeEdgeB[i]] = 1;
      treeGraph->Edge[TreeEdgeB[i]][TreeEdgeA[i]] = 1;
   }
   //////////////////////////////////////////////
   //alocate the memory for the star components//
   //////////////////////////////////////////////
   StarComp = new int*[nTreeEdges]; CheckPointer(StarComp);
   for(i=0;i<nTreeEdges;i++)
   {
      StarComp[i] = new int[nVertices]; CheckPointer(StarComp[i]);
      memset(StarComp[i],0,nVertices*sizeof(int)); 
   }
   /////////////////////////////////////////////////////////
   //take out one edge at a time and get the two connected//
   //components of the resulting tree graph               //
   /////////////////////////////////////////////////////////
   for(i=0;i<nTreeEdges;i++)
   {
      ///////////////////////////
      //eliminate the i-th edge//
      ///////////////////////////
      treeGraph->Edge[TreeEdgeA[i]][TreeEdgeB[i]] = 0;
      treeGraph->Edge[TreeEdgeB[i]][TreeEdgeA[i]] = 0;
      /////////////////////////////////
      //find the connected components//
      /////////////////////////////////
      treeGraph->GenerateLabels();
      //////////////////////////////////////////////////////////////////
      //mark with a "1" the vertices in the first connected component,//
      //with a "3" the vertices in the separator and with a "2"       //
      //the vertices in the second connected component                //
      //////////////////////////////////////////////////////////////////
      for(j=0;j<nCliques;j++)
      {
	 for(k=0;k<CliquesDimens[j];k++)
	    StarComp[i][Cliques[j][k]] = treeGraph->Labels[j];
      }
      for(k=0;k<SeparatorsDimens[i];k++)
      {
	 StarComp[i][Separators[i][k]] = 3;
      }
      ///////////////////////////////////
      //put back the edge we eliminated//
      ///////////////////////////////////
      treeGraph->Edge[TreeEdgeA[i]][TreeEdgeB[i]] = 1;
      treeGraph->Edge[TreeEdgeB[i]][TreeEdgeA[i]] = 1;      
   }
   ////////////////  
   //clean memory//
   ////////////////
   delete treeGraph;
   delete[] LeafConComp;
   delete[] RootConComp;
   delete[] LabelEdges;
   delete tempGraph;
   return;
}
//class Graph::Ends

//class SectionGraph::Begins
SectionGraph::SectionGraph(LPGraph InitialGraph,int* velim) : Graph(InitialGraph)
{
	int i,j;
	
	Eliminated = new int[nVertices];
	CheckPointer(Eliminated);
	memset(Eliminated,0,nVertices*sizeof(int));
	nEliminated = 0;
	for(i=0;i<nVertices;i++)
	{
		if(velim[i])
		{	
			Eliminated[i] = 1;
			nEliminated++;
		}	
	}
	//delete all the edges corresponding to the vertices
	//we eliminated
	for(i=0;i<nVertices;i++)
	{		
		if(Eliminated[i])
		{
			for(j=0;j<nVertices;j++)
			{
				if(1==Edge[i][j])
				{
					Edge[i][j] = Edge[j][i] = 0;
				}	
			}	
		}	
	}
	return;
}

SectionGraph::~SectionGraph()
{
	delete[] Eliminated;
	nEliminated = 0;
	return;
}	

int SectionGraph::IsChain(int u,int v)
{
	if(nLabels==0)
	{	
		GenerateLabels();
	}	
	if(Eliminated[u] || Eliminated[v])
	{
		printf("One of the vertices %d,%d has been eliminated...\n",u,v);
		exit(1);
	}		
	if(Labels[u]==Labels[v]) return 1;
	return 0;
}	
//class SectionGraph::Ends

//class EliminationGraph::Begins
EliminationGraph::EliminationGraph(LPGraph InitialGraph,int vertex) : Graph(InitialGraph)
{	
	Eliminated = new int[nVertices];
	CheckPointer(Eliminated);
	memset(Eliminated,0,nVertices*sizeof(int));
	nEliminated = 0;
	EliminateVertex(vertex);
	return;
}

EliminationGraph::~EliminationGraph()
{
	delete[] Eliminated;
	nEliminated = 0;
	return;
}

void EliminationGraph::EliminateVertex(int x)
{
	int i,j;
	
	//adding edges in Def(Adj(x)) so that Adj(x) becomes a clique
	for(i=0;i<nVertices;i++)
	{
		if((i!=x)&&(!Eliminated[i])&&(Edge[x][i]==1))
		{
			for(j=i+1;j<nVertices;j++)
			{
				if((j!=x)&&(!Eliminated[j])&&(Edge[x][j]==1)&&(Edge[i][j]==0))
				{
					Edge[i][j] = Edge[j][i] = 1;
				}	
			}	
		}	
	}	
	
	//eliminate all edges incident to x
	for(i=0;i<nVertices;i++)
	{
		if((i!=x)&&(!Eliminated[i])&&(Edge[x][i]==1))
		{
			Edge[x][i] = Edge[i][x] = 0;
		}	
	}	
	
	//eliminate vertex x
	Eliminated[x] = 1;
	nEliminated++;
	return;
}

int EliminationGraph::SearchVertex()
{
	int x, u, v;
	int okay;
	int* sxAdj = new int[nVertices];
	CheckPointer(sxAdj);
	memset(sxAdj,0,nVertices*sizeof(int));
	
	for(x=0;x<nVertices;x++)
	{
		if(Eliminated[x]) continue;
		memmove(sxAdj,Edge[x],nVertices*sizeof(int));	
		sxAdj[x] = 1;
		okay = 1;
		for(u=0;u<nVertices;u++)
		{
			if(Eliminated[u]) continue;
			if((u!=x)&&(Edge[x][u]==1))
			{
				sxAdj[u] = 0; //we take u out
				for(v=u+1;v<nVertices;v++)
				{
					if(Eliminated[v]) continue;
					if((v!=x)&&(Edge[x][v]==1)&&(Edge[u][v]==0))
					{
						sxAdj[v] = 0;//we take v out
						SectionGraph sgraph(this,sxAdj);
						okay = sgraph.IsChain(u,v);
						sxAdj[v] = 1;//now put v back in the adjacency list of x
					}
					if(!okay) break;
				}
				sxAdj[u] = 1; //we put u back
			}
			if(!okay) break;
		}
		if(okay) break;
	}
	delete[] sxAdj;
	if(x==nVertices) x = -1;
	return x;
}
//class EliminationGraph::Ends

LPGraph MakeFillInGraph(LPGraph graph)
{
	int u,v;
	int i,j;
	
	LPGraph gfill = new Graph(graph);
	CheckPointer(gfill);
	//if the graph is decomposable, there is no need to do anything
	if(gfill->IsDecomposable()) return gfill;
	
	int v1 = gfill->SearchVertex();
	//printf("v1 = %d\n",v1);
	//add edges to Def(Adj(x)) so that Adj(x) becomes a clique	
	for(u=0;u<gfill->nVertices;u++)
	{
		if(gfill->Edge[v1][u]==1)
		{
			for(v=u+1;v<gfill->nVertices;v++)
			{
				if((gfill->Edge[v1][v]==1)&&(gfill->Edge[u][v]==0))
				{
					gfill->Edge[v][u] = gfill->Edge[u][v] = 1;
					//printf("u = %d, v = %d\n",u,v);
				}	
			}	
		}	
	}		
	EliminationGraph egraph(graph,v1);
	for(i=1;i<graph->nVertices-1;i++)
	{
		v1 = egraph.SearchVertex();
		//printf("v1 = %d\n",v1);
		for(u=0;u<egraph.nVertices;u++)
		{
			if(egraph.Eliminated[u]) continue;
			if(egraph.Edge[v1][u]==1)
			{
				for(v=u+1;v<egraph.nVertices;v++)
				{
					if(egraph.Eliminated[v]) continue;
					if((egraph.Edge[v1][v]==1)&&(egraph.Edge[u][v]==0))
					{
						gfill->Edge[v][u] = gfill->Edge[u][v] = 1;
						//these are the edges that are added
						//to the initial graph
						//printf("u = %d, v = %d\n",u,v);
					}	
				}	
			}	
		}
		egraph.EliminateVertex(v1);
	}	
	return gfill;
}


void Graph::WriteInfo1(FILE* out)
{
   int i, j;
	
   fprintf(out,"\n");
   for(i=0; i<nCliques; i++)
   {
      fprintf(out, "Clique %d :: ", i+1);
      for(j=0; j<CliquesDimens[i]; j++)
      {
         fprintf(out, " %d", Cliques[i][j]+1);
      }
      fprintf(out, "\n");
   }
   	
   for(i=0; i<nSeparators; i++)
   {
      fprintf(out, "Separator %d :: ", i+1);
      for(j=0; j<SeparatorsDimens[i]; j++)
      {
         fprintf(out, " %d", Separators[i][j]+1);
      }
      fprintf(out, "\n");
   }
   if(NULL!=StarComp)
   {
      for(i=0;i<nTreeEdges;i++)
      {
	 fprintf(out,"StarComp %d :: ",i);
         for(j=0;j<nVertices;j++)
	    fprintf(out," %d",StarComp[i][j]);
         fprintf(out,"\n");
      }
   }
   return;
}


/////////////////////////////////////////////////////////////////////
/*
int main()
{
   int i, j, k;
   LPGraph graph = new Graph;
	
   graph->ReadMss("tarjan.dat");
   graph->InitGraphFromMss();
   //graph->WriteInfo(stdout);
   for(i=0; i<graph->nVertices; i++)
   {
      for(j=0; j<graph->nVertices; j++)
      {
         printf("%d  ", graph->Edge[i][j]);
      }
      printf("\n");
   }
   printf("\n");	
   graph->GetMPSubgraphs();
   graph->FindCliqueTree();
   graph->WriteInfo(stdout);
   delete graph;
   return 1;
}
*/








