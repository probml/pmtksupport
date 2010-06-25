#include "matrix.h"
typedef class Graph* LPGraph;

class Graph
{
public:
//data members
	int nVertices; //number of vertices in the graph
	int** Edge;  //matrix containing the edges of the graph
	int* Labels; //identifies the connected components of the graph
	int nLabels; //number of labels or connected components
	int** Cliques; //storage for cliques
	int* CliquesDimens; //number of vertices in each clique
	int nCliques; //number of cliques
       
public:
        int** ConnectedComponents;
        int* ConnectedComponentsDimens;
        int nConnectedComponents; //number of connected components
public:
        int** StarComp;
        
public:	
	int* TreeEdgeA; //edges of the clique tree
	int* TreeEdgeB;
	int nTreeEdges; //number of edges in the generated clique tree
public:
	int   nMss; //the number of MSS that define the graph
	int** Mss; //storage for the MSS
	int*  MssDimens; //number of vertices in each MSS
public:	
	int* ordering;
	int** Separators; //storage for separators
	int* SeparatorsDimens;
	int nSeparators;
private:
	int* localord;
	
        //methods
public:
	Graph(); //constructor
	Graph(LPGraph InitialGraph); //constructor 
	~Graph(); //destructor	
public:	
	int  SearchVertex(); //identifies the next vertex to be eliminated
	
public:
        int  EdgeListToMss(int* edgelist, int s); //read the MSS from an sx2 edge list
	int  ReadMss(char* sFileName); //read the MSS from file
	void InitGraphFromMss(); //initialize the graph based on the MSS
        void InitConnectedComponents();
public:	
	//the MSS (Minimal Sufficient Statistics) are the maximal cliques for our graph
	void InitGraph(int n);
	int  ReadGraph(char* sFileName);
	void WriteInfo(FILE* out);
	void WriteInfo1(FILE* out);
        void GenerateCliques(int label);
	int  CheckCliques(); //checks whether each generated component
                             //is complete in the given graph
	int  IsClique(int* vect,int nvect); //checks if the vertices in vect
	                                    //form a clique in our graph
        int  IsSubsetMss(int* vect,int nvect);
	void GenerateSeparators();
	void AttachLabel(int v, int label);
	void GenerateLabels();	
	int  GenerateAllCliques();
	int  IsDecomposable();
	void GetMPSubgraphs(); 
        //if the graph is decomposable, the mp-subgraphs will be the maximal cliques (the MSS)
	//otherwise, the minimum fill-in graph is generated
	//the mp-subgraphs will be stored in the Clique arrays
        void FindCliqueTree();
        //this method should be called after calling GetMPSubgraphs(),
        //to init the clique tree if the graph is not decomposable
};

//////////////////////////////////////////////////////////////////////

typedef class SectionGraph* LPSectionGraph;

class SectionGraph : public Graph
{
public:
	int* Eliminated; //shows which vertices were eliminated from
	//the initial graph
	int nEliminated; //number of vertices we eliminated
	
	//methods
public:
	SectionGraph(LPGraph InitialGraph,int* velim); //constructor
	~SectionGraph(); //destructor
	
public:
	int IsChain(int u,int v);//see if there is a chain between u and v
	//or, equivalently, checks if u and v are in the same connected component
};

////////////////////////////////////////////////////////////////////////

typedef class EliminationGraph* LPEliminationGraph;

class EliminationGraph : public Graph
{
public:
	int* Eliminated; //shows which vertices were eliminated from
	//the initial graph
	int nEliminated; //number of vertices we eliminated
	
	//methods
public:
	EliminationGraph(LPGraph InitialGraph,int vertex); //constructor
	~EliminationGraph(); //destructor
public:	
	int  SearchVertex(); //identify a vertex to be eliminated		
public:
	void EliminateVertex(int x); //eliminates an extra vertex
};

//////////////////////////////////////////////////////////////////////////

//constructs the minimum fill-in graph for a nondecomposable graph
LPGraph MakeFillInGraph(LPGraph graph);

