#include <vector>
#include <iostream>
#include <cmath> 
#include <string>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <iomanip>

using namespace std;

#include "Interface.h"

// Declaration of a struct named IE
IE ie;

// Declaration of spaces and time variables, diffusion
double Nx1_pom;
double Nx2_pom;
double Nx3_pom;
double t_0;
double t_N;
double dt;
double diffusion_c;
int count; 


// Function to generate a random number between 0 and 1
double epsilon() {
    const int c = 10000;
    int n = rand() % c;
    return static_cast<double>(n) / static_cast<double>(c - 1);
}

// Function to generate a random number that is either 1 or -1
int randomPlusMinusOne() {
    // Generate a random number: 0 or 1
    int randomValue = rand() % 2;

    // Return 1 if the random value is 0, otherwise return -1
    return (randomValue == 0) ? 1 : -1;
}

void graph1_branch(int step, int level, int parent_index, double A) {
	// Base case: If the current level exceeds the maximum level, return
    if (step > level+1) return; 
    // Extract coordinates of the parent node
    double x1 = ve1.V[parent_index].x1;
    double x2 = ve1.V[parent_index].x2;
    double x3 = ve1.V[parent_index].x3;
    // Generate a random number of child nodes to create
    int number = rand() % (step + 3) + 2;
 	// Loop to create child nodes
    for (int i = 0; i < number; i++) {
        // Calculate the cross-section of new edge ei
        double A_ei = A / number;
        // Generate random coordinates for the new node
        double x1_new = x1 + 0.5 + 0.5 * epsilon() / (level + 1);
        double x2_new = x2 + randomPlusMinusOne() * 4 * epsilon() / (level + 2);
        double x3_new = x3 + randomPlusMinusOne() * 2 * epsilon() / (level + 1);
        // Adjust x3 coordinate for a specific case
        if (i == 1) {
            x3_new = x3 + 1 * epsilon() / (level + 5); 
        }
        // Check if the new node is within the defined bounds of extravascular environment
        if (x1_new >= ee.x1_0 && x1_new <= ee.x1_N && x2_new >= ee.x2_0 && x2_new <= ee.x2_N && x3_new >= ee.x3_0 && x3_new <= ee.x3_N) {
            // Add the new node to the graph
            int v_i= ve1.addNode(x1_new, x2_new, x3_new, 0, 0, 0);
            // Add edge e_i connecting v_i and parent_node
            ve1.addEdge(parent_index, v_i, A_ei);
            // Recursively build the graph from the new node
            graph1_branch(step+1,level,v_i,A_ei);
        }
    }
}
   
// Function to create a graph with a defect
void graph1(double d, int level, double reflection_point = 7.) { 
    // Number of step
    int step = 0;
    // Seed the random number generator
    srand(11);
    // Number of entry nodes in x2 and x3 axis
    int number_x2 = 3;
    int number_x3 = 3;
    // Set x1 coordinates for vertices x1_in and x1_1 of edge e1
    double x1_in = 0;
    double x1_1 = 0.5;
    // Initial area of edge e1
    double A_e1 = 1;
	// Coordinates x2 and x3 of entry nodes(EN)
    double EN_x2[3] = {2.5, 5, 7.5};
    double EN_x3[3] = {2.5, 5, 7.5};
	// Integer for indexing entry and help nodes	
    int count = 0;
    // Creating vascular bed until step = level + 1
    for (int k = 0; k < number_x2; k++) {
        for (int l = 0; l < number_x3; l++) {
            // Add entry and help nodes to the graph
            ve1.EntryNode[count] = ve1.addNode(x1_in, EN_x2[k], EN_x3[l],0,0,Dirichlet);
            ve1.HelpNode[count] = ve1.addNode(x1_1,EN_x2[k],EN_x3[l],0,0,0);
            step = 1;
            // Add adge e1 
            ve1.addEdge(ve1.EntryNode[count],ve1.HelpNode[count],A_e1);
            // Build the graph starting from the help node with level
            graph1_branch(step,level,ve1.HelpNode[count],A_e1);
            count++;
        }
    } 
    // Step - mirroring
    // Build matrix of adjacent matrix
    ve1.buildAdjacent();
    // Current number of vertices and edges
    int mNodes = ve1.V.size();
    int mEdges = ve1.E.size();
    // Add reflected nodes and edges to the graph
    for (int i = 0; i < mNodes; i++)
        ve1.addNode(1.5*reflection_point-ve1.V[i].x1,ve1.V[i].x2,ve1.V[i].x3);
    
    for (int i = 0; i < mEdges; i++)
        int ei = ve1.addEdge(ve1.E[i].v1+mNodes,ve1.E[i].v2+mNodes,ve1.E[i].A); 

    // Step - connecting
    for (int i = 0; i < mNodes; i++){
        // Conditions for adding edges in step connecting: 1 neighbour and it can not be v1_1
        if (ve1.V[i].x1 > x1_1 / 2.0 && ve1.V[i].adjacentNodes.size() < 2) {
            int j;
            for (int k = 0; k < mEdges; k++) 
                if (i == ve1.E[k].v1 || i == ve1.E[k].v2) {
                    j = k; 
                    break;
                }
                ve1.addEdge(i,i+mNodes,ve1.E[j].A);
        }
    }
    // Set Neumann boundary conditions 
    count = 0;
    for (int i = 0; i < ve1.V.size(); i++) {
        // Only for nodes with one neighbour and not the entry node
        if (ve1.V[i].x1 >= ee.x1_N / 2.0 && ve1.V[i].adjacentNodes.size() < 2) {
            ve1.V[i].boundary_condition = Neumann;
            // Index of final node
            ve1.FinalNode[count] = i;
            count++;
        }
    }
		
	// Defect: Modification of edge cross-sections A based on a defect factor d
    // Center point of the defect area - circle with r radius
	double x1_D = (ee.x1_N + ee.x1_0) / 2.0;
	double x2_D = (ee.x2_N + ee.x2_0) / 2.0;
	double x3_D = (ee.x3_N + ee.x3_0) / 2.0;
	double r = 2;   
	// 
	mNodes = ve1.V.size();
	mEdges = ve1.E.size();
	for (int e = 0; e < mEdges; e++) {
		int v1 = ve1.E[e].v1;
		int v2 = ve1.E[e].v2;
		if(fabs((ve1.V[v1].x1 + ve1.V[v2].x1)/2.0-x1_D) < r && fabs((ve1.V[v1].x2 + ve1.V[v2].x2)/2.0-x2_D) < r && fabs((ve1.V[v1].x3 + ve1.V[v2].x3)/2.0-x3_D) < r) {
			ve1.E[e].A *= d; 
		}
	}	
	// Rebuild adjacency matrix
	ve1.buildAdjacent();	
	// Solve pressure
	ve1.solvePressure();
}


// Function for discretization with N divisions
void graph2(int N) {
    // Add nodes from ve1 to ve2
	for (int i = 0; i < ve1.V.size(); i++) {
		ve2.addNode(ve1.V[i].x1, ve1.V[i].x2, ve1.V[i].x3, ve1.V[i].pressure, 0, ve1.V[i].boundary_condition);
	}   
    int nodeIndex = 0;
    int nodeIndex2 = 0;
    int numNodes = ve1.V.size();
    int numEdges = ve1.E.size();
		
    // Adding discretization nodes
    for (int e = 0; e < numEdges; e++) {
        nodeIndex = ve1.E[e].v1;
        nodeIndex2 = ve1.E[e].v2;

        for (int k = 1; k < N - 1; k++) { 
            ve2.addNode(ve1.V[nodeIndex].x1 + k * (ve1.V[nodeIndex2].x1 - ve1.V[nodeIndex].x1) / (N - 1), ve1.V[nodeIndex].x2 + k * (ve1.V[nodeIndex2].x2 - ve1.V[nodeIndex].x2) / (N - 1), ve1.V[nodeIndex].x3 + k * (ve1.V[nodeIndex2].x3 - ve1.V[nodeIndex].x3) / (N - 1), ve1.V[nodeIndex].pressure + k * (ve1.V[nodeIndex2].pressure - ve1.V[nodeIndex].pressure) / (N - 1));
		}		
    }
        
    // Adding discretization edges
    for (int e = 0; e < numEdges; e++) {
        nodeIndex = ve1.E[e].v1;
        nodeIndex2 = ve1.E[e].v2;
            
        ve2.addEdge(nodeIndex, numNodes + (N - 2) * e, ve1.E[e].A, ve1.E[e].L / (N - 1), ve1.E[e].u);
            
        for (int k = 1; k < N - 2; k++) {
            ve2.addEdge(numNodes - 1 + (N - 2) * e + k , numNodes - 1 + e * (N - 2) + k + 1, ve1.E[e].A, ve1.E[e].L / (N - 1), ve1.E[e].u);
        }
		
		ve2.addEdge(numNodes - 1 + (N - 2) * (e + 1), nodeIndex2, ve1.E[e].A, ve1.E[e].L / (N - 1), ve1.E[e].u);			
    }
}


// Function to write data to a file
void write_to_file_M(string filename, double sum, double t, int N, double k, double defect)
{
    // Forming the filename with parameters
    ostringstream name;
    name << "RESULTS/" << filename.c_str() << '_'  << k << "_" << N << "_" << defect;

    // Opening the output file
    ofstream output;
    if (t == 0)
    {
        output.open(name.str());
    }
    else
    {
        output.open(name.str(), ios_base::app); // Append to existing file if t != 0
    }

    // Checking if file opened successfully
    if (!output)
    {
        cout << "Unable to open file" << endl;
    }

    // Writing time and the calculated value to the file
    output << t <<"\t" << sum << endl;

    // Closing the file
    output.close();

    return;
}

// Setting up the vascular environment
void setup_VE(int N, double dt, double D_c, double d) {
    // Creating the graph and discretizing it
    graph1(d, 2);    // Corresponds to structure ve1
    graph2(N);             // Corresponds to structure ve2 - discretized ve1
    // Initializing the concentration values and constructing the concentration matrix
    ve2.initialConditions(); 
    ve2.buildConcentrationMatrix();
    // Filling the concentration constant matrix
    ve2.fillConcentrationMatrix(D_c, dt);
}

// Calculating mass in the pipeline
double calculate_mass_VE() {
    double mass_c = 0;

    // Calculating mass contribution of each edge and summing up
    for (int e = 0; e < ve2.E.size(); e++) {
        int node1 = ve2.E[e].v1;
        int node2 =  ve2.E[e].v2;
        mass_c += ((ve2.V[node1].concentration + ve2.V[node2].concentration) * ve2.E[e].L * ve2.E[e].A) / 2.0;
    }
    return mass_c;
}

double calculate_M_bil() {
    double M_in = 0;
    double M_out = 0;
    double M_bil = 0;

    for (int index = 0; index < 9; index++) {
        M_in += ve2.V[ve1.EntryNode[index]].concentration * ve2.E[0].A * ve2.E[0].u;
        M_out += ve2.V[ve1.FinalNode[ve1.FinalNode[index]]].concentration * ve2.E[0].A * ve2.E[0].u;
    }
    M_bil = M_in - M_out;
    return M_bil;
}


void inicialization(int N, double k) {
    Nx1_pom = N;
    Nx2_pom = N;
    Nx3_pom = N;
    t_0 = 0;
    t_N = 100;
    dt = 0.01;
    diffusion_c = 0.0025;
    
    ee.dt = dt;
    ee.Nx1 = Nx1_pom;
    ee.Nx2 = Nx2_pom;
    ee.Nx3 = Nx3_pom;
    ee.dx1 = (ee.x1_N - ee.x1_0) / (ee.Nx1 - 1);
    ee.dx2 = (ee.x2_N - ee.x2_0) / (ee.Nx2 - 1);
    ee.dx3 = (ee.x3_N - ee.x3_0) / (ee.Nx3 - 1);
    ie.k = k;
    

    return;
}

// Main function
int main(int argc, char **argv) {
    // Checking if the required number of parameters is provided
    const int params = 1;
    if (argc <= params)
    {
        printf("error: required %d parameters:\n %s res[1,...]\n", params, argv[0]);
        return 1;
    }
    // Retrieving parameters from command line arguments
    int N = atoi(argv[1]);          // Number of discretization points
    double k = atof(argv[2]);       // Transfer coefficient
    double d = atof(argv[3]);       // Defect coefficient
    // Setting up space steps, diffusion, number of nodes in each axis
    inicialization(N, k);
	// Setting up the vascular environment and setting up constant matrix for concentration
	setup_VE(N, dt, diffusion_c, d);
    // Setting up the extravascular environment and setting up constant matrix for concentration
	ee.solve( );
	// Allocating memory for source terms g_c and g_m
    ie.allocate_g();
    // Convolution approximation - choose of type psi_4b
    convolution.n = 4;
    // Initial mass calculation
    double Mc = 0;
    double Mm = 0;
    double M = 0;
    double M_bil = 0;
    Mc = calculate_mass_VE();
    Mm = ee.calculate_mass_EE(ee.c_m);
    M = Mc + Mm;
	// Initial time t_0 plus time step dt
    double time = t_0 + dt;
	// Setting up neighbours of each vertex	in both environments
    ie.neighbours_in_VE();
    ie.neighbours_in_EE();	
    // Time loop for concentration calculation
    while (time <= t_N + dt) {
        // Calculating source terms
        ie.function_g_c(time);
        ie.function_g_m(time);
        // Solving concentration
        ve2.solveConcentration(dt, time, ie.g_c);
        ee.solve_concentration_m(ie.g_m);
		// Calculating mass quantities
        Mc = calculate_mass_VE();
        Mm = ee.calculate_mass_EE(ee.c_m);
        M= Mc + Mm;
        M_bil = calculate_M_bil();
        // Setting up new time
        time = time + dt;
	}    
   return 0;     
}
