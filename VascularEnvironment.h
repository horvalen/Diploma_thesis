#include <vector>
#include <iostream>
#include <cmath> 
#include <string>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>

#include "NodesEdges.h"             // Include header for node and edge structures
#include "SparseMatrixUMFPACK.h"    // Include header for sparse matrix implementation
#include "Convolution.h"            // Include header for convolution operations

using namespace std;

extern "C" {
#include "/usr/include/suitesparse/umfpack.h"    // External C library for sparse matrix operations
}

// Define aliases for vector of nodes and vector of edges
using NodeV = std::vector< Node >;
using EdgeE = std::vector< Edge >;

// Numbering of boundary conditions
int Dirichlet = 1;    // Constant representing Dirichlet boundary condition
int Neumann = 2;      // Constant representing Neumann boundary condition

// VE = Vascular Environment
struct VE{
    NodeV V;                         // Vector of nodes
    EdgeE E;                         // Vector of edges
    SparseMatrix<int> A;             // Sparse matrix for system
    SparseMatrix<double> VE_SparseMatrix;         // Sparse matrix for vascular environment
    
    // Optional constants
    double K_source = 0.5;            // Source constant
    double t_source = 5;              // Time source constant
 
    int EntryNode[9];                // Array for entry nodes
    int HelpNode[9];                 // Array for help nodes
    int FinalNode[9];                // Array for final nodes

    // Function to add a node
    int addNode(double x, double y, double z, double p = 0, double c = 0, int bc = 0, double Vv = 0) {
        Node node(x, y, z);          // Create a new node
        node.pressure = p;           // Assign pressure
        node.concentration = c;      // Assign concentration
        node.boundary_condition = bc; // Assign boundary condition
        node.volume = Vv;            // Assign volume
            
        V.push_back(node);           // Add the node to the vector
        return V.size() - 1;         // Return the index of the added node
    }
    
    // Function to add an edge
    int addEdge(int v1, int v2, double A = 1.0, double length1 = -1, double velocity = 0) {
        // Calculate the length between nodes
        double length = sqrt((V[v1].x1 - V[v2].x1) * (V[v1].x1 - V[v2].x1) + 
                             (V[v1].x2 - V[v2].x2) * (V[v1].x2 - V[v2].x2) + 
                             (V[v1].x3 - V[v2].x3) * (V[v1].x3 - V[v2].x3));
        
        Edge edge(v1, v2, A, length, velocity); // Create a new edge
        edge.u = velocity;                           // Assign velocity
        
        E.push_back(edge);  // Add the edge to the vector
        
        // Connect nodes with the edge
        V[v1].adjacentNodes.push_back(v2);
        V[v2].adjacentNodes.push_back(v1);
        
        V[v1].adjacentEdges.push_back(E.size() - 1);
        V[v2].adjacentEdges.push_back(E.size() - 1);
        
        return E.size() - 1;  // Return the index of the added edge
    }

    // Function to print nodes
    void printNodes() {
        std::cout << "nodes: " << std::endl;
        for (int i = 0; i < V.size(); i++) {
            std::cout << i << ": " << V[i].x1 << "," << V[i].x2 << "," << V[i].x3 << "  pressure: " << V[i].pressure;
            cout << endl;
        }
    }
    
    // Function to print edges
    void printEdges() {
        std::cout << "edges: " << std::endl;
        for (int i = 0; i < E.size(); i++) {
            std::cout << i << ": " << E[i].v1 << "--" << E[i].v2 << "  velocity: " << E[i].u << std::endl;
        }
    }
    
    // Function to write data to a file
    void writeToFile(string filename, int p, double t) {
        ostringstream nazev, dir1;
        
        dir1 << "VYSLEDKY/Paraview_cevni";
        mkdir(dir1.str().c_str(), 0755);
        
        nazev << dir1.str().c_str() << "/" << filename.c_str() << "_" << p << ".vtk";
        
        ofstream output;
        output.open(nazev.str());
        
        if (!output) {
            cout << "Nelze otevrit soubor" << endl;
        }
        
        output << "# vtk DataFile Version 3.1" << endl;
        output << "Autor: Lenka Horvátová" << endl;
        output << "ASCII" << endl;
        output << "DATASET UNSTRUCTURED_GRID" << endl << endl;
        
        output << "POINTS " << V.size() << " FLOAT" << endl;
        
        for (int i = 0; i < V.size(); i++) {
            output << V[i].x1 << " " << V[i].x2 << " " << V[i].x3 << endl;
        }
        output << endl;
        
        output << "CELLS " << E.size() << " " << 3 * E.size() << endl;
        
        for (int i = 0; i < E.size(); i++) {
            output << 2 << " " << E[i].v1 << " " << E[i].v2 << endl;
        }        
        output << endl;
        
        output << "CELL_TYPES " << E.size() << endl;
        
        for (int i = 0; i < E.size(); i++) {
            output << 3 << " ";
        }
        output << endl << endl;
        
        output << "POINT_DATA " << V.size() << endl;
        output << "SCALARS pressure FLOAT" << endl;
        output << "LOOKUP_TABLE default" << endl;
        
        for (int i = 0; i < V.size(); i++) {
            output << V[i].pressure << endl;
        }
        output << endl;
        
        output << "CELL_DATA " << E.size() << endl;
        output << "SCALARS velocity FLOAT" << endl;
        output << "LOOKUP_TABLE default" << endl;
        
        for (int i = 0; i < E.size(); i++) {
            output << fabs(E[i].u) << endl;
        }
        output << endl;
        
        output << "POINT_DATA " << V.size() << endl;
        output << "SCALARS concentration FLOAT" << endl;
        output << "LOOKUP_TABLE default" << endl;
        
        for (int i = 0; i < V.size(); i++) {
            output << V[i].concentration << endl;
        }
        
        output.close();
    }

    // Function to build the adjacency matrix
 
   
//...........VYPOCET ADJACENCNI MATICE......................... 
    void buildAdjacent() {
        A.N = V.size();
        A.NCC = E.size() * 2;
        A.allocate();
        
        int p = 0;
        A.Ap[0] = 0;
        
        for (int c = 0; c < A.N; c++) {
            int col = 0;
            for (int r = 0; r < A.N; r++) {
                if (c != r) {
                    for (int i = 0; i < E.size(); i++) {
                        if ((E[i].v1 == r && E[i].v2 == c) || (E[i].v1 == c && E[i].v2 == r)) {
                            A.Ax[p] = i + 1;
                            A.Ai[p] = r;
                            p++;
                            col++; 
                        }
                    }    
                }
                A.Ap[c + 1] = A.Ap[c] + col;
            }
        }
    }
    
    // Function to calculate velocity
    double velocity(int e) {
        int v1 = E[e].v1;
        int v2 = E[e].v2;
        return -1 / (E[e].L) * (V[v2].pressure - V[v1].pressure);
    }

    // Function to calculate velocity with specified nodes
    double velocity(int v1, int v2, int e) {
        return -1 / (E[e].L) * (V[v2].pressure - V[v1].pressure);
    }
    

	





    // Function to solve pressure using sparse matrix
    void solvePressure() {
        SparseMatrix<double> P;  // Matrix for computing pressure
        
        P.N = V.size();
        
        int NCC = 0;
        
        // Determine the number of non-zero elements
        for (int r = 0; r < P.N; r++) {
            for (int c = 0; c < P.N; c++) {
                if (A(r, c) != 0 && V[r].boundary_condition == 0)
                    NCC++;
            }
        }
        
        P.NCC = P.N + NCC;
        P.allocate();
        
        int p = 0;
        P.Ap[0] = 0;
        
        for (int c = 0; c < P.N; c++) {
            int col = 0;
            for (int r = 0; r < P.N; r++) {
                if (c == r) {
                    P.Ax[p] = 0;
                    P.Ai[p] = r;
                    p++;
                    col++;
                } else {
                    if (A(r, c) != 0 && V[r].boundary_condition == 0) {
                        P.Ax[p] = 0;
                        P.Ai[p] = r;
                        p++;
                        col++;
                    }
                }
                P.Ap[c + 1] = P.Ap[c] + col;
            }
        }
        
        double *b = new double[P.N];
        double *x = new double[P.N];
        
        for (int i = 0; i < P.N; i++) { 
            x[i] = 0;   
            
            if (V[i].boundary_condition == 0) {
                b[i] = 0;
            } else if (V[i].boundary_condition == 1) {
                b[i] = 50;
            } else if (V[i].boundary_condition == 2) {
                b[i] = 1;
            }
        }
        
        double edgeP = 0; 
        for (int r = 0; r < P.N; r++) {
            P(r, r) = 1.0;
            
            if (V[r].adjacentEdges.size() > 1) {
                P(r, r) = 0;
                for (int s = 0; s < V[r].adjacentEdges.size(); s++) {
                    edgeP = V[r].adjacentEdges[s];
                    
                    P(r, r) +=  E[edgeP].A / E[edgeP].L;
                    P(r, V[r].adjacentNodes[s]) = -E[edgeP].A / E[edgeP].L;
                }
            }
        }
        
        P.solve(b, x);
        
        for (int i = 0; i < P.N; i++) {
            V[i].pressure = x[i]; 
        }
        
        for (int i = 0; i < E.size(); i++) {
            E[i].u = velocity(i); 
        }
        
        delete [] x;
        delete [] b;
    }

    // Function to build the matrix for concentration
    void buildConcentrationMatrix() {
        VE_SparseMatrix.N = V.size();
        
        int NCCC = 0;
        
        for (int i = 0; i < V.size(); i++) {
            if (V[i].adjacentNodes.size() < 2) {
                if (V[i].boundary_condition == Dirichlet)
                    NCCC += 1;
                else
                    NCCC += 2;
            } else {
                NCCC += V[i].adjacentNodes.size() + 1;
            }
        }

        VE_SparseMatrix.NCC = NCCC;
            
        VE_SparseMatrix.allocate();
        
        int p = 0;
        VE_SparseMatrix.Ap[0] = 0;
        
        for (int c = 0; c < VE_SparseMatrix.N; c++) {
            int col = 0;
            
            for (int r = 0; r < VE_SparseMatrix.N; r++) {
                if (c == r) {
                    VE_SparseMatrix.Ax[p] = 0;
                    VE_SparseMatrix.Ai[p] = r;
                    p++;
                    col++;
                } else {
                    for (int s = 0; s < V[c].adjacentNodes.size(); s++) {
                        if (r == V[c].adjacentNodes[s] && V[r].adjacentNodes.size() >= 2) {
                            VE_SparseMatrix.Ax[p] = 0;
                            VE_SparseMatrix.Ai[p] = r;
                            p++;
                            col++; 
                        } else if (r == V[c].adjacentNodes[s] && V[r].adjacentNodes.size() == 1 && (V[V[c].adjacentNodes[s]].boundary_condition == 2)) {
                            VE_SparseMatrix.Ax[p] = 0;
                            VE_SparseMatrix.Ai[p] = r;
                            p++;
                            col++;  
                        }
                    }
                }
                
                VE_SparseMatrix.Ap[c + 1] = VE_SparseMatrix.Ap[c] + col;
            }
        }
    }

    // Function to fill the concentration matrix
    void fillConcentrationMatrix(double diffusion, double dt) {
        int edge = 0;
        // Loop through each node
        for (int r = 0; r < VE_SparseMatrix.N; r++) {
            // Set diagonal matrix to 1
            VE_SparseMatrix(r, r) = 1.0;
            // If node has more neighbours
            if (V[r].adjacentEdges.size() > 1) {	
                // Calculate the total volume of adjacent edges
                for (int s = 0; s < V[r].adjacentEdges.size(); s++) {
                    edge = V[r].adjacentEdges[s];
                    V[r].volume += E[edge].A * E[edge].L / 2.0;
                }
                // Loop through the neighbors of each node
                for (int s = 0; s < V[r].adjacentEdges.size(); s++) {
                    edge = V[r].adjacentEdges[s];
                    // Velocity direction for upwind scheme
                    double vel = velocity(r, V[r].adjacentNodes[s], edge); 
                    // Diffusion part
                    VE_SparseMatrix(r, r) += E[edge].A * dt * diffusion / (E[edge].L * V[r].volume); 
                    VE_SparseMatrix(r, V[r].adjacentNodes[s]) += -E[edge].A * diffusion * dt / (E[edge].L * V[r].volume);
                    // Advection part
                    if (vel < 0) {
                        VE_SparseMatrix(r, r) += -E[edge].A * vel * dt / V[r].volume; 
                        VE_SparseMatrix(r, V[r].adjacentNodes[s]) +=  E[edge].A * vel * dt / V[r].volume;            
                    }
                }
            }
            else {
                // Neumann boundary condition
                if (V[r].boundary_condition == 2) {
                    VE_SparseMatrix(r, V[r].adjacentNodes[0]) = -1;
                }
            }
        }
    }

    // Function to set initial concentration conditions
    void initialConditions() {
        for (int k = 0; k < V.size(); k++) {
            V[k].concentration = 0;
        }
    }

    // Function to compute the source term
    double source(double t) {
        return exp(-K_source * (t - t_source) * (t - t_source));
        // Alternatively, you can use a constant source term:
        // return 1;
    }
 
    // Function to solve the concentration equation
    void solveConcentration(double dt, double t, double gc[]) {
        double *b = new double[VE_SparseMatrix.N];
        double *x = new double[VE_SparseMatrix.N];

        // Initialize the vectors
        for (int i = 0; i < VE_SparseMatrix.N; i++) {
            x[i] = 0; 
            b[i] = 0;
        }

        // Populate the right-hand side vector
        for (int r = 0; r < VE_SparseMatrix.N; r++) {
            if (V[r].boundary_condition == 0) {
                b[r] = V[r].concentration + dt * gc[r];
            } else {
                if (V[r].boundary_condition == 1) {
                    b[r] = source(t);
                } else if (V[r].boundary_condition == 2) {
                    b[r] = 0;
                }
            }
        }

        // Solve the concentration equation
        VE_SparseMatrix.solve(b, x);

        // Update the concentrations
        for (int i = 0; i < V.size(); i++) {
            V[i].concentration = x[i];
        }

        // Clean up memory
        delete[] b;
        delete[] x;
    }
};


