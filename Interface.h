#include <vector>
#include <iostream>
#include <cmath> 
#include <string>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <iomanip>

#include "VascularEnvironment.h"	
#include "ExtravascularEnvironment.h"

using namespace std;

// Declaration of the Vascular and Extravascular Environment (VE) objects
VE ve1; // Non-discretized
VE ve2; // ve1 + discretized vertices and edges
EE ee; // Extravascular environment

// Convolution object
Convolution convolution;

// Structure of external sources
struct IE { 
    double *g_m = nullptr; // Source function for the extravascular environment
    double *g_c = nullptr; // Source function for the vascular environment
    
	double k = 0; // Transfer coefficient
   
    
	~IE() {
        // Destructor to deallocate memory
        if (g_m != nullptr)
            delete[] g_m;
        if (g_c != nullptr)
            delete[] g_c;
    }
	
    // Allocate memory for source functions
    void allocate_g() {
        g_m = new double[ee.Nx1 * ee.Nx2 * ee.Nx3];
        g_c = new double[ve2.V.size()];
    }
	
	// Compute source function for the extravascular environment
    void function_g_m(double t) {
		int X;
		// Iterate over all vertices in the extravascular environment
    	for ( int x = 0; x < ee.Nx1*ee.Nx2*ee.Nx3; x++ ) {
			g_m[x] = 0; // Initialize the source function value for the current vertex

			// Iterate over adjacent vascular nodes for the current extravascular vertex
            for (int s = 0; s < ee.vertices[x].adjacent_VE_Nodes.size() ; s++ ){
				X= ee.vertices[x].adjacent_VE_Nodes[s];
            				
            	g_m[x] += k*(ve2.V[X].concentration - ee.vertices[x].vertex_c)*
						convolution.diracDelta(ee.vertices[x].vertex_x1 - ve2.V[X].x1,ee.dx1)*
						convolution.diracDelta(ee.vertices[x].vertex_x2 - ve2.V[X].x2,ee.dx2)*
						convolution.diracDelta(ee.vertices[x].vertex_x3 - ve2.V[X].x3,ee.dx3)*
						ve2.V[X].volume;				
            }
    	}
    	return;
    }
	
	// Compute source function for the vascular environment
	void function_g_c(double t) {
		int x;
		// Iterate over all vertices in the vascular environment
		for (int X = 0; X < ve2.V.size(); X++) {
			g_c[X] = 0; // Initialize the source function value for the current vertex

			// Iterate over adjacent extravascular nodes for the current vascular vertex
			for (int ijk = 0; ijk < ve2.V[X].adjacent_EE_Nodes.size(); ijk++) {
				x = ve2.V[X].adjacent_EE_Nodes[ijk];
				
				g_c[X] += k * (ee.vertices[x].vertex_c - ve2.V[X].concentration) *
					convolution.diracDelta(ee.vertices[x].vertex_x1 - ve2.V[X].x1, ee.dx1) *
					convolution.diracDelta(ee.vertices[x].vertex_x2 - ve2.V[X].x2, ee.dx2) *
					convolution.diracDelta(ee.vertices[x].vertex_x3 - ve2.V[X].x2, ee.dx3) *
					ee.dx1*ee.dx2*ee.dx3;
			}
		}
		return;
	}

	// Neighbours of extravascular vertices in vascular environment
	void neighbours_in_VE() {
	    // Iterate over the discretized space of the extravascular environment
		for (int i = 0; i < ee.Nx1; i++) {
			for (int j = 0; j < ee.Nx2; j++) {
				for (int k = 0; k < ee.Nx3; k++) {
					// Iterate over all vertices in the vascular environment
					for (int l = 0; l < ve2.V.size(); l++) {
						// Check if the vascular vertex is adjacent to the current extravascular node
						if ((convolution.diracDelta(ve2.V[l].x1 - ee.x1[i], ee.dx1) > 0) &&
							(convolution.diracDelta(ve2.V[l].x2 - ee.x2[j], ee.dx2) > 0) &&
							(convolution.diracDelta(ve2.V[l].x3 - ee.x3[k], ee.dx3) > 0)) {
							// Add the index of the adjacent vascular vertex to the current extravascular node's list of adjacent vertices
							ee.vertices[ee.index(i, j, k)].adjacent_VE_Nodes.push_back(l);
						}
					}
				}
			}
		}
	}
	
	// Neighbours of vascular vertices in etravascular environment
	void neighbours_in_EE() {
		// Iterate over all vertices in the vascular environment
		for (int l = 0; l < ve2.V.size(); l++) {
			// Iterate over the discretized space of the extravascular environment
			for (int k = 0; k < ee.Nx3; k++) {
				for (int j = 0; j < ee.Nx2; j++) {
					for (int i = 0; i < ee.Nx1; i++) {
						// Check if the extravascular node is adjacent to the current vascular vertex
						if ((convolution.diracDelta(ee.x1[i] - ve2.V[l].x1, ee.dx1) > 0) &&
							(convolution.diracDelta(ee.x2[j] - ve2.V[l].x2, ee.dx2) > 0) &&
							(convolution.diracDelta(ee.x3[k] - ve2.V[l].x3, ee.dx3) > 0)) {
							// Add the index of the adjacent extravascular node to the current vascular vertex's list of adjacent vertices
							ve2.V[l].adjacent_EE_Nodes.push_back(ee.index(i, j, k));
						}
					}
				}
			}
		}
	}
};

