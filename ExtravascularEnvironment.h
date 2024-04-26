#include <vector>
#include <iostream>
#include <cmath> 
#include <string>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <iomanip>

using namespace std;

// Structure for each vertex
struct vertex {
    std::vector<int> adjacent_VE_Nodes; // Neighbors in the concentration field
    double weight;
    double vertex_x1 = 0; // x-coordinate of the vertex
    double vertex_x2 = 0; // y-coordinate of the vertex
    double vertex_x3 = 0; // z-coordinate of the vertex
    double vertex_c = 0; // Concentration at the vertex
};
    
extern "C" {
#include "/usr/include/suitesparse/umfpack.h"
}

// EE = Extravascular Environment
struct EE {
	SparseMatrix<double> EE_SparseMatrix; // Matrix for calculating concentration in the extravascular environment
    // Parameters for the external environment
    int Nx1; // Number of points in x-direction
    int Nx2; // Number of points in y-direction
    int Nx3; // Number of points in z-direction
    
    double dt; // Time step
    
    // Domain limits
    double x1_0 = 0;
    double x1_N = 10;
    double dx1;
    
    double x2_0 = 0;
    double x2_N = 10;
    double dx2;
    
    double x3_0 = 0;
    double x3_N = 10;
    double dx3;
    
    double diffusion_m = 0.003; // Diffusion of the medium

    double *x1 = 0; // Array for x-coordinates
    double *x2 = 0; // Array for y-coordinates
    double *x3 = 0; // Array for z-coordinates
    double *c_m = 0; // Array for concentrations
    
    vertex *vertices = 0; // Array of vertices

    // Destructor to deallocate memory
    ~EE() {
        if (x1 != 0)
            delete[] x1;
        if (x2 != 0)
            delete[] x2;
        if (x3 != 0)
            delete[] x3;
        if (c_m != 0)
            delete[] c_m;
        if (vertices != 0)
            delete[] vertices;
    }

    // Function to calculate index in 3D array
    int index(int i, int j, int k) {
        // Handle periodic boundary conditions
        if (i < 0)
            return index(Nx1 - 1, j, k);
        if (j < 0)
            return index(i, Nx2 - 1, k);
        if (k < 0)
            return index(i, j, Nx3 - 1);
        if (i >= Nx1)
            return index(0, j, k);
        if (j >= Nx2)
            return index(i, 0, k);
        if (k >= Nx3)
            return index(i, j, 0);
        return i + Nx1 * (j + Nx2 * k);
    }

    // Function to allocate memory
    void allocate() {
        x1 = new double[Nx1];
        x2 = new double[Nx2];
        x3 = new double[Nx3];
        c_m = new double[Nx1 * Nx2 * Nx3];
        vertices = new vertex[Nx1 * Nx2 * Nx3];
    }

    // Function to initialize spatial coordinates
    void initialize_coordinates() {
        for (int i = 0; i < Nx1; i++) {
            x1[i] = ((Nx1 - 1 - i) * x1_0 + i * x1_N) / (Nx1 - 1);
        }
        for (int j = 0; j < Nx2; j++) {
            x2[j] = ((Nx2 - 1 - j) * x2_0 + j * x2_N) / (Nx2 - 1);
        }
        for (int k = 0; k < Nx3; k++) {
            x3[k] = ((Nx3 - 1 - k) * x3_0 + k * x3_N) / (Nx3 - 1);
        }
    }
	
    // Function to write data to file
    void write_to_file(string filename, int timestep) {
        ostringstream name_stream, dir_stream;
        
        // Create directory for results
        dir_stream << "RESULTS/Paraview_external";
        mkdir(dir_stream.str().c_str(), 0755);
        
        // Construct filename
        name_stream << dir_stream.str().c_str() << "/" << filename.c_str() << "_" << timestep << ".vtk";
        
        // Open file for writing
        ofstream output;
        output.open(name_stream.str());
        
        // Check if file opened successfully
        if (!output) {
            cout << "Unable to open file" << endl;
        }
        
        // Write header information
        output << "# vtk DataFile Version 3.1" << endl;
        output << "Author: Lenka Horvátová" << endl;
        output << "ASCII" << endl;
        output << "DATASET STRUCTURED_GRID" << endl << endl;
        output << "DIMENSIONS " << Nx1 << " " << Nx2 << " " << Nx3 << endl << endl;
        
        // Write grid points
        output << "POINTS " << Nx1 * Nx2 * Nx3 << " FLOAT" << endl;
        for (int k = 0; k < Nx3; k++) {
            for (int j = 0; j < Nx2; j++) {
                for (int i = 0; i < Nx1; i++) {
                    output << x1[i] << " " << x2[j] << " " << x3[k] << endl;
                }
            }
        }
        output << endl;

        // Write concentration values
        output << "POINT_DATA " << Nx1 * Nx2 * Nx3 << endl;
        output << "SCALARS concentration FLOAT" << endl;
        output << "LOOKUP_TABLE default" << endl;
        for (int r = 0; r < EE_SparseMatrix.N; r++) {
            output << c_m[r] << endl;
        }
        output << endl;

        // Close the file
        output.close();
    }
		
    // Function to set initial conditions
    void set_initial_conditions() {
        for (int i = 0; i < Nx1; i++) {
            for (int j = 0; j < Nx2; j++) {
                for (int k = 0; k < Nx3; k++) {
                    c_m[index(i, j, k)] = 0;
                    vertices[index(i, j, k)].vertex_c = c_m[index(i, j, k)];
                }
            }
        }
    }
	
	
	// Function to calculate the integral of the concentration field in the external environment
	double calculate_mass_EE(double c_m[]) {
		double M = 0;
		// Calculate contributions from internal points
		for (int i = 0; i < Nx1 ; i++) {
			for (int j = 0; j < Nx2 ; j++) {
				for (int k = 0; k < Nx3 ; k++) {
					M += dx1 * dx2 * dx3 * c_m[index(i, j, k)];
				}
			}
		}

		return M;
	}	
		
	// Function to assemble the matrix for the external environment
	void assemble_matrix_m() {
		EE_SparseMatrix.N = Nx1 * Nx2 * Nx3; // Total number of grid points
		EE_SparseMatrix.NCC = 7 * EE_SparseMatrix.N; // Number of non-zero elements (7 neighbors in 3D)

		// Allocate memory for the matrix
		EE_SparseMatrix.allocate();

		int p = 0; // Counter for the current position in the matrix arrays
		int boundary_helper = 0; // Helper variable for boundary conditions

		// Loop through all grid points
		for (int k = 0; k < Nx3; k++) {
			for (int j = 0; j < Nx2; j++) {
				for (int i = 0; i < Nx1; i++) {
					int r = index(i, j, k); // Current grid point index
					EE_SparseMatrix.Ax[p] = 0; // Set the corresponding value in the matrix to zero
					EE_SparseMatrix.Ai[p] = r; // Store the diagonal element index
					p++; // Move to the next position in the matrix arrays

					// Add neighbors to the matrix
					// Right neighbor
					if (i == Nx1 - 1) {
						boundary_helper = -1;
						EE_SparseMatrix.Ai[p] = index(boundary_helper + 1, j, k);
					} else {
						EE_SparseMatrix.Ai[p] = index(i + 1, j, k);
					}
					EE_SparseMatrix.Ax[p] = 0; // Set the corresponding value in the matrix to zero
					p++; // Move to the next position in the matrix arrays

					// Bottom neighbor
					if (j == Nx2 - 1) {
						boundary_helper = -1;
						EE_SparseMatrix.Ai[p] = index(i, boundary_helper + 1, k);
					} else {
						EE_SparseMatrix.Ai[p] = index(i, j + 1, k);
					}
					EE_SparseMatrix.Ax[p] = 0; // Set the corresponding value in the matrix to zero
					p++; // Move to the next position in the matrix arrays

					// Front neighbor
					if (k == Nx3 - 1) {
						boundary_helper = -1;
						EE_SparseMatrix.Ai[p] = index(i, j, boundary_helper + 1);
					} else {
						EE_SparseMatrix.Ai[p] = index(i, j, k + 1);
					}
					EE_SparseMatrix.Ax[p] = 0; // Set the corresponding value in the matrix to zero
					p++; // Move to the next position in the matrix arrays

					// Left neighbor
					if (i == 0) {
						boundary_helper = Nx1;
						EE_SparseMatrix.Ai[p] = index(boundary_helper - 1, j, k);
					} else {
						EE_SparseMatrix.Ai[p] = index(i - 1, j, k);
					}
					EE_SparseMatrix.Ax[p] = 0; // Set the corresponding value in the matrix to zero
					p++; // Move to the next position in the matrix arrays

					// Top neighbor
					if (j == 0) {
						boundary_helper = Nx2;
						EE_SparseMatrix.Ai[p] = index(i, boundary_helper - 1, k);
					} else {
						EE_SparseMatrix.Ai[p] = index(i, j - 1, k);
					}
					EE_SparseMatrix.Ax[p] = 0; // Set the corresponding value in the matrix to zero
					p++; // Move to the next position in the matrix arrays

					// Back neighbor
					if (k == 0) {
						boundary_helper = Nx3;
						EE_SparseMatrix.Ai[p] = index(i, j, boundary_helper - 1);
					} else {
						EE_SparseMatrix.Ai[p] = index(i, j, k - 1);
					}
					EE_SparseMatrix.Ax[p] = 0; // Set the corresponding value in the matrix to zero
					p++; // Move to the next position in the matrix arrays
				}
			}
		}

		// Set the starting positions for each row in the matrix
		EE_SparseMatrix.Ap[0] = 0;
		for (int c = 0; c < EE_SparseMatrix.N; c++) {
			EE_SparseMatrix.Ap[c + 1] = EE_SparseMatrix.Ap[c] + 7;
		}

		// Sort the elements in each row of the matrix
		int num_parts = EE_SparseMatrix.NCC / 7.0;
		for (int part = 0; part < num_parts; part++) {
			for (int s = 0 + 7 * part; s < 6 + 7 * part; s++) {
				int min_index = s;
				for (int t = s + 1; t < 7 + 7 * part; t++) {
					if (EE_SparseMatrix.Ai[t] < EE_SparseMatrix.Ai[min_index]) {
						min_index = t;
					}
				}

				if (min_index != s) {
					swap(EE_SparseMatrix.Ai[s], EE_SparseMatrix.Ai[min_index]);
				}
			}
		}
	}
            	
	// Function to fill the matrix M
	void fill_matrix_M() {
		// Loop through all grid points
		for (int k = 0; k < Nx3; k++) {
			for (int j = 0; j < Nx2; j++) {
				for (int i = 0; i < Nx1; i++) {
					int r = index(i, j, k); // Current grid point index

					// Fill the matrix entries with coefficients for the diffusion equation
					EE_SparseMatrix(r, index(i, j, k)) = 1 + 2 * dt * diffusion_m / (dx1 * dx1) + 2 * dt * diffusion_m / (dx2 * dx2) + 2 * dt * diffusion_m / (dx3 * dx3);
					EE_SparseMatrix(r, index(i - 1, j, k)) = -dt * diffusion_m / (dx1 * dx1);
					EE_SparseMatrix(r, index(i + 1, j, k)) = -dt * diffusion_m / (dx1 * dx1);
					EE_SparseMatrix(r, index(i, j - 1, k)) = -dt * diffusion_m / (dx2 * dx2);
					EE_SparseMatrix(r, index(i, j + 1, k)) = -dt * diffusion_m / (dx2 * dx2);
					EE_SparseMatrix(r, index(i, j, k - 1)) = -dt * diffusion_m / (dx3 * dx3);
					EE_SparseMatrix(r, index(i, j, k + 1)) = -dt * diffusion_m / (dx3 * dx3);
				}
			}
		}
	}

	// Function to solve for concentration in the exterior environment
	void solve_concentration_m(double gm[]) {
		// Allocate memory for the vectors b_m and x_m
		double *b_m = new double[EE_SparseMatrix.N];
		double *x_m = new double[EE_SparseMatrix.N];

		// Initialize vectors b_m and x_m to zeros
		for (int i = 0; i < EE_SparseMatrix.N; i++) {
			x_m[i] = 0; 
			b_m[i] = 0;
		}

		// Update the right-hand side vector b_m using the current concentration and source term
		for (int k = 0; k < Nx3; k++) {
			for (int j = 0; j < Nx2; j++) {
				for (int i = 0; i < Nx1; i++) {
					int r = index(i, j, k); // Calculate the 1D index from 3D (i, j, k)
					b_m[r] = c_m[r] + dt * gm[r]; // Update b_m with concentration and source term
				}
			}
		}

		// Solve the linear system EE_SparseMatrix * x_m = b_m
		EE_SparseMatrix.solve(b_m, x_m);

		// Update concentration values and related data using the solution x_m
		for (int k = 0; k < Nx3; k++) {
			for (int j = 0; j < Nx2; j++) {
				for (int i = 0; i < Nx1; i++) {
					int r = index(i, j, k); // Calculate the 1D index from 3D (i, j, k)
					c_m[r] = x_m[r]; // Update concentration values
					vertices[r].vertex_c = c_m[r]; // Update associated vertex concentration
				}
			}
		}

		// Deallocate memory for b_m and x_m
		delete [] b_m;
		delete [] x_m;
	}
	
	// Function to resolve the problem
	void solve() {
		// Allocate memory for necessary arrays
		allocate();
		// Initial condition
		set_initial_conditions();
		// Generate grid points
		initialize_coordinates();
		// Assemble the matrix for the exterior environment
		assemble_matrix_m();
		// Constructing the concentration matrix
		fill_matrix_M();

		// Update the coordinates of vertices with grid points
		for (int k = 0; k < Nx3; k++) {
			for (int j = 0; j < Nx2; j++) {
				for (int i = 0; i < Nx1; i++) {
					int r = index(i, j, k); // Calculate the 1D index from 3D (i, j, k)

					// Update the coordinates of the vertex with grid points
					vertices[r].vertex_x1 = x1[i];
					vertices[r].vertex_x2 = x2[j];
					vertices[r].vertex_x3 = x3[k];
				}
			}
		}
	}
};
