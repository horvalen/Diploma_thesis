extern "C" {
#include "/usr/include/suitesparse/umfpack.h"
}

// Sparse matrix structure
template<typename T>
struct SparseMatrix {
    int N;     // Size of the matrix
    int NCC;   // Non-zero count (number of non-zero entries)
    
    T node_error = 0;               // Default value for non-existing entries
    double *null = (double *)NULL;  // Null pointer for UMFpack functions
    void *Symbolic = NULL;    // Symbolic factorization structure
    void *Numeric = NULL;     // Numeric factorization structure

    int *Ai = 0;              // Array of row indices of non-zero entries
    int *Ap = 0;              // Array of column pointers indicating the start of each column's non-zero entries
    T *Ax = 0;                // Array of non-zero entries
    
    // Function to solve a system of linear equations Ax = b
    void solve(T *b, T *x) {
        int status;
        // Perform symbolic factorization if not already done
        if (Symbolic == NULL)
        	status = umfpack_di_symbolic(N, N, Ap, Ai, Ax, &Symbolic, null, null);

        // Perform numeric factorization if not already done
        if (Numeric == NULL)
        	status = umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, null, null);

        // Solve the system of linear equations
        status = umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null);
    }
    
    // Function to allocate memory for the matrix
    void allocate() {
        Ai = new int[NCC];
        Ap = new int[N + 1];
        Ax = new T[NCC];
    }
    
    // Function to get the position of a non-zero element in the internal storage
    int xpos(int i, int j) {
        for (int k = Ap[j]; k < Ap[j + 1]; k++)
            if (Ai[k] == i)
                return k;
        return -1;  // Return -1 if the element does not exist
    }

    // Function to get the value of an element at position (i, j)
	T& get(int i, int j) {
        int k = xpos(i, j);
        if (k == -1) {
            return node_error;  // Return the default value if the element does not exist
        }
        return Ax[k];
    }
	
    // Overloaded operator for accessing matrix elements
	inline T& operator()(int i, int j) {
        return get(i, j);
    }

    // Function to print the matrix
    void printMatrix() {
        printf("Matrix print \n");
        for (int i = 0; i < N; i++) {
            printf("row %d: ", i);
            for (int j = 0; j < N; j++) {
                std::cout << get(i, j) << "\t";
            }
            printf(" \n");
        }
    }
    
    // Destructor to free memory and UMFpack structures
    ~SparseMatrix() {
        if (Ai)
            delete[] Ai;
        if (Ap)
            delete[] Ap;
        if (Ax)
            delete[] Ax;
        if (Symbolic)
            umfpack_di_free_symbolic(&Symbolic);
        if (Numeric)
            umfpack_di_free_numeric(&Numeric);
    } 
};


