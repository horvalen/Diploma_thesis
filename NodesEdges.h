// Structure for nodes
struct Node {
    double x1 = 0;
    double x2 = 0;
    double x3 = 0;

    std::vector<int> adjacentNodes;  // Indices of adjacent vertices
    std::vector<int> adjacentEdges;     // Indices of adjacent edges
    std::vector<int> adjacent_EE_Nodes;  // Indices of adjacent nodes in the extravascular environment(EE)

    double pressure = 0;          // Pressure at the node
    double concentration = 0;     // Concentration at the node
    int boundary_condition = 0;   // Boundary condition
    double volume = 0;            // Volume of the node

    // Constructor for Node struct
    Node(double x, double y, double z, double pressure1 = 0, double concentration1 = 0, int boundary_condition1 = 0, double volume1 = 0) {
        x1 = x;
        x2 = y;
        x3 = z;
        pressure = pressure1;
        concentration = concentration1;
        boundary_condition = boundary_condition1;
        volume = volume1;
    }
};

// Structure for edges
struct Edge {
    int v1, v2;               // Indices of vertices connected by the edge
    double A;     // Cross-sectional area of the edge
    double L;            // Length of the edge
    double u;          // Velocity property of the edge

    // Constructor for Edge struct
    Edge(int u, int v, double cross_section1, double length1, double velocity1) {
        v1 = u;
        v2 = v;
        A = cross_section1;
        L = length1;
        u = velocity1;
    }
};
