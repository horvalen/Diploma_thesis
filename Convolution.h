# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>
# include <string>

using namespace std;

// Global constant for pi
double pi = 3.14;

// Structure for convolution functions
struct Convolution {
    int n; // Parameter for the convolution function

    // Function for Dirac's delta approximation
	double diracDelta(double z, double d) {
    	return (1.0 / d) * approximationFunction(n, z / d);
	}

// Function for the approximation of the convolution function
    double approximationFunction(int n, double r) {
        switch (n) {
            case 1:
                if (r > -1 && r < 1) {
                    return 1 - fabs(r);
                } else {
                    return 0;
                }
                break;

            case 2:
                if (r > -0.5 && r < 0.5) {
                    return (1.0 / 3.0) * (1 + sqrt(1 - 3 * r * r));
                } else if ((r > -1.5 && r <= -0.5) || (r >= 0.5 && r < 1.5)) {
                    return (1.0 / 6.0) * (5 - 3 * fabs(r) - sqrt(-2 + 6 * fabs(r) - 3 * r * r));
                } else {
                    return 0;
                }
                break;

            case 3:
                if (r > -1 && r < 1) {
                    return (1.0 / 8.0) * (3 - 2 * fabs(r) + sqrt(1 + 4 * fabs(r) - 4 * r * r));
                } else if ((r > -2 && r <= -1) || (r >= 1 && r < 2)) {
                    return (1.0 / 8.0) * (5 - 2 * fabs(r) - sqrt(-7 + 12 * fabs(r) - 4 * r * r));
                } else {
                    return 0;
                }
                break;

            case 4:
                if (r > -2 && r < 2) {
                    return (1.0 / 4.0) * (1 + cos(pi * r / 2.0));
                } else {
                    return 0;
                }
        }
        return 0; // Default return statement
    }
};