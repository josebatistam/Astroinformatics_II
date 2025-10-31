#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

using namespace std;

// Function to calculate the luminosity distance DL at redshift z
double calculateDL(double z){
    // Using the ΛCDM standard cosmology parameters
    const double c = 299792.458;                  // Speed of light in km/s
    const double H0 = 70.0;                       // Hubble constant in km/s/Mpc
    const double h = H0 / 100;                    // Hubble parameter
    const double OmegaM = 0.3;                    // Matter density
    const double OmegaR = 4.165E-5 / pow(h, 2);   // Radiation density
    const double OmegaV = 0.7;                    // Vacuum density
    const double OmegaK = 1 - OmegaM - OmegaR - OmegaV;   // Curvature density
    double az = 1 / (1 + z);                      // Scale factor at redshift z
    // Integrate to get the radial comoving distance (DCMR)
    double DCMR = 0.0;
    int n = 10000;                               // number of points in integral
    // Integrate DCMR = ∫ da / (a * sqrt(Ω_K + Ω_M/a + Ω_R/a² + Ω_Λ a²))
    for (int i = 1; i < n; ++i) {
        double a = az + (1.0 - az) * (i + 0.5) / n;
        double adot = sqrt(OmegaK + (OmegaM / a) + (OmegaR / pow(a,2)) + (OmegaV * pow(a,2)));
        DCMR += 1.0 / (a * adot);
    }
    DCMR = (1 - az) * DCMR / n;
    // Calculate the tangential comoving distance in the FLRW coordinate
    double DCMT = 1.0;
    if (OmegaK < 0.0) {
        double SOmegaK = sqrt(fabs(OmegaK));
        DCMT = sin(SOmegaK * DCMR) / SOmegaK;
    } else if (OmegaK > 0.0) {
        double SOmegaK = sqrt(abs(OmegaK));
        DCMT = sinh(SOmegaK * DCMR) / SOmegaK;
    } else {
        DCMT = DCMR;
    }
    double DA = (c/H0) * (az * DCMT);             // Angular distance in Mpc
    double DL = DA / pow(az,2);                   // Luminosity distance in Mpc
    return DL;
}

int main() {
    ifstream infile("coma_members_Jim2025.dat");           // Read the data file
    if (!infile.is_open()) {
        cerr << "Error: cannot open file: coma_members_Jim2025.dat" << endl;
        return 1;
    }
    // Skip the header line starting with '#'
    string header_line;                
    getline(infile, header_line);
    double ra, dec, z;                      // Variables for the three columns
    int galaxy_id = 0;                      // Initialize the galaxy id variable
    cout.setf(ios::fixed);
    cout.precision(3);                      // Set output precision
    cout << "ID\tDistance(Mpc)" << endl;    // Print the output header
    while (infile >> ra >> dec >> z) {
        galaxy_id++;                        // Use a line counter as the ID
        double DL = calculateDL(z);         // Calculate the luminosity distance
        cout << galaxy_id << "\t" << DL << endl; // Print the galaxy ID and DL
    }
    infile.close();                         // Close the file
    return 0;
}