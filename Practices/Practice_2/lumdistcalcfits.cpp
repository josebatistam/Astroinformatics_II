#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <CCfits/CCfits>

using namespace std;
using namespace CCfits;

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
    try {
        // Read the FITS file
        string inputFile = "Coma_members_Jim2025.fits";
        FITS* pInFits = new FITS(inputFile, Read, true);
        
        // Get the table
        ExtHDU& table = pInFits->extension(1);
        
        // Read data columns
        vector<double> ra, dec, z;
        int nrows = table.rows();
        
        table.column("RAJ2000").read(ra, 1, nrows);
        table.column("DEJ2000").read(dec, 1, nrows);
        table.column("z").read(z, 1, nrows);
        
        // Calculate distances
        vector<double> distances;
        for (int i = 0; i < nrows; ++i) {
            distances.push_back(calculateDL(z[i]));
        }
        
        // Create output FITS file
        string outputFile = "Coma_members_Jim2025_v2.fits";
        
        // Remove if exists
        remove(outputFile.c_str());
        
        FITS* pOutFits = new FITS(outputFile, Write);
        
        // Define columns
        vector<string> colNames = {"RAJ2000", "DEJ2000", "z", "Distance"};
        vector<string> colFormats = {"D", "D", "D", "D"};
        vector<string> colUnits = {"deg", "deg", "", "Mpc"};
        
        // Create table
        Table* newTable = pOutFits->addTable("GALAXIES", nrows, colNames, colFormats, colUnits);
        
        // Write data
        newTable->column("RAJ2000").write(ra, 1);
        newTable->column("DEJ2000").write(dec, 1);
        newTable->column("z").write(z, 1);
        newTable->column("Distance").write(distances, 1);
        
        // Add some header info
        newTable->addKey("H0", 70.0, "Hubble constant");
        newTable->addKey("OMEGA_M", 0.3, "Matter density");
        newTable->addKey("OMEGA_L", 0.7, "Dark energy density");
        
        cout << "Created " << outputFile << " with " << nrows << " galaxies" << endl;
        cout << "Distance column added in Mpc" << endl;
        
        // Clean up
        delete pInFits;
        delete pOutFits;
        
    } catch (FITS::CantOpen&) {
        cerr << "Error: cannot open FITS file" << endl;
        return 1;
    } catch (CCfits::FitsException& e) {
        cerr << "FITS Error: " << e.message() << endl;
        return 1;
    }
    
    return 0;
}