#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <stdexcept>

using namespace std;

// Galaxy class definition
class Galaxy{
    public:
        Galaxy(double ra, double dec, double z){
            setRA(ra);
            setDEC(dec);
            setZ(z);
        }
        ~Galaxy(){}
        void setRA(double ra){
            if (ra < 0.0 || ra > 360.0){
                throw invalid_argument("RA must be between 0.0 and 360.0 degrees.");
            }
            coordRA = ra;
        }
        void setDEC(double dec) {
            if (dec < -90.0 || dec > 90.0) {
                throw invalid_argument("Dec must be between -90.0 and 90.0 degrees.");
            }
            coordDEC = dec;
        }
        void setZ(double z) {
            if (z < 0.0) {
                throw invalid_argument("Redshift must be non-negative");
            }
            redshift = z;
        }
        double getRA() const{
            return coordRA;
        }
        double getDEC() const{
            return coordDEC;
        }
        double getZ() const{
            return redshift;
        }
        // Method to calculate the on-sky angular distance to another galaxy
        double angularDistanceTo(const Galaxy& other) const {
            // Convert degrees to radians
            double radRA1 = coordRA * M_PI / 180.0;
            double radDEC1 = coordDEC * M_PI / 180.0;
            double radRA2 = other.coordRA * M_PI / 180.0;
            double radDEC2 = other.coordDEC * M_PI / 180.0;
            // Calculate angular distance using the haversine formula
            double deltaDEC = radDEC2 - radDEC1;
            double deltaRA = radRA2 - radRA1;
            double havTheta = sin(deltaDEC/2) * sin(deltaDEC/2)
                   + cos(radDEC1) * cos(radDEC2) * sin(deltaRA/2) * sin(deltaRA/2);
            double angDist = 2 * asin(sqrt(havTheta));
            return angDist * 180.0 / M_PI;            // Convert back to degrees
        }
    private:
        double coordRA;
        double coordDEC;
        double redshift;
};

int main() {
    ifstream infile("coma_members_Jim2025.dat");           // Read the data file
    if (!infile.is_open()){
        cerr << "Error: cannot open file: coma_members_Jim2025.dat" << endl;
        return 1;
    }
    // Skip the header line starting with '#'
    string header_line;      
    getline(infile, header_line);
    double ra, dec, z;                      // Variables for the three columns
    int galaxyID = 0;                       // Initialize the galaxy id counter
    vector<Galaxy> galaxies;                // Store the galaxies in a vector
    cout.setf(ios::fixed);
    cout.precision(3);                      // Set output precision
    cout << "ID\tRA\tDEC\tz" << endl;
    // Read data and create Galaxy objects
    try {
        while (infile >> ra >> dec >> z){
            if (infile.fail()){
                throw runtime_error("Error reading galaxy data from file");
            }
            galaxyID++;                         // Use a line counter as the ID
            galaxies.emplace_back(ra, dec, z);  // Store galaxy object in vector
            cout << "Galaxy object created and stored: ID = " << galaxyID
                 << ", RA = " << galaxies.back().getRA() << ", DEC = "
                 << galaxies.back().getDEC() << ", z = "<< galaxies.back().getZ()
                 << endl;
        }
        // Variables to track the closest galaxy pair
        double minDistance = 180.0; // Start with the maximum angular distance possible
        size_t minI = 0, minJ = 0;
        // Calculate angular distances between all pairs
        size_t galaxyPairs = (galaxies.size() * (galaxies.size() - 1)) / 2;
        cout << "\nCalculating angular distances between " << galaxies.size() 
             << " galaxies (" << galaxyPairs << " pairs)...\n" << endl;
        cout.precision(10);
        for (size_t i = 0; i < galaxies.size(); i++) {
            for (size_t j = i + 1; j < galaxies.size(); j++) {
                double currentDist = galaxies[i].angularDistanceTo(galaxies[j]);
                if (currentDist < minDistance) {
                    minDistance = currentDist;        // Update minimum distance
                    minI = i;                         // Update indices
                    minJ = j;
                    // Print alert when new minimum is found
                    cout << "New minimum found: " << minDistance
                         << " degrees (between galaxies #" << i+1 << " and #" << j+1
                         << ")" << endl;
                }
            }
        }
        // Display results
        cout << "\n=== FINAL RESULTS ===\n" << endl;
        cout << "Total galaxies: " << galaxies.size() << endl;
        cout << "Closest galaxies: " << endl;
        cout << "- Galaxy #" << minI+1 << ": RA = " << galaxies[minI].getRA() 
             << ", Dec = " << galaxies[minI].getDEC() << ", z = "
             << galaxies[minI].getZ() << endl;
        cout << "- Galaxy #" << minJ+1 << ": RA = " << galaxies[minJ].getRA() 
             << ", Dec = " << galaxies[minJ].getDEC() << ", z = "
             << galaxies[minJ].getZ() << endl;
        cout << "Minimum angular distance: " << minDistance << " degrees" << endl;
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    infile.close();
    return 0;
}