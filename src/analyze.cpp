#include <vector>
#include <list>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <math.h>

namespace fs = boost::filesystem;

typedef double floatT;
typedef float floatTdisk;
typedef std::vector<floatT> Vector;

template<class T>
T sqr(T x) {
    return x*x;
}

floatT sum(const Vector &values) {
    floatT result = 0.;
    for (int i=0; i<values.size(); i++)
        result += values[i];
    return result;
}

floatT sumsqr(const Vector &values) {
    floatT result = 0.;
    for (int i=0; i<values.size(); i++)
        result += values[i]*values[i];
    return result;
}


//! This runs a jackknife error analysis on an array of magnetisations. All other
//! routines in this programm are just code to read in data, etc., so this is the
//! interesting one you might want to read:
void jackknife(
        Vector &values, 	// Array and number of data points
        int num_blocks,   		// Number of jackknife blocks
        floatT p,                       // this is 1-exp(-4beta)
        floatT beta, int extent)	// Beta and extent of the lattice
{
    // Size of jackknife blocks
    const int blocksize = values.size() / num_blocks;

    // Cut off the last datapoints so all blocks are filled
    const int count = num_blocks*blocksize;
    values.resize(count); 

    // Check we have a reasonably big data set
    if (blocksize < 4) {
        std::cerr << " - skipped.";
        return; 
    }

    // Calculate the standard average and variance
    // These are called "T" (without index) in your script
    const floatT std_sum_mag = sum(values);
    const floatT std_sum_mag2 = sumsqr(values);

    const floatT std_avr_mag = std_sum_mag / count;
    const floatT std_var_mag = std_sum_mag2 / count - sqr(std_avr_mag);

    // Calculate the T_i and the pseudo-jackknife values J_i
    Vector Ji_avr_mag(num_blocks),
           Ji_var_mag(num_blocks);

    // Loop all jackknife blocks
    for (int j=0;j<num_blocks; j++) {

        // Calculate the average, leaving out data points from the current JK block
        // From the overall sum the jackknife block is substracted. This is significantly
        // faster than summing over all but the blocked data points each iteration
        floatT Ti_avr_mag = std_sum_mag;
        for (int i=j*blocksize; i<(j+1)*blocksize; i++)
            Ti_avr_mag -= values[i];
        Ti_avr_mag /= (count-blocksize);


        // Calculate the variance, leaving out data points from the current JK block
        floatT Ti_var_mag = std_sum_mag2;
        for (int i = j*blocksize; i<(j+1)*blocksize; i++)
            Ti_var_mag -= sqr(values[i]);
        Ti_var_mag /= count-blocksize;
        Ti_var_mag -= sqr(Ti_avr_mag);

        // Calculate the pseudo-jackknife values for average (mag) and variance (sus):
        Ji_avr_mag[j] = num_blocks * std_avr_mag - (num_blocks-1) * Ti_avr_mag;
        Ji_var_mag[j] = num_blocks * std_var_mag - (num_blocks-1) * Ti_var_mag;
    }

    // Do a standard average and variance on the pseudo jackknife values for both
    // the average/magnetisation and the variance/susceptibility:

//     double res_mag_avr=0, res_mag_var=0, res_sus_avr=0, res_sus_var=0;

    // Averages
    const floatT res_mag_avr = sum(Ji_avr_mag)/num_blocks;
    floatT res_sus_avr = sum(Ji_var_mag)/num_blocks;

    // Variances
    const floatT res_mag_var = (sumsqr(Ji_avr_mag)/num_blocks - sqr(res_mag_avr))/(num_blocks-1);
    floatT res_sus_var = (sumsqr(Ji_var_mag)/num_blocks - sqr(res_sus_avr))/(num_blocks-1);

    // To get the susceptibility, volume and beta have to be multiplied
    res_sus_avr *= extent*extent*beta;
    res_sus_var *= sqr(extent*extent*beta);


    // Data is outputted to file
    
    std::cout << std::fixed<<
        ' ' << std::setw(10) << count << 
        ' ' << std::setw(10) << extent << 
        ' ' << std::setw(10) << p << 
        ' ' << std::setw(10) << beta << 
        ' ' << std::setw(10) << std_avr_mag << 
        ' ' << std::setw(10) << std_var_mag << 
        ' ' << std::setw(10) << res_mag_avr << 
        ' ' << std::setw(10) << sqrt(res_mag_var) << 
        ' ' << std::setw(10) << res_sus_avr << 
        ' ' << std::setw(10) << sqrt(res_sus_var) << std::endl;
}


floatT get_p(std::string sig) {
    floatT result = 0;
    floatT val = 1;
    for (int i=0; i<sig.size(); i++) {
        val /= 2;
        if (sig[i] == '1')
            result += val;
    }
    return result;
}

//! Run the analysis on a beta value. This reads in the file
//! and calls the jackknife analysis on an array of magnetization
//! values. It takes the file, size of JK blocks, the lattice
//! extent and the number of thermalisation steps as parameters
// void analyse_beta(std::ofstream &out, std::string filename, int blocksize, int extent, int therm) 
int main(int argc, char** argv)
{ 
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " blockcount therm filename" << std::endl;
        return 1;
    }
    int blockcount = atoi(argv[1]);
    float therm = atof(argv[2]);
    std::string filename = argv[3];

    fs::path path(filename);
    const floatT p = get_p(path.filename().c_str());
    const floatT beta = -0.25*log(1-p);
    const int extent = atoi(path.parent_path().filename().c_str());


    //some checks
    if ((extent<=0)||(extent%2)) {
        std::cerr << "# Extent is not positive and even" << std::endl;
        return -1;
    }

    if ((therm<0)||(therm>1)) {
        std::cerr << "# percentage of thermalisation should be >0, <1, is " << therm << std::endl;
        return -1;
    }
    if (blockcount<1) {
        std::cerr << "# JK-blockcount shoud be > 0, is " <<  blockcount << std::endl;
        return -1;
    }
    // Open file
    std::ifstream f(filename.c_str());


    // magnetization values
    Vector mag;

    // check for file open error
    if (f.fail()) {
        std::cerr << "# Opening file " << filename << " failed!" << std::endl;
        return 1; 
    }
    while (true) {
        floatTdisk value;
        f.read((char*)&value, sizeof(value));
        if (f.eof())
            break;
        mag.push_back(value);
    }
    const int skip = mag.size()*therm;
    const int count = ((mag.size()-skip)/blockcount)*blockcount;
    for (int i=0; i<count; i++)
        mag[i] = mag[i+skip];

    // Run analysis
    std::cerr << "# At beta " << std::fixed << beta << 
        ": read in " << std::setw(7) << mag.size() << 
        " values, analysing " << std::setw(7) << count;
    mag.resize(count);

    jackknife( mag, blockcount, p, beta, extent );
    std::cerr << std::endl;
}

