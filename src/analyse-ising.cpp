#include <vector>
#include <list>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

typedef std::vector<double> Vector;

double sqr(double x) {
    return x*x;
}

double sum(const Vector &values) {
    double result = 0.;
    for (auto v : values)
        result += v;
    return result;
}

double sumsqr(const Vector &values) {
    double result = 0.;
    for (auto v : values)
        result += v*v;
    return result;
}


//! This runs a jackknife error analysis on an array of magnetisations. All other
//! routines in this programm are just code to read in data, etc., so this is the
//! interesting one you might want to read:
void jackknife(
        std::ostream &out,              // output stream
        std::vector<double> &values, 	// Array and number of data points
        int blocksize,   		// Blocksize for the jackknife
        double beta, int extent)	// Beta and extent of the lattice
{
    // Number of jackknife block
    const int num_blocks = values.size() / blocksize;

    // Cut off the last datapoints so all blocks are filled
    const int count = num_blocks*blocksize;
    values.resize(count); 

    // Check we have a reasonably big data set
    if (num_blocks < 4) {
        std::cerr << " - skipped.";
        return; 
    }

    // Calculate the standard average and variance
    // These are called "T" (without index) in your script
    const double std_sum_mag = sum(values);
    const double std_sum_mag2 = sumsqr(values);

    const double std_avr_mag = std_sum_mag / count;
    const double std_var_mag = std_sum_mag2 / count - sqr(std_avr_mag);

    // Calculate the T_i and the pseudo-jackknife values J_i
    Vector Ji_avr_mag(num_blocks),
           Ji_var_mag(num_blocks);

    // Loop all jackknife blocks
    for (int j=0;j<num_blocks; j++) {

        // Calculate the average, leaving out data points from the current JK block
        // From the overall sum the jackknife block is substracted. This is significantly
        // faster than summing over all but the blocked data points each iteration
        double Ti_avr_mag = std_sum_mag;
        for (int i=j*blocksize; i<(j+1)*blocksize; i++)
            Ti_avr_mag -= values[i];
        Ti_avr_mag /= (count-blocksize);


        // Calculate the variance, leaving out data points from the current JK block
        double Ti_var_mag = std_sum_mag2;
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
    const double res_mag_avr = sum(Ji_avr_mag)/num_blocks;
    double res_sus_avr = sum(Ji_var_mag)/num_blocks;

    // Variances
    const double res_mag_var = (sumsqr(Ji_avr_mag)/num_blocks - sqr(res_mag_avr))/(num_blocks-1);
    double res_sus_var = (sumsqr(Ji_var_mag)/num_blocks - sqr(res_sus_avr))/(num_blocks-1);

    // To get the susceptibility, volume and beta have to be multiplied
    res_sus_avr *= extent*extent*beta;
    res_sus_var *= sqr(extent*extent*beta);


    // Data is outputted to file
    
    out << std::fixed<<
        ' ' << std::setw(10) << count << 
        ' ' << std::setw(10) << extent << 
        ' ' << std::setw(10) << beta << 
        ' ' << std::setw(10) << std_avr_mag << 
        ' ' << std::setw(10) << std_var_mag << 
        ' ' << std::setw(10) << res_mag_avr << 
        ' ' << std::setw(10) << sqrt(res_mag_var) << 
        ' ' << std::setw(10) << res_sus_avr << 
        ' ' << std::setw(10) << sqrt(res_sus_var) << std::endl;
}



//! Run the analysis on a beta value. This reads in the file
//! and calls the jackknife analysis on an array of magnetization
//! values. It takes the file, size of JK blocks, the lattice
//! extent and the number of thermalisation steps as parameters
void analyse_beta(std::ofstream &out, std::string filename, int blocksize, int extent, int therm) 
{ 
    // Open file
    std::ifstream f(filename);


    // magnetization values
    Vector mag;
    float beta;
    f.read((char*)&beta,sizeof(beta));

    // check for file open error
    if (f.fail()) {
        std::cerr << "# Opening file " << filename << " failed!" << std::endl;
        return; 
    }
    int v_count = 0; //number of read values
    while (true) {
        float value;
        f.read((char*)&value, sizeof(value));
        if (f.eof())
            break;
        v_count++;
        if (v_count > therm)
            mag.push_back(value);
    }


    // Run analysis
    std::cout << "# At beta " << std::fixed << beta << 
        ": read in " << std::setw(7) << v_count << 
        " values, analysing " << std::setw(7) << mag.size() ;

    jackknife( out, mag, blocksize, beta, extent );
    std::cout << std::endl;
}


//! The main routine looks for simulation results in 'data/NNN' and passes
//! all beta values to read in to "analyse_beta"
int main(int argc, char** argv)
{

    int blocksize = 20;	// The (default) size of jackknife blocks

    // Command line arguments
    if ((argc<3)||(argc>4)) {
        std::cerr << "# Usage: " << argv[0] << " [extent] [#therm] {[JK-blocksize (defaults to " << blocksize << ")]}" << std::endl;
        return -1;
    }

    // Get extent
    std::stringstream args;
    for (int i=1; i<argc; i++)
        args << argv[i] << ' ';
    int extent, therm;
    args >> extent >> therm;
    if (argc == 4)
        args >> blocksize;

    if ((extent<0)||(extent%2)||(extent>1024)) {
        std::cerr << "# Extent out of range: [0-1024], even" << std::endl;
        return -1;
    }

    if ((therm<0)||(therm>100000)) {
        std::cerr << "# Number of thermalisation steps should be >0, <100000, is " << therm << std::endl;
        return -1;
    }
    if (blocksize<1) {
        std::cerr << "# JK-blocksize shoud be > 0, is " <<  blocksize<< std::endl;
        return -1;
    }

    // open the directory
    std::stringstream dirname;
    dirname << "data/" << std::setfill('0') << std::setw(3) << extent;

    std::list<std::string> filenames;
    for (fs::directory_iterator itr(dirname.str()); itr != fs::directory_iterator(); itr++)
        if (fs::is_regular_file(itr->status()))
            filenames.push_back(itr->path().string());
    filenames.sort();
//     for (auto f: filenames)
//         std::cout << f << std::endl;

    // check if files where found
    if (filenames.empty()) {
        std::cerr << "# No files found in result/" << extent << "/" << std::endl;
        return -2;
    }

    // Start the analysis on each beta
    std::cout << "# Found " << filenames.size() << " different result files for lattice size " << extent << ", analysing..." << std::endl;

    fs::create_directory("result");

    // Open output file
    dirname.str(""); //reset
    dirname << "result/" << std::setw(3) << std::setfill('0') << extent << ".txt";
    std::ofstream f(dirname.str());

    f << "#Jackknife error analysis, blocksize " << blocksize << ", # of therm. steps " <<  therm << std::endl;
    f << "#" << 
        std::setw(10) << "NData" << 
        std::setw(11) << "Extent"<< 
        std::setw(11) << "Beta"<< 
        std::setw(11) << "Avr"<< 
        std::setw(11) << "Var"<< 
        std::setw(11) << "JK_Mag"<< 
        std::setw(11) << "JKerrM" << 
        std::setw(11) << "JK_Sus"<< 
        std::setw(11) << "JKErrS" << std::endl;
//     f << "#NData\tExtent\tBeta\t\tAvr\t\tVar\t\tJK_Mag\t\tJKerrM\t\tJK_Sus\t\tJKerrS\n" << std::endl;

    // loop available betas
    for (auto filename : filenames)
        analyse_beta(f, filename, blocksize, extent, therm);

    return 0;
}

