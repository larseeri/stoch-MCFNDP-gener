// Copyright (c) 2025 Eric Larsen
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the “Software”), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
// ****************************************************************************
// ****************************************************************************
//
//
#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include <cstring>
#include <filesystem>
// #include <boost/filesystem.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "../HOYLAND_KAUT_WALLACE/matrix.h"

using namespace std;

// prefix of all parameters, except configuration file related
#define PARAM_PREFIX '-'

// prefix of parameters related to configuration file and to reading
// of configuration file
#define CONFIG_PREFIX '+'

// prefix of non-default configuration file name
#define CONFIG_FILE_KEY "+F"

// requesting verbose mode when looking for configuration file
#define CONFIG_VERBOSE_MODE "+v"

// maximum length of lines in configuration file
const int lineLength = 256;

// input and output formats available
enum IO_FORMAT
{
    DOW = 0,
    GSTD, // generic STD
    RSTD, // restricted STD
};

// distributional characteristics of moments available
enum MOMENTS_DIST_CHARACT
{
    UNIFORM = 0,
    TRIANGULAR,
    UNKNOWN
};

// Types of stochastic elements (i.e. commodity demands, total arc
// capacities, capacities per arc and commodity, fixed costs per arc,
// variable costs per arc and commodity) that can be randomized.
// MAX_TYPE will be used as an index upper bound when looping over
// STOCH_ELEMS_TYPE. ALL_TYPES is used to identify arrays containing
// stochastic elements from all types.
enum STOCH_ELEMS_TYPE
{
    DEMAND = 0,
    CAP_ARC,
    CAP_COMMOD,
    FIXED_COST,
    VAR_COST,
    MAX_TYPE,
    ALL_TYPES = MAX_TYPE
};

const string formats[] = {"dow", "std"};

const string elements[] = {"demand", "capacity of arc", "capacity of comm. on arc", "fixed cost of arc",
                           "variable cost of comm. on arc"};

// Binary flags signaling which elements are to be randomized (i.e. are to
// vary across generated scenarios). Their sum is assigned to configuration
// parameter X.
const int flagStochDemand = 1;
const int flagStochArcCapac = 2;
const int flagStochCommodOnArcCapac = 4;
const int flagStochFixedCost = 8;
const int flagStochVarCost = 16;

// Holding the elements belonging to a particular type in STOCH_ELEMS_TYPE.
struct StochElem
{
    // pointer of array containing the elements
    double *ptr;

    // type of the stochastic element
    STOCH_ELEMS_TYPE type;

    // size of array
    uint size;

    // operators retrieving individual elements stored in the array
    double &operator[](int index)
    {
        return ptr[index];
    }
    const double &operator[](int index) const
    {
        return ptr[index];
    }
};

typedef vector<StochElem> VectElem;
typedef vector<STOCH_ELEMS_TYPE> VectStoch;
typedef vector<double> VectDbl;
typedef vector<VectDbl> MatDbl;

// Options for running the program can be specified directly in the command
// line or in a configuration file. If the name of a configuration file is
// specified in the command line with
//
// ./exe +F nameOfConfigFile
//
// and if this file exists, then an attempt to process its content is made.
// Otherwise, the file whose name is specified in the constructor in line
// 367 of main.cpp is processed, provided it exists. Verbosity option +v can
// be used in the command line to obtain information about all options that
// have been set in the parameter file or in the command line.
//
// The configuration file should have this format:
//
// optionKey1   value1
// optionKey2   value2
// .
// .
// optionKeyN   valueN
//
// It may specify at most one option per line. Option keys may be prefixed or
// not with '-'. Comments may appear either at the beginning of a line, in
// which case they must start with a non-alphabetic character (but not with
// '-'), or at the end of a line. Once the configuration file has been
// processed, options appearing in the command line (other than +F and +v)
// are processed. Option values specified in the command line supersede those
// specified in the configuration file. That is, if present, they replace,
// modify and supplement the latter. Option keys appearing in the command line
// must be prefixed with '-' as follow:
//
// -optionkeyA valueA -optionKeyB valueB...
//
// By default, boolean values are set to false and presence of option key is
// sufficient to set the corresponding value to true. For other options
// expecting string, integer or double values, values should be specified.
//
// (See main.cpp for details on configuration parameters and construction
// of command lines and configuration files.)
struct Configuration
{

    // number of arguments appearing in command line
    int nbArg;

    // char array containing the command line
    char **vectArg;

    // flags displaying parameter values when parsing
    bool verbose;

    // name of configuration file
    string configFileName;

    // format of input and output instances: DOW, GSTD or RSTD
    IO_FORMAT instanceFormat;

    // flag indicating that correlations are either directly specified in
    // blocks or read from file
    bool specCorrBlocks;

    // distributional characteristics of moments:
    // uniform or triangular
    MOMENTS_DIST_CHARACT momentsCharacts;

    // flag indicating to either generate or read moments
    bool generMoments;

    // alpha parameter used in generating moments
    double alpha;

    // beta parameter used in generating moments
    double beta;
    //
    // When MOMENTS_DIST_CHARACT == UNIFORM, then the stochastic elements
    // of the network that are changed according to scenarios are each assumed
    // to follow a uniform(a, b) distribution where,
    // a = D - (alpha * D),
    // b = D + (beta * D),
    // D is the base value of the stochastic element taken from the base
    // deterministic network,
    // alpha in [0, 1),
    // beta in [0, \infty).
    // The distribution is symmetric around D when alpha = beta.
    //
    // When MOMENTS_DIST_CHARACT == TRIANGULAR, then the stochastic elements
    // of the network that are changed according to scenarios are each assumed
    // to follow a triangular(a, b, c) distribution where,
    // where,
    // a = D - (alpha * D),
    // b = D + (beta * D),
    // c = D,
    // D is the base value of the stochastic element taken from the base
    // deterministic network,
    // alpha in [0, 1),
    // beta in [0, \infty).
    // The distribution is symmetric around D when alpha = beta.
    //
    // name of file containing base determistic MCFND network (default or
    // modified with option key I)
    string inputInstFileName;

    // name of file containing generated MCFND networks (default or modified
    // with option key O)
    string outputInstFileName;

    // name of file containing detailed moments (default or modified with
    // option M); read from when option key G is absent; written to when
    // option key G is present
    string momentsFileName;

    // name of file containing detailed correlations (default or modified
    // with option key C); read from when option key K is absent; written
    // to when option key K is present
    string correlsFileName;

    // vector listing the types specified as stochastic (i.e. elements
    // belonging to these types will vary across generated scenarios)
    VectStoch stochTypes;

    // vector listing the types specified as non-stochastic (i.e. elements
    // belonging to these types will remain equal to their value in the base
    // instance across generated scenarios)
    VectStoch nonStochTypes;

    // matrix of size (MAX_TYPE, MAX_TYPE) specifying uniform correlation
    // coefficient between distinct elements within a type and between elements
    // across two types; correlation between elements of two types is zero when
    // at least one of the types is non-stochastic
    MatDbl correlCoeffsMatrix;
    //
    // The following parameters are supplied to the Hoyland-Kaut-Wallace
    // algorithm:
    //
    // number of scenarios
    int numbScenarios;

    // characteristics of moments (sum of bits)
    // To ensure expected working of program, distCharactMoments(7), 7=1+2+4,
    // SHOULD be supplied to Hoyland-Kaut-Wallace algorithm, where:
    // 1 -> use population estimators (as in spreadsheets),
    // 2 -> 2nd moment is Var instead of StDev,
    // 4 -> 4th moment is Kurtosis - 3,
    // Note: by default, 3rd and 4th moments are scaled by StDev.
    int distCharactMoments;

    // max error allowed when matching moments (scaled to var=1)
    double maxErrorMoment;

    // max error allowed when matching correlations (scaled to var=1)
    double maxErrorCorrel;

    // amount of information displayed by the program
    int displayLevel;

    // max number of attempts to generate scenarios
    int maxTrial;

    // max number of iterations in HKW algorithm
    int maxIter;

    // random seed (pcg_random)
    int randomSeed;

    // random stream (pcg_random)
    int randomStream;

    // name of file where resulting matrix of scenarios are written in HKW
    // format (default or modified with option key o)
    string outScenariosFileName;

    // flag indicating to read probabilites from file
    bool readProbs;

    // name of input file containing probabilities (set with option key p,
    // otherwise equiprobable scenarios)
    string probsFileName;

    // flag indicating to read starting scenarios from file
    bool readStartScenarios;

    // name of file where a matrix of scenarios saved in HKW format in a
    // previous run (default or modified with option key o) are read and used
    // as starting values (default or modified with option key s)
    string startScenariosFileName;

    // for validation purposes
    bool paramX;
    bool paramT;

    Configuration(int argc, char **argv, string configFileName);

    void printHelpAndExit();

    void printUsageAndExit(char ExecName[], bool anonymous);

    // While traversing command line, attempts to find (i) name of
    // of a configuation file and (ii) request of verbosity.
    void parseConfigFileNameFromCmdLine();

    // Processes keys or (key, value) pairs specified on non-comment line of
    // configuration file.
    bool processConfigFile();

    // Parses command line.
    void parseCmdLine();

    // Parses key or (key, value) pair found on a line of the configuration
    // file.
    void parseKeyValue(string key, string value);

    // Removes quotes from value in (key, value) pair on a non-comment line
    // of the configuration file.
    string removeQuotes(istringstream &configLine);

    // Verifies configuration values.
    void verify();
	
	// Verifies that output path is valid.
	void verifyOutputPath(string);
	
};

#endif
