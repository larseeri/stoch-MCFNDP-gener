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
#include "configuration.h"

// Constructor
Configuration::Configuration(int argc, char **argv, string _configFileName)
    : nbArg(argc), vectArg(argv), verbose(false), configFileName(_configFileName), instanceFormat(DOW),
      specCorrBlocks(false), momentsCharacts(UNIFORM), generMoments(false), alpha(0.25), beta(0.25),
      inputInstFileName("./base_network.dow"), outputInstFileName("./generated_networks.txt"),
      momentsFileName("./target_moms.dat"), correlsFileName("./target_corrs.dat"), numbScenarios(1000),
      //
      // distCharactMoments = 7 (7=1+2+4), ensures that the target moments supplied
      // to HKW_ScenGen() in HKW_sg.c are properly interpreted.
      // 1 -> use population estimators (as in spreadsheets),
      // 2 -> 2nd moment is Var instead of StDev,
      // 4 -> 4th moment is Kurtosis - 3,
      // Note: by default, 3rd and 4th moments are scaled by StDev.
      // This MUST be set to 7.
      distCharactMoments(7),
      //
      maxErrorMoment(1e-3), maxErrorCorrel(1e-3),
      //
      // displayLevel = 3 ensures that sufficient information is reported about
      // the generation process by HKW_ScenGen() in HKW_sg.c.
      // This should be set to 3.
      displayLevel(3),
      //
      maxTrial(20), maxIter(50), randomSeed(7654321), randomStream(1000), outScenariosFileName("raw_scens.dat"),
      readProbs(false), probsFileName(""), readStartScenarios(false), startScenariosFileName(""), paramX(false),
      paramT(false)
{
    correlCoeffsMatrix.resize(MAX_TYPE, VectDbl(MAX_TYPE, 0.0));

    if (nbArg == 2 && (strcmp(vectArg[1], "-h") == 0 || strcmp(vectArg[1], "-H") == 0))
        printHelpAndExit();

    parseConfigFileNameFromCmdLine();
    processConfigFile();
    parseCmdLine();
    verify();
}

// Parses command line.
void Configuration::parseCmdLine()
{
    int elemStoch = DEMAND;

    bool reinitCorrMatrix = true;

    int a = 1;
    const char *key;
    const char *value;

    cout << "****** Parsing configuration parameters from the command line.\n\n";

    while (a < nbArg)
    {
        if ((vectArg[a][0] == PARAM_PREFIX))
        {

            key = &vectArg[a][1];

            if (a < nbArg - 1)
            {
                value = vectArg[a + 1];

                if (((vectArg[a + 1][0] == PARAM_PREFIX) || (vectArg[a + 1][0] == CONFIG_PREFIX)) &&
                    isalpha(vectArg[a + 1][1]))
                    value = "";
                else
                    a++;
            }
            else
                value = "";

            if (verbose)
                cout << "****** command line -- parameter [" << key << "] = [" << value << "]\n" << endl;

            switch (*key)
            {
            case 'F':
                if (*value == 'D')
                {
                    instanceFormat = DOW;
                    break;
                }
                else if (*value == 'G')
                {
                    instanceFormat = GSTD;
                    break;
                }
                else if (*value == 'R')
                {
                    instanceFormat = RSTD;
                    break;
                }
                printUsageAndExit(vectArg[0], false);
                break;
            case 'T':
                if (*value == 'U')
                {
                    momentsCharacts = UNIFORM;
                    paramT = true;
                    break;
                }
                else if (*value == 'T')
                {
                    momentsCharacts = TRIANGULAR;
                    paramT = true;
                    break;
                }
                printUsageAndExit(vectArg[0], false);
                break;
            case 'S':
                elemStoch = atoi(value);

                if (elemStoch & flagStochDemand)
                    stochTypes.push_back(DEMAND);
                else
                    nonStochTypes.push_back(DEMAND);
                if (elemStoch & flagStochArcCapac)
                    stochTypes.push_back(CAP_ARC);
                else
                    nonStochTypes.push_back(CAP_ARC);
                if (elemStoch & flagStochCommodOnArcCapac)
                    stochTypes.push_back(CAP_COMMOD);
                else
                    nonStochTypes.push_back(CAP_COMMOD);
                if (elemStoch & flagStochFixedCost)
                    stochTypes.push_back(FIXED_COST);
                else
                    nonStochTypes.push_back(FIXED_COST);
                if (elemStoch & flagStochVarCost)
                    stochTypes.push_back(VAR_COST);
                else
                    nonStochTypes.push_back(VAR_COST);
                break;
            case 'X':
                // whenever any information about correlations is supplied in
                // command line, any other information previously supplied
                // in configuration file is disregarded
                if (reinitCorrMatrix)
                {
                    correlCoeffsMatrix.resize(MAX_TYPE, VectDbl(MAX_TYPE, 0.0));
                    reinitCorrMatrix = false;
                }
                STOCH_ELEMS_TYPE elem[2];
                if (strlen(vectArg[a]) == 4)
                {
                    for (int j = 0; j < 2; j++)
                        switch (vectArg[a][j + 2])
                        {
                        case 'D':
                            elem[j] = DEMAND;
                            break;
                        case 'A':
                            elem[j] = CAP_ARC;
                            break;
                        case 'C':
                            elem[j] = CAP_COMMOD;
                            break;
                        case 'F':
                            elem[j] = FIXED_COST;
                            break;
                        case 'V':
                            elem[j] = VAR_COST;
                            break;
                        default:
                            printUsageAndExit(vectArg[0], false);
                            break;
                        }
                    paramX = true;
                    correlCoeffsMatrix[elem[0]][elem[1]] = atof(value);
                    correlCoeffsMatrix[elem[1]][elem[0]] = atof(value);
                    break;
                }
                else
                {
                    printUsageAndExit(vectArg[0], false);
                    break;
                }
            case 'G':
                generMoments = true;
                //--i;
                break;
            case 'M':
                momentsFileName = value;
                break;
            case 'K':
                specCorrBlocks = true;
                //--i;
                break;
            case 'C':
                correlsFileName = value;
                break;
            case 'A':
                alpha = atof(value);
                break;
            case 'B':
                beta = atof(value);
                break;
            case 'I':
                inputInstFileName = value;
                break;
            case 'O':
                outputInstFileName = value;
                break;

            // the following arguments are supplied to the
            // Hoyland-Kaut-Wallace algorithm
            case 'n':
                numbScenarios = atoi(value);
                break;
            case 't':
                maxTrial = atoi(value);
                break;
            case 'i':
                maxIter = atoi(value);
                break;
            case 'm':
                maxErrorMoment = atof(value);
                break;
            case 'c':
                maxErrorCorrel = atof(value);
                break;
            case 'l':
                displayLevel = atoi(value);
                break;
            case 'r':
                randomSeed = atoi(value);
                break;
            case 'R':
                randomStream = atoi(value);
                break;
            case 'o':
                outScenariosFileName = value;
                break;
            case 'p':
                readProbs = true;
                probsFileName = value;
                break;
            case 's':
                readStartScenarios = true;
                startScenariosFileName = value;
                break;

            // help
            case 'H':
                printHelpAndExit();
                break;
            case 'h':
                printHelpAndExit();
                break;
            default:
                ////printUsageAndExit(vectArg[0], false);
                break;
            }
        }

        a++;
    }

    if (verbose)

        cout << endl;
}

// Parses key or (key, value) pair found on a line of the configuration file.
inline void Configuration::parseKeyValue(string key, string value)
{
    int elemStoch = DEMAND;

    switch (key[0])
    {
    case 'F':
        if (value == "D")
        {
            instanceFormat = DOW;
            break;
        }
        else if (value == "G")
        {
            instanceFormat = GSTD;
            break;
        }
        else if (value == "R")
        {
            instanceFormat = RSTD;
            break;
        }
        printUsageAndExit(0, true);
        break;
    case 'T':
        if (value == "U")
        {
            momentsCharacts = UNIFORM;
            paramT = true;
            break;
        }
        else if (value == "T")
        {
            momentsCharacts = TRIANGULAR;
            paramT = true;
            break;
        }
        printUsageAndExit(0, true);
        break;
    case 'S':
        elemStoch = stoi(value);

        if (elemStoch & flagStochDemand)
            stochTypes.push_back(DEMAND);
        else
            nonStochTypes.push_back(DEMAND);
        if (elemStoch & flagStochArcCapac)
            stochTypes.push_back(CAP_ARC);
        else
            nonStochTypes.push_back(CAP_ARC);
        if (elemStoch & flagStochCommodOnArcCapac)
            stochTypes.push_back(CAP_COMMOD);
        else
            nonStochTypes.push_back(CAP_COMMOD);
        if (elemStoch & flagStochFixedCost)
            stochTypes.push_back(FIXED_COST);
        else
            nonStochTypes.push_back(FIXED_COST);
        if (elemStoch & flagStochVarCost)
            stochTypes.push_back(VAR_COST);
        else
            nonStochTypes.push_back(VAR_COST);
        break;
    case 'X':
        STOCH_ELEMS_TYPE elem[2];
        if (key.length() == 3)
        {
            for (int j = 0; j < 2; j++)
                switch (key[j + 1])
                {
                case 'D':
                    elem[j] = DEMAND;
                    break;
                case 'A':
                    elem[j] = CAP_ARC;
                    break;
                case 'C':
                    elem[j] = CAP_COMMOD;
                    break;
                case 'F':
                    elem[j] = FIXED_COST;
                    break;
                case 'V':
                    elem[j] = VAR_COST;
                    break;
                default:
                    printUsageAndExit(0, true);
                    break;
                }
            paramX = true;
            correlCoeffsMatrix[elem[0]][elem[1]] = stof(value);
            correlCoeffsMatrix[elem[1]][elem[0]] = stof(value);
            break;
        }
        else
        {
            printUsageAndExit(0, true);
            break;
        }
    case 'G':
        generMoments = true;
        break;
    case 'M':
        momentsFileName = value;
        break;
    case 'K':
        specCorrBlocks = true;
        break;
    case 'C':
        correlsFileName = value;
        break;
    case 'A':
        alpha = stof(value);
        break;
    case 'B':
        beta = stof(value);
        break;
    case 'I':
        inputInstFileName = value;
        break;
    case 'O':
        outputInstFileName = value;
        break;

    // the following arguments are supplied to the
    // Hoyland-Kaut-Wallace algorithm
    case 'n':
        numbScenarios = stoi(value);
        break;
    case 't':
        maxTrial = stoi(value);
        break;
    case 'i':
        maxIter = stoi(value);
        break;
    case 'm':
        maxErrorMoment = stof(value);
        break;
    case 'c':
        maxErrorCorrel = stof(value);
        break;
    case 'l':
        displayLevel = stoi(value);
        break;
    case 'r':
        randomSeed = stoi(value);
        break;
    case 'R':
        randomStream = stoi(value);
        break;
    case 'o':
        outScenariosFileName = value;
        break;
    case 'p':
        readProbs = true;
        probsFileName = value;
        break;
    case 's':
        readStartScenarios = true;
        startScenariosFileName = value;
        break;

    // help
    case 'H':
        printHelpAndExit();
        break;
    case 'h':
        printHelpAndExit();
        break;
    default:
        ////printUsageAndExit(0, true);
        break;
    }
}

// Verifies configuration values.
void Configuration::verify()
{
	verifyOutputPath(outputInstFileName);
	
	verifyOutputPath(outScenariosFileName);
	
    if (inputInstFileName == "")
    {
        cerr << "\nERROR:\n";
        cerr << "No deterministic base instance file was specified.\n";
        cerr << "Use option -I" << endl;
        exit(1);
    }
    else if (!filesystem::exists(inputInstFileName))
    {
        cerr << "\nERROR:\n";
        cerr << "The specified deterministic base instance file could not be found.\n";
        cerr << "Use option -I" << endl;
        exit(1);
    }
    if (stochTypes.size() == 0)
    {
        cerr << "\nERROR:\n";
        cerr << "No stochastic type was specified." << endl;
        cerr << "Use option -S" << endl;
        exit(1);
    }

    if (alpha < 0.0 || alpha > 0.99)
    {
        cerr << "\nERROR:\n";
        cerr << "alpha must fall in [0, 0.99]." << endl;
        exit(1);
    }

    if (beta < 0.0)
    {
        cerr << "\nERROR:\n";
        cerr << "beta must be non-negative." << endl;
        exit(1);
    }

    if (!specCorrBlocks)
    {		
		if (paramX)
		{		
        cerr << "\nERROR:\n";
        cerr << "Option(s) X?? cannot be used without option K." << endl;
        exit(1);
		}
	}	
	else
		verifyOutputPath(correlsFileName);

    if (!generMoments)
    {		
		if (paramT)
		{		
			cerr << "\nERROR:\n";
			cerr << "Option T cannot be used without option G." << endl;
			exit(1);
		}
	}
	else
		verifyOutputPath(momentsFileName);

    for (int row = 0; row < MAX_TYPE; row++)
    {
        for (int col = 0; col < MAX_TYPE; col++)
        {
            if (correlCoeffsMatrix[row][col] < -1.0 || correlCoeffsMatrix[row][col] > 1)
            {
                cerr << "\nERROR:\n";
                cerr << "Correlations must fall in [-1, 1]." << endl;
                exit(1);
            }
        }
    }
	
}

void Configuration::verifyOutputPath(string fullPath)
{
	size_t last_sep = fullPath.find_last_of("/");

	if(last_sep != string::npos)
	{
		string folderPath = fullPath.substr(0, last_sep + 1);	
	
		if (!filesystem::exists(folderPath) || !filesystem::is_directory(folderPath))			
		{			
			cerr << "\nERROR: " << endl;    
			cerr << "Output file path " << fullPath << " was specified.\n" 
				<< "However, folder " << folderPath << " is unreachable.\n"
				<< "Make sure to create individual folders on output file path and to grant sufficient access." << endl;
			cerr << "EXECUTION ABORTED\n" << endl;
			exit(1);
		}
	}	
}

// While traversing command line, attempts to find (i) name of
// of a configuation file and (ii) request of verbosity.
void Configuration::parseConfigFileNameFromCmdLine()
{
    bool foundFile = false;

    for (int a = 1; a < nbArg; a++)
    {
        if (vectArg[a][0] == CONFIG_PREFIX)
        {
            if (strcmp(vectArg[a], CONFIG_FILE_KEY) == 0)
            {
                // checks if configuration file name is specified on cmd line
                if (a < nbArg - 1)
                {
                    configFileName = vectArg[a + 1];
                    foundFile = true;
                }
            }
            else if (strcmp(vectArg[a], CONFIG_VERBOSE_MODE) == 0)
            {
                verbose = true;
            }
            //
            else if (isalpha(vectArg[a][1]))
                cerr << "****** configuration prefix [" << vectArg[a] << "] unknown.\n";
        }
    }
    if (verbose)
    {
        cout << "\n********* VERBOSE MODE *********\n\n";

        if (foundFile)
            cout << "****** configuration file name [" << CONFIG_FILE_KEY << "] = [" << configFileName << "]\n";

        cout << endl;
    }
}

// Processes keys or (key, value) pairs specified on a non-comment lines of
// configuration file.
bool Configuration::processConfigFile()
{
    char line[lineLength];
    ifstream inStr;

    // checking if configuration file exists
    if (!filesystem::is_regular_file(configFileName))
    {
        cerr << "\nERROR:" << endl;
        cerr << "The specified configuration file " << configFileName << " could not be found." << endl;
        cerr << "EXECUTION ABORTED\n" << endl;
        exit(1);
    }

    inStr.open(configFileName);

    if (!inStr.is_open())
    {
        cerr << "\nERROR:" << endl;
        cerr << "The specified configuration file " << configFileName << " could not be opened." << endl;
        cerr << "EXECUTION ABORTED\n" << endl;
        exit(1);
    }

    bool found = 0;

    if (verbose)
        cout << "\n****** Parsing this configuration file: " << configFileName << "]\n\n";

    inStr.getline(line, lineLength);
    while (inStr)
    {
        istringstream configLine(line);

        string key;
        configLine >> key;
        string value = removeQuotes(configLine);

        int pos = (key[0] == PARAM_PREFIX) ? 1 : 0;

        // ensures that it is not a comment
        if (isalpha(key[pos]))
        {
            found = 1;
            if (verbose)
                cout << "****** configuration file -- parameter [" << key.substr(pos) << "] = [" << value << "]\n";
            parseKeyValue(key.substr(pos), value);
        }
        inStr.getline(line, lineLength);
    }
    if (verbose)
        cout << endl;
    return found;
}

// Removes quotes from value in (key, value) pair on a non-comment line
// of the configuration file.
inline string Configuration::removeQuotes(istringstream &configLine)
{
    char caract = '\0';
    string outString("");

    configLine >> caract;

    if (configLine)
    {
        if (caract == '"')
        {
            char temp0[lineLength];
            temp0[0] = '\0';
            configLine.getline(temp0, lineLength, '"');
            outString = temp0;
        }
        else
        {
            string temp1;
            configLine >> temp1;
            outString = caract + temp1;
        }
    }

    return outString;
}

void Configuration::printHelpAndExit()
{

    printf("\nHELP REQUESTED: DISPLAYING DOCUMENTATION\n");

    printf("\nSCENARIOS GENERATOR FOR THE LINEAR, TWO-STAGE, MULTI-COMMODITY, FIXED\n");
    printf("CHARGE, CAPACITATED NETWORK DESIGN PROBLEM WITH STOCHASTIC RECOURSE\n");
    //
    printf("\nThis application randomly generates finite sets of random realizations,\n");
    printf("i.e., scenarios, while matching target moments and correlations, for\n");
    printf("specified subsets of parameters of a given instance of the deterministic\n");
    printf("multi-commodity, fixed charge, capacitated network design (MCFND) problem.\n");
    printf("Details of this process are governed by user-specified configuration\n");
    printf("settings. The code implementing the algorithm presented in Hoyland, Kaut\n");
    printf("and Wallace (2003) (HKW) and built by Michal Kaut and Diego Mathieu is used\n");
    printf("to compute the scenarios.\n");
    //
    printf("\nFor detailed background, see\n");
    printf("Larsen, E., Bisaillon, S., Cordeau, J.F. and Frejinger, E.\n");
    printf("Pseudo-random Instance Generators in C++ for Deterministic and Stochastic\n");
    printf("Multi-commodity Network Design Problems.\n");
    //
    printf("\nThe two following paragraphs are partly adapted from\n");
    printf("../HOYLAND_KAUT_WALLACE/sg_HKW.c:\n");
    //
    printf("\nScenario generation code is based on 'A Heuristic for Moment-matching\n");
    printf("Scenario Generation' by K. Høyland, M. Kaut & S.W. Wallace; Computational\n");
    printf("Optimization and Applications, 24 (2-3), pp. 169–185, 2003;\n");
    printf("doi:10.1023/A:1021853807313. Code by Michal Kaut (michal.kaut@iot.ntnu.no)\n");
    printf("& Diego Mathieu.\n");
    //
    printf("\nThe HKW code generates scenarios for multivariate random variables by\n");
    printf("matching specified first four moments and correlations. The algorithm\n");
    printf("is iterative and gradually closes the gap from target values. Although\n");
    printf("convergence has not been demonstrated theoretically, extensive empirical\n");
    printf("testing and application have shown that convergence occurs provided the\n");
    printf("requested number of scenarios is sufficiently large. Several trials may be\n");
    printf("required, each one starting from new random initializations automatically.\n");
    printf("(Cf. Section 6 below.)\n");
    //
    printf("\n1 Building and running the program directly from the Linux command line\n");
    printf("***********************************************************************\n");
    //
    printf("\n1.1 Building the program from the Linux command line\n");
    printf("****************************************************\n");
    //
    printf("\nUnder Linux, follow these steps to build the program:\n");
    //
    printf("\na- Install the CPLEX solver and/or the HIGHS solver (if not already done)\n");
    printf("and note the installation folder(s). According to our verifications, CPLEX\n");
    printf("and HIGHS work equally well in our application. HIGHS is freely available\n");
    printf("under the MIT license at https://highs.dev/. See Huangfu, Q., Hall, J.A.J.\n");
    printf("Parallelizing the dual revised simplex method. Math. Prog. Comp. 10, 119-142\n");
    printf("(2018).\n");
    printf("b- Either (i) use Git to clone locally the remote repository\n");
	printf("https://github.com/larseeri/stoch-MCFNDP-gener, or (ii) generate, download\n");
	printf("locally to an appropriate folder and uncompress a .zip archive containing\n");
	printf("the remote repository file tree.\n");	
    printf("c- Edit the file named 'makefile' to set the variables indicating the\n");
    printf("locations relevant to the CPLEX solver and/or the HIGHS solver, as noted in\n");
    printf("step a-. These locations are chiefly CPLEX_HOME and HIGHS_HOME. Save the\n");
    printf("modifications.\n");
    printf("d- In a Linux terminal, go to the /program folder just created.\n");
    printf("e- Make sure that the gcc C++ compiler and the make C++ build tool are\n");
    printf("available.\n");
    printf("f- In the terminal, run ONLY ONE of these commands to build the program\n");
    printf("equipped either with the CPLEX solver or with the HIGHS solver:\n");
    //
    printf("\nmake CPLEX_SOLVER=1\n");
    //
    printf("\nmake HIGHS_SOLVER=1\n");
    //
    printf("\nAn executable file named 'stoch_gen' will be created inside the program\n");
    printf("folder. Building with versions 11.5.0, 12.2.0, 13.2.0, 14.2.0 and 15.1.0 of\n");
    printf("the gcc compiler has been verified.\n");
    //
    printf("\nRemarks:\n");
    printf("- To build once again with HIGHS if the program has already been built with\n");
    printf("CPLEX, run the following:\n");
    //
    printf("\nmake clean CPLEX_SOLVER=1\n");
    printf("make HIGHS_SOLVER=1\n");
    //
    printf("\n- To build once again with CPLEX if the program has already been built with\n");
    printf("HIGHS, run the following:\n");
    //
    printf("\nmake clean HIGHS_SOLVER=1\n");
    printf("make CPLEX_SOLVER=1\n");
	//
    printf("\n1.2 Running the program with Linux command line instructions\n");
	printf("************************************************************\n");
	//
	printf("\nThe program can be run as follows from the /program folder just created\n");
	printf("along with default option parameter values:\n");
	//
	printf("\n./stoch_gen\n");
	//
	printf("\nOptions for running the program can be specified directly in the command\n");
	printf("line or in a configuration file. By default, the program will attempt to\n");
	printf("read option values from file stoch_config.par. (This default setting can be\n");
	printf("updated by modifying line 1197 of file main.cpp.) The default file can\n");
	printf("be overridden with option +F:\n");
	//
	printf("\n./stoch_gen +F newParamfile.par\n");
	//
	printf("\nThis assumes that stoch_config.par and newParamfile.par are located at the\n");
	printf("root of the /program folder. Otherwise, their absolute paths or their paths\n");
	printf("relative to the /program folder must be explicited.\n");
	//
	printf("\nThe default configuration file stoch_config.par has initially been filled\n");
	printf("with illustrative values and should be edited appropriately if no\n");
	printf("replacement configuration file is indicated on the command line with option\n");
	printf("key +F.\n");
	//
	printf("\nVerbosity option +v should be used to obtain more information about all\n");
	printf("options that have been set in the parameter file or in the command line:\n");
	//
	printf("\n./stoch_gen +v +F newParamfile.par\n");
	//
	printf("\nThe configuration file should have this format:\n");
	//
	printf("\noptionKey1   value1\n");
	printf("optionKey2   value2\n");
	printf(".\n");
	printf(".\n");
	printf("optionKeyN   valueN\n");
	//
	printf("\nIt may specify at most one option per line. Option keys may be prefixed or\n");
	printf("not with '-'. Comments may appear either at the beginning of a line, in\n");
	printf("which case they must start with a non-alphabetic character (but not with\n");
	printf("'-'), or at the end of a line.\n");
	//
	printf("\nThe general command line for running the program is as follows:\n");
	//
	printf("\n./stoch_gen [+v] [+F nameOfConfigFile] [main and HKW options]\n");
	//
	printf("\nOnce the configuration file has been processed, options appearing in the\n");
	printf("command line (other than +F and +v) are processed. Option values specified\n");
	printf("in the command line supersede those specified in the configuration file.\n");
	printf("That is, if included, they replace, modify and supplement the latter.\n");
	printf("Option keys appearing in the command line must be prefixed with '-' as\n");
	printf("follow:\n");
	//
	printf("\n-optionkeyA valueA -optionKeyB valueB...\n");
	//
	printf("\nBy default, boolean values are set to false and presence of option key is\n");
	printf("sufficient to set the corresponding value to true. For other options\n");
	printf("expecting string, integer or double values, values should be specified.\n");
	//
	printf("\nTo display the complete documentation about the application and the\n");
	printf("available options, run\n");
	//
	printf("\n./stoch_gen -H\n");
	//
	printf("\nor\n");
	//
	printf("\n./stoch_gen -h\n");
	//
	printf("\nSumming up, option values are successively updated from these locations in\n");
	printf("this order:\n");
	printf("(i) hard-coded initialization values in lines 27-51 of configuration.cpp\n");
	printf("(ii) values appearing in configuration file whose name is specified in\n");
	printf("command line with +F (if it is specified and if it exists); if no valid\n");
	printf("configuration file is specified with +F in command line, values will be\n");
	printf("taken from default configuration file whose name is specified in line 1197\n");
	printf("of main.cpp; at this point, the latter should exist\n");
	printf("(iii) values specified in the command line.\n");
	printf("Notice that location (i) is insufficiently transparent to the user and\n");
	printf("should not be relied upon. Parameter values should be set through (ii) and \n");
	printf("(iii).\n");
	//
	printf("\nRemark:\n");
	printf("See Section 8 below for examples of configuration files and instructions\n");
	printf("when the program is run with command line instructions.\n");
	//
	printf("\n2 Running the program through Docker instructions\n");
	printf("*************************************************\n");
	//
	printf("\nAn advantage of running the executable stoch_gen inside a Docker container\n");
	printf("is the possibility to avoid selection and installation of suitable C++\n");
	printf("compiler, C++ make build tool and required C++ libraries. It suffices that\n");
	printf("the Docker engine be installed on the host where the computation will be\n");
	printf("run. For this, see https://docs.docker.com/engine/install/.\n");
	//
	printf("\nOnce the Docker engine is installed, follow these steps:\n");
	//
	printf("\na- Open a terminal and retrieve from the GitHub container registry the\n");
	printf("Docker image ghcr.io/larseeri/stoch_gen_image:latest as follows:\n");
	//
	printf("\ndocker pull ghcr.io/larseeri/stoch_gen_image:latest\n");
	//
	printf("\n(From this image, Docker containers running the executable 'stoch_gen' will\n");
	printf("be generated.)\n");
	printf("b- For simplicity, define the alias 'stoch_gen_image' through\n");
	//
	printf("\ndocker tag ghcr.io/larseeri/stoch_gen_image:latest stoch_gen_image\n");
	//
	printf("\nRemark:\n");
	printf("If the Docker image stoch_gen_image were instead available as the tar ball\n");
	printf("stoch_gen_image.tar, then it could be loaded through\n");
	//
	printf("\ndocker load --input stoch_gen_image.tar\n");
	//
	printf("\nc- Either (i) use Git to clone locally the remote repository\n");
	printf("https://github.com/larseeri/stoch-MCFNDP-gener, or (ii) generate, download\n");
	printf("locally to an appropriate folder and uncompress a .zip archive containing\n");
	printf("the remote repository file tree. The /docker folder contains (i) the\n");
	printf("primitive files named 'Dockerfile' and 'makefile' based on which the\n");
	printf("Docker image was originally generated, (ii) a subfolder /inout used\n");
	printf("for exchanging files between the program running inside a Docker container\n");
	printf("and the directory structure of the host, (iii) some auxiliary files, and\n");
	printf("(iv) the subfolder /examples containing a number of configuration files and\n");
	printf("their resulting outputs (examined in Section 5 below).\n");
	printf("d- In a terminal, go to the /docker folder just created.\n");
	//
	printf("\nThen, instructions similar to those appearing in Section 1.2 can be run from\n");
	printf("the terminal on AMD and Intel CPUs (other CPUs are addressed further below).\n");
	printf("For example, the default Docker instruction\n");
	//
	printf("\ndocker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image /stoch_gen\n");
	//
	printf("\nis equivalent to the default Linux command line instruction\n");
	//
	printf("\n./stoch_gen\n");
	//
	printf("\nIn general, Docker instructions can be formulated from those appearing in\n");
	printf("Section 1.2 by (i) adding the prefix\n"); 
	//
	printf("\ndocker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image\n");
	//
	printf("\nand by (ii) replacing the executable location ./stoch_gen with /stoch_gen\n");
	printf("(i.e. removing the initial dot), since the executable stoch_gen is located\n");
	printf("at the root of the Docker container being generated and since locations\n"); 
	printf("inside a Docker container must be accessed through absolute paths.\n");
	//
	printf("\nHence, the following command line instructions:\n");
	//
	printf("\n./stoch_gen +v +F newParamfile.par\n");
	//
	printf("\n./stoch_gen -h\n");
	//
	printf("\n./stoch_gen +v +F newParamFile.par -r 1234 -R 4567\n");
	//
	printf("\nrespectively correspond to the following Docker instructions:\n");
	//
	printf("\ndocker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image \\\n");
	printf("/stoch_gen +v +F /inout/newParamfile.par\n");
	//
	printf("\ndocker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image \\\n");
	printf("/stoch_gen -h\n");
	//
	printf("\ndocker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image \\\n");
	printf("/stoch_gen +v +F /inout/newParamFile.par -r 1234 -R 4567\n");
	//
	printf("\nImportant:\n"); 
	printf("When running the program through Docker instructions, all input and output\n");
	printf("files must be located inside the docker/inout subfolder. Within the Docker\n"); 
	printf("instructions and the configuration files, this subfolder must be denoted\n");
	printf("/inout since it is located at the root of the Docker container. Hence, a\n");
	printf("replacement configuration file located at docker/inout/newParamFile.par\n");
	printf("will be specified in the Docker instruction with +F /inout/newParamFile.par.\n");
	printf("Similarly, a base network file at docker/inout/baseNetworkFile.std will be\n");
	printf("specified in the Docker instruction or inside the configuration file with\n");
	printf("-I /inout/baseNetworkFile.std.\n");		
	//	
	printf("\nThe default configuration file is located at docker/inout/stoch_config.par.\n");
	printf("It has initially been filled with illustrative values and should be edited\n");
	printf("appropriately if no replacement configuration file is indicated on the\n");
	printf("Docker instruction line with option key +F.\n");
	//
	printf("\nOption parameter values are successively updated from these locations in\n");
	printf("this order:\n");
	printf("(i) initialization values hard-coded in lines 27-51 of configuration.cpp,\n");
	printf("(ii) either replacement configuration file specified in the Docker\n");
	printf("instruction line with option key +F, or default configuration file\n");
	printf("./inout/stoch_config.par if no replacement configuration file is\n")
	;printf("specified in the Docker instruction line with option key +F,\n");
	printf("(iii) values specified in the Docker instruction.\n");
	printf("Notice that location (i) is insufficiently transparent to the user and\n");
	printf("should not be relied upon. Parameter values should be set through (ii) and\n");
	printf("(iii).\n");
	//
	printf("\nImportant:\n");
	printf("The supplied Docker image conforms to the amd64 CPU instructions standard\n");
	printf("and will run natively on AMD and Intel CPUs. On CPUs conforming instead to\n");
	printf("the arm64 standard, e.g., on recent Apple computers, the supplied Docker\n");
	printf("image can be run through emulation. Essentially, this may be achieved as\n");
	printf("follows:\n");
	printf("a- Running the following preliminary instruction once. (For details, see\n");
	printf("https://hub.docker.com/r/tonistiigi/binfmt.)\n");
	//
	printf("\ndocker run --privileged --rm tonistiigi/binfmt --install all\n");
	printf("b- Adding the option --platform linux/amd64 to the Docker instruction line\n");
	printf("prefix as follows:\n");
	//
	printf("\ndocker run --platform linux/amd64 --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image\n");
	//
	printf("\nThe following characteristics are specific to the operation of the program\n");
	printf("through Docker instructions:\n");
	printf("- An empty folder absolutePathToDockerFolder/inout must already exist on\n");
	printf("the host for the Docker instructions to work correctly.\n");
	printf("- Option --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout\n");
	printf("connects subfolder docker/inout on the host to folder /inout\n");
	printf("of the Docker container. Hence, the connection between\n");
	printf("absolutePathToDockerFolder/inout on the host and folder /inout in the\n");
	printf("container is used to supply input files to and retrieve output files from\n");
	printf("the container.\n");
	printf("- Option --rm ensures suppression of the container implicitly generated by\n");
	printf("docker run. Note that the container created by docker run for running the\n");
	printf("executable stoch_gen stops running immediately after the latter finishes\n");
	printf("and that it cannot be reused afterwards, for instance through a docker exec\n");
	printf("instruction. Immediate suppression with --rm avoids uselessly cluttering\n");
	printf("Docker's container repository.\n");
	//
	printf("\nThe following Docker instructions may be useful:\n");
	//
	printf("\ndocker image ls (list all currently existing docker images)\n");
	//
	printf("\ndocker image rm --force imageName (forcefully remove imageName)\n");
	//
	printf("\ndocker container ls (list all running containers)\n");
	//
	printf("\ndocker container ls --all (list all containers, running or not)\n");
	//
	printf("\nSee Section 8 below for examples of configuration files and instructions\n");
	printf("when the program is run through Docker.\n");
    //
    printf("\n3 Randomized elements, target moments and target correlations\n");
    printf("*************************************************************\n");
    //
    printf("\nThe generator rests on three primary components: a subset of randomized\n");
    printf("parameters, target moments and target correlations. We proceed to describe\n");
    printf("how each one is specified.\n");
    //
    printf("\nRemarks:\n");
    printf("- Section 7 below describes detailed usage of every option. Section 8\n");
    printf("presents systematically some examples of configuration options specified in\n");
    printf("configuration files.\n");
    printf("- For the purposes of the current discussion, <...> indicates that some\n");
    printf("value or values of the type specified inside the brackets is expected. For\n");
    printf("instance, <string> means that a string is expected whereas <> means that no\n");
    printf("value is expected.\n");
    //
    printf("\n3.1 Base deterministic MCFND instance and randomized elements\n");
    printf("*************************************************************\n");
    //
    printf("\nOption S <int> must always be specified, either in the instruction line or\n");
    printf("in the configuration file. It indicates which subset of parameters of a\n");
    printf("given base deterministic MCFND instance will be assigned a joint probability\n");
    printf("distribution and will thus vary between scenarios. We call these parameters\n");
    printf("randomized elements. Option S expects an integer from 1 to 31 that is the\n");
    printf("sum of the following, where symbol -> is a shorthand notation for 'calls\n");
    printf("for randomization of':\n");
    printf("1 -> demand volumes\n");
    printf("2 -> total capacities of arcs\n");
    printf("4 -> commodity-specific capacities of arcs\n");
    printf("8 -> fixed costs\n");
    printf("16 -> variable costs (either commodity-specific or not)\n");
    printf("Hence, specifying -S 3 in the instruction line indicates that total\n");
    printf("capacities of arcs (2) and demand volumes (1) must be randomized\n");
    printf("(2 + 1 = 3).\n");
    //
    printf("\nBy default, the name of the file describing the base deterministic MCFND\n");
    printf("network is './base_network.dow' (Linux instructions),\n");
    printf("'./inout/base_network.dow' (Docker instructions) and the default input\n");
    printf("format is DOW. This name may be changed with option I <string> and the input\n");
    printf("format can be indicated accordingly with option F <char>. For example,\n");
    printf("specifying -S 3 -I instB.std -F G in the instruction line indicates that\n");
    printf("demands and total capacities of arcs will be randomized based on the\n");
    printf("deterministic network stored in file instB.std using the generic STD format.\n");
    printf("(Cf. Section 5.1 below for details.)\n");
    //
    printf("\nImportant:\n");
    printf("It is critical that the content of the file describing the base\n");
    printf("deterministic MCFND instance and the indicated format of its content match\n");
    printf("together.\n");
    //
    printf("\nRemark:\n");
    printf("While the generator is sufficiently general to randomize fixed costs, the\n");
    printf("latter are viewed as being non-stochastic in the context of two-stage\n");
    printf("stochastic programming. Hence, in standard applications, randomization of\n");
    printf("the MCFND will be limited to a subset of the following: demand volumes,\n");
    printf("total capacities of arcs, commodity-specific capacities of arcs and variable\n");
    printf("costs (either commodity-specific or not).\n");
    //
    printf("\n3.2 Target moments\n");
    printf("******************\n");
    //
    printf("\nTarget moments can be set and saved as follows.\n");
    printf("a- If option key G is included, target moments will be generated according\n");
    printf("to a distributional assumption and will be written individually in a file\n");
    printf("whose default name is './target_moms.dat' (Linux instructions) or\n");
    printf("'./inout/target_moms.dat' (Docker instructions). This name can be specified\n");
    printf("otherwise with option key M. The distributional assumption can be specified\n");
    printf("with option key T <char> where <char> stands for either U (uniform) or D\n");
    printf("(triangular) and two parameters, alpha and beta, that can be specified using\n");
    printf("options A <double> and B <double>.\n");
    //
    printf("\nb- If option key G is not included, individual target moments will be read\n");
    printf("from the file whose default name is is './target_moms.dat' (Linux\n");
    printf("instructions) or './inout/target_moms.dat' (Docker instructions). The latter\n");
    printf("can be specified otherwise with option key M. Format of the target moments\n");
    printf("file is specified in Section 5.2.1 below. Without option G, option T will be\n");
    printf("rejected and options A and B will be ignored.\n");
    //
    printf("\nImportant:\n");
    printf("In case b-, that is when option G is omitted and target moments are read\n");
    printf("from a file, it is crucial that they are compatible with a proper\n");
    printf("probability distribution. If this condition fails, an attempt at\n");
    printf("randomization is meaningless and convergence of the HKW algorithm is\n");
    printf("expected to fail. In case a-, that is when option G is included, the target\n");
    printf("moments that are generated are by construction compatible with a proper\n");
    printf("probability distribution.\n");
    //
    printf("\n3.2 Target correlations\n");
    printf("***********************\n");
    //
    printf("\nTarget correlations can be set and saved as follows.\n");
    printf("a- If option key K is included, target correlations can be directly\n");
    printf("specified in blocks either (i) in the instruction line, or (ii) in the\n");
    printf("configuration file. A common target correlation coefficient within subset of\n");
    printf("parameters U or between subsets of parameters U and V is specified with\n");
    printf("X<UU> <double> and X<UV> <double> respectively, where <double> is the target\n");
    printf("correlation coefficient and U and V, can take these values, where symbol ->\n");
    printf("is a shorthand notation for 'stands for':\n");
    printf("D -> demands\n");
    printf("A -> total capacities of arcs\n");
    printf("C -> commodity-specific capacities of arcs\n");
    printf("F -> fixed costs\n");
    printf("V -> variable costs\n");
    printf("Target correlations that are left unspecified are assumed by default to be\n");
    printf("zero. The resulting individual correlations will be written to a file whose\n");
    printf("default name is './target_corrs.dat' (Linux instructions) or\n");
    printf("'./inout/target_corrs.dat' (Docker instructions). This name can be\n");
    printf("specified otherwise with option key C.\n");
    //
    printf("\nImportant:\n");
    printf("Cases (i) and (ii) are not coexistent: if any correlation information is\n");
    printf("supplied in the instruction line, any correlation information supplied in\n");
    printf("the configuration file is ignored.\n");
    //
    printf("\nb- If option key K is not included, individual correlations will be read\n");
    printf("from a file whose default name is './target_corrs.dat' (Linux instructions)\n");
    printf("or './inout/target_corrs.dat' (Docker instructions). The latter can be\n");
    printf("specified otherwise with option key C. Format of the target correlations\n");
    printf("file is specified in Section 5.2.2 below. Without option K, options X?? will\n");
    printf("be rejected.\n");
    //
    printf("\nImportant:\n");
    printf("In case a-, when option key K is included, it is crucial that the block\n");
    printf("correlations that are specified result in a positive definite correlation\n");
    printf("matrix. Similarly, in case b-, when option key K is omitted, the supplied\n");
    printf("matrix must define a positive definite correlation matrix. Whenever the\n");
    printf("alleged correlation matrix, either constructed (case a-) or read from file\n");
    printf("(case b-), fails positive definitess, its Cholesky decomposition will fail\n");
    printf("in the HKW algorithm. A suggestion to verify and correct the alleged\n");
    printf("correlation matrix will be displayed and the program will be interrupted.\n");
    //
    printf("\n3.4 Example\n");
    printf("***********\n");
    //
    printf("\nConsider this hypothetical Linux command line instruction:\n");
    //
    printf("\n./stoch_gen -S 3 -I instB.std -F G -K -XDD 0.5 -XDA -0.3 -XAA 0.7 -G -T U -A 0.25 -B 0.25\n");
    //
    printf("It is interpreted as follows:\n");
    printf("-S 3: demands and total capacities of arcs will be randomized\n");
    printf("-I instB.std: base deterministic network is to be read from file instB.std\n");
    printf("-F G: format of the latter is generic STD\n");
    printf("-K: target correlations are specified in blocks (based on explicit and\n");
    printf("default values)\n");
    printf("-XDD 0.5: target correlations among demands are 0.5 for all dissimilar\n");
    printf("pairs and self-correlations are 1.0\n");
    printf("-XDA -0.3: target correlations between demands and total capacities of\n");
    printf("arcs are -0.3 for all pairs\n");
    printf("-XAA 0.7: target correlations among total capacities of arcs are 0.7\n");
    printf("for all dissimilar pairs and self-correlations are 1.0\n");
    printf("-G: target moments are generated (based on explicit and default\n");
    printf("specifications)\n");
    printf("-T U: target moments are generated under the assumption that\n");
    printf("randomizations each obey a uniform distribution\n");
    printf("-A 0.25: parameter alpha defining the latter distribution is 0.25\n");
    printf("-B 0.25: parameter beta defining the latter distribution is 0.25\n");
    //
    printf("\n4 Output file\n");
    printf("*************\n");
    //
    printf("\nBy default, the name of the file containing the MCFND instances generated\n");
    printf("by the application is './generated_networks.txt' (Linux instructions) or\n");
    printf("'./inout/generated_networks.txt' (Docker instructions). This can be\n");
    printf("changed with option O. This text file superposes representations of complete\n");
    printf("deterministic MCFND instances, one for each scenario realization. Each\n");
    printf("representation is preceded by a separator as follows:\n");
    //
    printf("\nSCENARIO   <int>\n");
    //
    printf("\nImmediately after the separator, the generated MCFND instance corresponding\n");
    printf("to the current scenario is described in its entirety with a format (DOW,\n");
    printf("restricted STD or generic STD; cf. Section 5.1 below), matching that used to\n");
    printf("supply the base deterministic MCFND instance. Stacking and separating full\n");
    printf("descriptions of the MCFND realizations in this fashion makes it easy to\n");
    printf("pinpoint and read the first and second stage information relevant to every\n");
    printf("individual scenario.\n");
    //
    printf("\nRemark:\n");
    printf("Feasibility of the second stage of each generated MCFND instance is verified\n");
    printf("with the selected solver (CPLEX or HIGHS) after opening all arcs. Only the\n");
    printf("generated networks that are demonstrated feasible are added to the output\n");
    printf("file. The application produces a report indicating how many scenarios were\n");
    printf("generated and how many among the latter are feasible.\n");
    //
    printf("\n5 Formats of data files\n");
    printf("***********************\n");
    //
    printf("\n5.1 Base deterministic MCFND instance and generated MCFND instances\n");
    printf("*******************************************************************\n");
    //
    printf("\nDOW and STD formats are commonly used to store deterministic MCFND instances\n");
    printf("in text files. These formats are described in files DOW-format-desc.txt and\n");
    printf("STD-format-desc.txt. The application can read from these two formats and, in\n");
    printf("addition, from a restricted version of STD. Thus, three formats are\n");
    printf("available for supplying the base deterministic MCFND instance: DOW, STD\n");
    printf("(calling it generic STD, or GSTD) and a restricted version of STD (calling\n");
    printf("it restricted STD, or RSTD). One of these three formats may be selected with\n");
    printf("option F <char>.\n");
    printf("a- DOW presupposes that (i) for each arc, capacities and variable\n");
    printf("costs are identical for all commodities and (ii) each commodity has a\n");
    printf("single source node (at which volume > 0) and a single sink node (at which\n");
    printf("volume < 0).\n");
    printf("b- Generic STD imposes neither (i), nor (ii).\n");
    printf("c- Restricted STD imposes only (ii).\n");
    //
    printf("\nImportant:\n");
    printf("When the base deterministic MCFND instance is received in the generic STD\n");
    printf("format, there are potentially multiple source and/or sink nodes for some or\n");
    printf("all commodities. In this case, whenever randomization of commodity volumes\n");
    printf("is requested through option S, it is handled through these steps:\n");
    printf("a- For each commodity, the total volume, the flows out of each source node\n");
    printf("and the flows into each sink node occurring in the base deterministic MCFND\n");
    printf("network are saved for ulterior computations while reading (see\n");
    printf("InstanceGStd::read() in generate.cpp).\n");
    printf("b- Total commodity volumes are randomized (along with other stochastic\n");
    printf("elements if there are any) as they are when there are single source and sink\n");
    printf("nodes for each commodity.\n");
    printf("c- For each commodity, the randomized flows out of each source node and\n");
    printf("randomized flows into each sink node are calculated by distributing the\n");
    printf("randomized total commodity volume using the original ratios of inflows and\n");
    printf("outflows over total volume obtained in step a-.\n");
    //
    printf("\nRemark:\n");
    printf("The MCFND instances generated by the application (one for each scenario) are\n");
    printf("expressed in the same format as that used for supplying the base\n");
    printf("deterministic MCFND instance. Their descriptions are stacked and returned in\n");
    printf("the output file whose default name './generated_networks.txt' (Linux\n");
    printf("instructions), './inout/generated_networks.txt' (Docker instructions) can be\n");
    printf("modified with option O. (Cf. Section 4 above.)\n");
    //
    printf("\n5.2 Input files supplied to HKW algorithm\n");
    printf("*****************************************\n");
    printf("(partly adapted from ../HOYLAND_KAUT_WALLACE/sg_HKW.c)\n");
    //
    printf("\n5.2.1 Target moments\n");
    printf("********************\n");
    //
    printf("\nWhen option G is omitted from instruction line and from configuration file,\n");
    printf("target moments are read from a file named './target_moms.dat' (Linux\n");
    printf("instructions) or './inout/target_moms.dat' (Docker instructions) unless this\n");
    printf("is overridden with option M. The file containing the target moments matrix\n");
    printf("is formatted as follows (see Section 8 for examples):\n");
    printf("a- 1st line shows integer 4 (as there are 4 target moments per randomized\n");
    printf("stochastic element and therefore 4 rows in the matrix of target moments),\n");
    printf("b- 2nd line states the total number of individual randomized stochastic\n");
    printf("elements,\n");
    printf("c- 3rd line states the first target moment of each randomized stochastic\n");
    printf("element, separated by spaces and listed in this order: demands (if\n");
    printf("randomized), total capacities of arcs (if randomized), commodity-specific\n");
    printf("capacities of arcs(if randomized), fixed costs (if randomized), variable\n");
    printf("costs (if randomized). The number of figures in this line is equal to the\n");
    printf("number specified in 2nd line.\n");
    printf("d- 4th to 6th lines respectively state second, third and fourth target\n");
    printf("moments, exactly as 3rd row does for first target moment.\n");
    printf("e- Remaining lines of the file are ignored.\n");
    //
    printf("\nImportant:\n");
    printf("Exact definitions of the four target moments must satisfy these\n");
    printf("conditions:\n");
    printf("a- Population estimators are used (as in spreadsheets).\n");
    printf("b- 2nd moment is defined as variance, that is second central\n");
    printf("statistical moment, and not standard deviation.\n");
    printf("c- 3rd moment is defined as third central statistical moment divided\n");
    printf("by third power of standard deviation (i.e., standardized skewness).\n");
    printf("d- 4th moment is defined as fourth central statistical moment divided\n");
    printf("by fourth power of standard deviation minus 3 (i.e., excess standardized\n");
    printf("kurtosis).\n");
    //
    printf("\n5.2.2 Target correlations\n");
    printf("*************************\n");
    //
    printf("\nWhen option K is omitted from instruction line and from configuration file,\n");
    printf("target correlations are read from a file named './target_corrs.dat' (Linux\n");
    printf("instructions) or './inout/target_corrs.dat' (Docker instructions), unless\n");
    printf("this is overridden with option C. The file containing the target\n");
    printf("correlations matrix is formatted as follows (see Section 8 for examples):\n");
    printf("a- 1st and 2nd lines are identical and show the total number of individual\n");
    printf("randomized elements.\n");
    printf("b- Remaining lines of the file contain the target correlations matrix, one\n");
    printf("row of the matrix per line of the file. Individual correlations are\n");
    printf("separated by spaces.\n");
    printf("c- The correlation pairs appearing in the matrix are ordered by row and by\n");
    printf("column. Order in a row and order in a column is similar to that in a row\n");
    printf("of the matrix of target moments (cf. Section 4.1.1).\n");
    printf("d- Any additional line in the file is ignored.\n");
    //
    printf("\n5.2.3 Probabilities of scenarios\n");
    printf("********************************\n");
    //
    printf("\nIf option p is omitted, then scenarios are assumed to be equiprobable.\n");
    printf("Option p can be used to supply the name of a file containing probabilities\n");
    printf("of all scenarios. The following format must be satisfied (see supplied\n");
    printf("example probs_100_eq_example.dat):\n");
    printf("a- 1st line states the number of scenarios\n");
    printf("b- 2nd line states all probabilities, separated by spaces\n");
    printf("c- Remaining lines of the file are ignored.\n");
    printf("This file would be supplied as './probs_100_eq_example.dat' (Linux\n");
    printf("instructions) or './inout/probs_100_eq_example.dat' (Docker instructions).\n");
    //
    printf("\n5.2.4 Raw scenarios in HKW format\n");
    printf("*********************************\n");
    //
    printf("\nRaw scenarios are saved in HKW format to './raw_scens.dat' (Linux\n");
    printf("instructions) or './inout/raw_scens.dat' (Docker instructions), unless this\n");
    printf("is overriden with option o. Raw scenarios saved in previous run can be used\n");
    printf("as starting values in a following run (using option s). These will only be\n");
    printf("used in the first trial of the current run of the HKW algorithm. If this\n");
    printf("trial does not converge and the maximum number of trials specified by option\n");
    printf("t is larger than 1, the next trial starts from random guesses.\n");
    //
    printf("\nImportant:\n");
    printf("Raw scenarios pertain to the individual stochastic elements. They are\n");
    printf("distinct from the MCFND instance scenarios produced by the application and\n");
    printf("described in Section 4 above.\n");
    //
    printf("\n6 Convergence\n");
    printf("*************\n");
    //
    printf("\nRemark:\n");
    printf("The level of on-screen reporting by the HKW algorithm should be set at\n");
    printf("no less than the default value of 3 (option l) to adequately follow the\n");
    printf("convergence process.\n");
    //
    printf("\nThe HKW algorithm might fail to satisfy the specified target moments and\n");
    printf("correlations within the default tolerances with the requested number of\n");
    printf("scenarios within the default number of trials and iterations. If\n");
    printf("convergence fails, the following course of action is recommended:\n");
    //
    printf("\n1- Verify that target moments and correlations are proper.\n");
    //
    printf("If option G has not been selected, target moments are read from a file. It\n");
    printf("should be ascertained they are compatible with a probability distribution.\n");
    printf("If this condition fails, convergence of the HKW algorithm is expected to\n");
    printf("fail. By construction, target moments that are generated when option G is\n");
    printf("selected are correctly specified.\n");
    //
    printf("\nIf option K has not been selected, a target correlation matrix is read from\n");
    printf("a file. it should be ascertained this matrix is well-specified (square,\n");
    printf("symmetric, all elements in [-1,1], diagonal elements equal to 1, positive\n");
    printf("definite). If this condition fails, the HKW algorithm is expected to fail.\n");
    //
    printf("\nIf option K has been selected, a target correlation matrix is built from\n");
    printf("block correlations specified with options X??. it should be ascertained the\n");
    printf("resulting alleged target correlation matrix is well-specified.\n");
    //
    printf("\nRemarks:\n");
    printf("- If the alleged correlation matrix fails positive definiteness, the HKW\n");
    printf("algorithm initially (i.e. before any iteration) reports a failure of its\n");
    printf("Cholesky decomposition and displays an explanation.\n");
    printf("- However, failure of the Cholesky decomposition during the iterations of\n");
    printf("the HKW algorithm might be due to an insufficient requested number of\n");
    printf("scenarios (see Paragraph 3 below).\n");
    //
    printf("\n2- Verify that slow convergence is not responsible.\n");
    //
    printf("\nIf both target moments and target correlations are known to be properly\n");
    printf("specified, consider launching a new run after increasing the maximum numbers\n");
    printf("of trials and iterations, respectively options t and i. Doubling their\n");
    printf("values from default is suggested.\n");
    //
    printf("\nRemark:\n");
    printf("Default values for the numbers of trials and iterations are usually\n");
    printf("sufficient to reach convergence under the default error tolerances. Failure\n");
    printf("of convergence is most frequently due to the causes examined in Paragraphs 1\n");
    printf("and 3 above and below.\n");
    //
    printf("\n3- Adjust the number of scenarios.\n");
    //
    printf("If both target moments and target correlations are properly specified and\n");
    printf("slow convergence is excluded, then the requested number of scenarios is too\n");
    printf("small to support the specified target moments and target correlations.\n");
    printf("Hence, the requested number of scenarios, option n, may be progressively\n");
    printf("increased by a moderate factor (e.g., between 1.5 and 5) until the algorithm\n");
    printf("converges. Once it has converged, the number of requested scenarios may be\n");
    printf("progressively adjusted downward, if desired, while the algorithm is still\n");
    printf("converging.\n");
    //
    printf("\nRemark:\n");
    printf("The lower bound on the number of requested scenarios leading to convergence\n");
    printf("is positively related to the following: number of random elements, magnitude\n");
    printf("of 3rd and 4th target moments, heterogeneity of target correlations,\n");
    printf("heterogeneity of probabilities.\n");

    printUsageAndExit(0, true);
    exit(1);
}

void Configuration::printUsageAndExit(char ExecName[], bool anonymous)
{

    if (!anonymous)
    {
        printf("\nUsage: %s [main options] [HKW options]\n", ExecName);
    }

    printf("\n7 Usage of options\n");
    printf("******************\n");
    //
    printf("\n(Some parts are adapted from ../HOYLAND_KAUT_WALLACE/sg_HKW.c.)\n");
    //
    printf("\nExcept for option key S (and leaving aside option key F when it is prefixed\n");
    printf("by '+' in the instruction line), non-boolean user-selected options that are\n");
    printf("specified neither in the instruction line nor the configuration file are set\n");
    printf("to the values that are hard-coded in the constructor in lines 27-51 of\n");
    printf("configuration.cpp. Boolean user-selected options that are specified neither\n");
    printf("in the instruction line nor in the configuration file are set to false by\n");
    printf("default.\n");
    //
    printf("\nNotation:\n");
    printf("For the purposes of this presentation, <...> indicates that some value or\n");
    printf("values of the type specified inside the brackets is expected. For instance,\n");
    printf("<string> means that a string is expected whereas <> means that no value is\n");
    printf("expected.\n");
    //
    printf("\n7.1 Options available in instruction line only (Linux or Docker)\n");
    printf("****************************************************************\n");
    //
    printf("\n+F <string>; name of configuration file supplying option parameter values;\n");
    printf("this replaces the default file whose name is hard-coded in line 1197 of\n");
    printf("main.cpp; any option value supplied in instruction line replaces value of\n");
    printf("same option supplied in configuration file\n");
    //
    printf("\n+v <> (no value supplied); requests verbosity while parsing options in\n");
    printf("instruction line or configuration file\n");
    //
    printf("\nRemark:\n");
    printf("Unknown option keys prefixed with '+' in the instruction line are\n");
    printf("disregarded.\n");
    //
    printf("\n7.2 Options available in instruction line (Linux or Docker) and in configuration file\n");
    printf("*************************************************************************************\n");
    //
    printf("\nThe following options must be prefixed by '-' when they are called from the\n");
    printf("instruction line. This prefix can be omitted when the options are called\n");
    printf("from a configuration file. To supply an option, a line of a configuration\n");
    printf("file must begin with a dash '-' or a letter. A line beginning with any other\n");
    printf("character will be considered a comment. The remainder of a line beyond the\n");
    printf("option key and its value (if any is expected) is disregarded and may\n");
    printf("therefore be used for adding comments.\n");
    //
    printf("\nRemark:\n");
    printf("Unknown option keys following '-' in the instruction line, stray characters\n");
    printf("or strings appearing without the prefix '-' in the instruction line and\n");
    printf("unknown option keys (following '-' or not) at beginning of line in a\n");
    printf("configuration file are disregarded.\n");
    //
    printf("\n7.2.1 Main generator options\n");
    printf("****************************\n");
    //
    printf("\nH <> (no value supplied); display complete documentation about the\n");
    printf("application and the available options\n");
    //
    printf("\nh <> (no value supplied); display complete documentation about the\n");
    printf("application and the available options\n");
    //
    printf("\nF <char>; format of file where deterministic MCFND base network is to be\n");
    printf("read; D: DOW, G: generic STD, R: restricted STD; default: D\n");
    //
    printf("\nI <string>; name of file where deterministic MCFND base network is to be\n");
    printf("read; default value is './base_network.dow' (Linux instructions) or\n");
    printf("'./inout/base_network.dow' (Docker instructions)\n");
    //
    printf("\nO <string>; name of file where MCFND instances resulting from\n");
    printf("scenario generation are to be written; default value is\n");
    printf("'./generated_networks.txt' (Linux instructions) or\n");
    printf("'./inout/generated_networks.txt' (Docker instructions)\n");
    //
    printf("\nS <int>; identifies which subsets of parameters should vary between\n");
    printf("scenarios; expects a number (from 1 to 31) which is a sum of the following:\n");
    printf("1 -> demands\n");
    printf("2 -> total capacities of arcs\n");
    printf("4 -> commodity-specific capacities of arcs\n");
    printf("8 -> fixed costs\n");
    printf("16 -> variable costs\n");
    //
    printf("\nK <> (no value supplied); if included, target correlations are directly\n");
    printf("specified in blocks using options X and written to file\n");
    printf("'./target_corrs.dat' (Linux instructions) or './inout/target_corrs.dat'\n");
    printf("(Docker instructions), unless this default name is changed with option C;\n");
    printf("if omitted, target correlations are read from file './target_corrs.dat'\n");
    printf("(Linux instructions) or './inout/target_corrs.dat' (Docker instructions),\n");
    printf("unless this default name is changed with option C\n");
    //
    printf("\nC <string>; name of file where target correlations are read if K is omitted\n");
    printf("and where generated target correlations are written if K is included;\n");
    printf("default: './target_corrs.dat' (Linux instructions) or\n");
    printf("'./inout/target_corrs.dat' (Docker instructions)\n");
    //
    printf("\nX<char><char> <double>; identifies target correlations within or between\n");
    printf("subsets of parameters that are randomized (i.e, that vary between\n");
    printf("scenarios); state once for each pair of subsets; target self-correlations\n");
    printf("are equal to 1.0; any use of X without G will be rejected;\n");
    printf("default: zero\n");
    printf("D -> demands\n");
    printf("A -> total capacities of arcs\n");
    printf("C -> commodity-specific capacities of arcs\n");
    printf("F -> fixed costs\n");
    printf("V -> variable costs\n");
    printf("Example: in the instruction line, -XDD 0.9 -XFD 0.1 -XFF 0.01 specifies\n");
    printf("target correlations between distinct demand parameters of 0.9\n");
    printf("(self-correlations being equal to 1.0), between fixed cost parameters and\n");
    printf("demand parameters of 0.1 and between distinct fixed cost parameters of 0.01\n");
    printf("(self-correlations being equal to 1.0); Remark: -XFD and -XDF are equivalent\n");
    //
    printf("\nG <> (no value supplied); if included, target moments are generated based on\n");
    printf("option T and written to file './target_moms.dat' (Linux instructions) or\n");
    printf("'./inout/target_moms.dat' (Docker instructions), unless this default name is\n");
    printf("changed with option M; if omitted, target moments are read from file\n");
    printf("'./target_moms.dat' (Linux instructions) or './inout/target_moms.dat'\n");
    printf("(Docker instructions) unless this default name is changed with option M\n");
    //
    printf("\nM <string>; name of file where target moments are read if option G is\n");
    printf("omitted or where generated target moments are written if option G is\n");
    printf("included; default: './target_moms.dat' (Linux\n");
    printf("instructions) or './inout/target_moms.dat' (Docker instructions)\n");
    //
    printf("\nT <char>; distributional characteristics of generated target moments\n");
    printf("(U : UNIFORM , T : TRIANGULAR); any use of option T without option G\n");
    printf("will be rejected; default: U\n");
    printf("A <double>; parameter alpha used in generating target moments (see\n");
    printf("discussion immediately below); without option G, any use of option A will be\n");
    printf("ignored; default: 0.25\n");
    //
    printf("\nB <double>; parameter beta used in generating target moments (see\n");
    printf("discussion immediately below); without option G, any use of option B will be\n");
    printf("ignored; default: 0.25\n");
    //
    printf("\nWhen option T is set to U, then the stochastic elements of the network that\n");
    printf("are changed according to scenarios are each assumed to follow a\n");
    printf("uniform(a, b) distribution, where\n");
    printf("a = D - (alpha * D),\n");
    printf("b = D + (beta * D),\n");
    printf("D is the base value of the stochastic element taken from the base\n");
    printf("deterministic network,\n");
    printf("alpha must be in [0, 0.99],\n");
    printf("beta >= 0.\n");
    printf("The distribution is symmetric around D when alpha = beta.\n");
    //
    printf("\nWhen option T is set to T, then the stochastic elements of the network that\n");
    printf("are changed according to scenarios are each assumed to follow a\n");
    printf("triangular(a, b, c) distribution, where\n");
    printf("a = D - (alpha * D),\n");
    printf("b = D + (beta * D),\n");
    printf("c = D,\n");
    printf("D is the base value of the stochastic element taken from the base\n");
    printf("deterministic network,\n");
    printf("alpha must be in [0, 0.99],\n");
    printf("beta >= 0.\n");
    printf("The distribution is symmetric around D when alpha = beta.\n");
    //
    printf("\n7.2.2 Options supplied to HKW algorithm\n");
    printf("***************************************\n");
    //
    printf("\nScenario generation code is based on paper 'A Heuristic for Moment-\n");
    printf("matching Scenario Generation' by K. Høyland, M. Kaut & S.W. Wallace;\n");
    printf("Computational Optimization and Applications, 24 (2-3), pp. 169–185, 2003;\n");
    printf("doi:10.1023/A:1021853807313. Code by Michal Kaut (michal.kaut@iot.ntnu.no)\n");
    printf("& Diego Mathieu. Additional details about HKW options available in\n");
    printf("../HOYLAND_KAUT_WALLACE/sg_HKW.c.\n");
    //
    printf("\nn <int>; number of scenarios to generate; default: 1000\n");
    //
    printf("\nt <int>; maximum number of trials (attempts to generate scenarios using\n");
    printf("alternative random starting values); default: 20\n");
    //
    printf("\ni <int>; maximum number of iterations in a trial; default: 50\n");
    //
    printf("\nm <double>; maximum error in matching moments (scaled to var=1);\n");
    printf("default: 0.001\n");
    //
    printf("\nc <double>; maximum error in matching correlations (scaled to var=1);\n");
    printf("default: 0.001\n");
    //
    printf("\nl <int>; level of on-screen reporting by HKW algorithm (between 0 and 11);\n");
    printf("default: 3\n");
    //
    printf("\no <string>; name of file where resulting matrix of raw scenarios in HKW\n");
    printf("format is written; default: './raw_scens.dat' (Linux instructions) or\n");
    printf("'./inout/raw_scens.dat' (Docker instructions)\n");
    //
    printf("\ns <string>; if included, name of file where a matrix of starting scenarios\n");
    printf("saved in HKW format in a previous run (possibly using option o) is read and\n");
    printf("used for starting values; no default\n");
    //
    printf("\np <string>; if included, read probabilities from this file; otherwise,\n");
    printf("scenarios are assumed equiprobable; no default\n");
    //
    printf("\nr <int>; random seed; default: 7654321\n");
    //
    printf("\nR <int>; random stream; default: 1000\n");
    //
    printf("\n8 Examples\n");
    printf("**********\n");
    //
    printf("\nThe /program/examples subfolder holds examples of configuration files and\n"); 
	printf("instruction lines occurring when the executable of the program is run\n");
	printf("directly from the Linux instruction line. It also contains the output files\n");
	printf("resulting from each pair of configuration file and instruction line.\n");
	//
	printf("\nSimilarly, the /docker/examples subfolder holds examples of configuration\n");
	printf("files and instruction lines occurring when the executable of the program is\n");
	printf("run through Docker instructions. These examples are identified by the\n");
	printf("addition of the prefix 'docker_'. The subfolder also contains the output\n");
	printf("files resulting from each pair of configuration file and instruction line.\n");
	printf("The former are identical to the files generated when the executable is run\n");
	printf("from the Linux instruction line.\n");
    //
    printf("\n8.1 One source and one sink for each commodity\n");
    printf("**********************************************\n");
    //
    printf("\n8.1.1 Folder ./examples/oneSrcOneSnk/DOW/\n");
    printf("*****************************************\n");
    //
    printf("\nIn the base network, there is one source and one sink for each commodity and\n");
    printf("no differentiation of variable costs and capacities according to\n");
    printf("commodities. The base network is received in the DOW format and the\n");
    printf("generated networks are saved in the DOW format.\n");
    //
    printf("\n8.1.1.1 Operation of the stochastic generator\n");
    printf("*********************************************\n");
    //
    printf("\nLinux command line instruction:\n");
	//
    printf("\n./stoch_gen +v +F stoch_config_oneSrcOneSnk_10oct2025-18h10_dow.par\n");
    //
    printf("\nConfiguration file (input) used with Linux command line:\n");
    printf("stoch_config_oneSrcOneSnk_10oct2025-18h10_dow.par\n");
    //
    printf("\nDocker instruction:\n");
	//
    printf("\ndocker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image \\\n");
	printf("/stoch_gen +v +F docker_stoch_config_oneSrcOneSnk_10oct2025-18h10_dow.par\n");
    //
    printf("\nConfiguration file (input) used with Docker instruction:\n");
    printf("docker_stoch_config_oneSrcOneSnk_10oct2025-18h10_dow.par\n");
    //
    printf("\nBase network file (input):\n");
    printf("determ_oneSrcOneSnk_10oct2025-18h10.dow\n");
    //
    printf("\nGenerated moments (output):\n");
    printf("stoch_moms_oneSrcOneSnk_10oct2025-18h10_dow.dat\n");
    //
    printf("\nGenerated correlations (output):\n");
    printf("stoch_corrs_oneSrcOneSnk_10oct2025-18h10_dow.dat\n");
    //
    printf("\nRaw HKW scenarios (output):\n");
    printf("stoch_rawscen_oneSrcOneSnk_10oct2025-18h10_dow.dat\n");
    //
    printf("\nStacked scenarios of the MCFND instance in DOW fomat (output):\n");
    printf("stoch_oneSrcOneSnk_10oct2025-18h10_dow.txt\n");
    //
    printf("\n8.1.1.2 Origin of the base deterministic network\n");
    printf("************************************************\n");
    //
    printf("\nThe base deterministic network was originally produced by the determistic\n");
    printf("generator as follows.\n");
    //
    printf("\nLinux command line instruction:\n");
	//
    printf("\n./determ_gen +v +F determ_config_oneSrcOneSnk_10oct2025-18h10.par\n");
    //
    printf("\nConfiguration file:\n");
    printf("determ_config_oneSrcOneSnk_10oct2025-18h10.par\n");
    //
    printf("\nOutputs in DOW, STD, LP, MPS, GRA (basic graph), TEX (code of graphical\n");
    printf("representation), pdf (compiled image) formats:\n");
    printf("determ_oneSrcOneSnk_10oct2025-18h10.dow\n");
    printf("determ_oneSrcOneSnk_10oct2025-18h10.std\n");
    printf("determ_oneSrcOneSnk_10oct2025-18h10.mps\n");
    printf("determ_oneSrcOneSnk_10oct2025-18h10.lp\n");
    printf("determ_oneSrcOneSnk_10oct2025-18h10.gra\n");
    printf("determ_oneSrcOneSnk_10oct2025-18h10.tex\n");
    printf("determ_oneSrcOneSnk_10oct2025-18h10.pdf\n");
    //
    printf("\n8.1.2 Folder ./examples/oneSrcOneSnk/RSTD/\n");
    printf("******************************************\n");
    //
    printf("\nSimilar to 8.1.1, except that the base network is received in the STD format\n");
    printf("and the generated networks are saved in the RSTD format.\n");
    //
    printf("\n8.1.2.1 Operation of the stochastic generator\n");
    printf("*********************************************\n");
    //
    printf("\nLinux command line instruction:\n");
	//
    printf("\n./stoch_gen +v +F stoch_config_oneSrcOneSnk_10oct2025-18h10_rstd.par\n");
    //
    printf("\nConfiguration file (input) used with Linux command line:\n");
    printf("stoch_config_oneSrcOneSnk_10oct2025-18h10_rstd.par\n");
    //
    printf("\nDocker instruction:\n");
	//
    printf("\ndocker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image \\\n");
    printf("/stoch_gen +v +F docker_stoch_config_oneSrcOneSnk_10oct2025-18h10_rstd.par\n");
    //
    printf("\nConfiguration file (input) used with Docker instruction:\n");
    printf("docker_stoch_config_oneSrcOneSnk_10oct2025-18h10_rstd.par\n");
    //
    printf("\nBase network file (input):\n");
    printf("determ_oneSrcOneSnk_10oct2025-18h10.std\n");
    //
    printf("\nGenerated moments (output):\n");
    printf("stoch_moms_oneSrcOneSnk_10oct2025-18h10_rstd.dat\n");
    //
    printf("\nGenerated correlations (output):\n");
    printf("stoch_corrs_oneSrcOneSnk_10oct2025-18h10_rstd.dat\n");
    //
    printf("\nRaw HKW scenarios (output):\n");
    printf("stoch_rawscen_oneSrcOneSnk_10oct2025-18h10_rstd.dat\n");
    //
    printf("\nStacked scenarios of the MCFND instance in RSTD fomat (output):\n");
    printf("stoch_oneSrcOneSnk_10oct2025-18h10_rstd.txt\n");
    //
    printf("\n8.1.2.2 Origin the base deterministic network\n");
    printf("*********************************************\n");
    //
    printf("\nSimilar to 8.1.1.2.\n");
    //
    printf("\n8.1.3 Folder ./examples/oneSrcOneSnk/GSTD/\n");
    printf("******************************************\n");
    //
    printf("\nSimilar to 8.1.1, except that the base network is received in the STD format\n");
    printf("and the generated networks are saved in the GSTD format (which is identical\n");
    printf("to RSTD when there are one source and one sink for each commodity).\n");
    //
    printf("\n8.1.3.1 Operation of the stochastic generator\n");
    printf("*********************************************\n");
    //
    printf("\nLinux command line instruction:\n");
	//
    printf("\n./stoch_gen +v +F stoch_config_oneSrcOneSnk_10oct2025-18h10_gstd.par\n");
    //
    printf("\nConfiguration file (input):\n");
    printf("stoch_config_oneSrcOneSnk_10oct2025-18h10_gstd.par\n");
    //
    printf("\nDocker instruction:\n");
	//
    printf("\ndocker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image \\\n");
    printf("/stoch_gen +v +F docker_stoch_config_oneSrcOneSnk_10oct2025-18h10_gstd.par\n");
    //
    printf("\nConfiguration file (input) used with Docker instruction:\n");
    printf("docker_stoch_config_oneSrcOneSnk_10oct2025-18h10_gstd.par\n");
    //
    printf("\nConfiguration file (input) used with Linux command line:\n");
    printf("determ_oneSrcOneSnk_10oct2025-18h10.std\n");
    //
    printf("\nGenerated moments (output):\n");
    printf("stoch_moms_oneSrcOneSnk_10oct2025-18h10_gstd.dat\n");
    //
    printf("\nGenerated correlations (output):\n");
    printf("stoch_corrs_oneSrcOneSnk_10oct2025-18h10_gstd.dat\n");
    //
    printf("\nRaw HKW scenarios (output):\n");
    printf("stoch_rawscen_oneSrcOneSnk_10oct2025-18h10_gstd.dat\n");
    //
    printf("\nStacked scenarios of the MCFND instance in GSTD fomat (output):\n");
    printf("stoch_oneSrcOneSnk_10oct2025-18h10_gstd.txt\n");
    //
    printf("\n8.1.3.2 Origin the base deterministic network\n");
    printf("*********************************************\n");
    //
    printf("\nSimilar to 8.1.1.2.\n");
    //
    printf("\n8.2 Multiple but equal sources and multiple but equal sinks\n");
    printf("***********************************************************\n");
    //
    printf("\n8.2.1 Folder ./examples/eqSrcsEqSnks/GSTD/\n");
    printf("******************************************\n");
    //
    printf("\nIn the base network, there are multiple sources and multiple sinks for each\n");
    printf("commodity. However, sources and sinks are identical across commodities.\n");
    printf("There is differentiation of variable costs and capacities according to\n");
    printf("commodities. The base network must be received in the STD format (as the DOW\n");
    printf("format is unsuitable in this situation) and the generated networks must be\n");
    printf("saved in the GSTD format.\n");
    //
    printf("\n8.2.1.1 Operation of the stochastic generator\n");
    printf("*********************************************\n");
    //
    printf("\nLinux command line instruction:\n");
	//
    printf("\n./stoch_gen +v +F stoch_config_eqSrcsEqSnks_10oct2025-18h10_gstd.par\n");
    //
    printf("\nConfiguration file (input) used with Linux command line:\n");
    printf("stoch_config_eqSrcsEqSnks_10oct2025-18h10_gstd.par\n");
    //
    printf("\nDocker instruction:\n");
	//
    printf("\ndocker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image \\\n");
    printf("/stoch_gen +v +F docker_stoch_config_eqSrcsEqSnks_10oct2025-18h10_gstd.par\n");
    //
    printf("\nConfiguration file (input) used with Docker instruction:\n");
    printf("docker_stoch_config_eqSrcsEqSnks_10oct2025-18h10_gstd.par\n");
    //
    printf("\nBase network file (input):\n");
    printf("determ_eqSrcsEqSnks_10oct2025-18h10.std\n");
    //
    printf("\nGenerated moments (output):\n");
    printf("stoch_moms_eqSrcsEqSnks_10oct2025-18h10_gstd.dat\n");
    //
    printf("\nGenerated correlations (output):\n");
    printf("stoch_corrs_eqSrcsEqSnks_10oct2025-18h10_gstd.dat\n");
    //
    printf("\nRaw HKW scenarios (output):\n");
    printf("stoch_rawscen_eqSrcsEqSnks_10oct2025-18h10_gstd.dat\n");
    //
    printf("\nStacked scenarios of the MCFND instance in DOW fomat (output):\n");
    printf("stoch_eqSrcsEqSnks_10oct2025-18h10_gstd.txt\n");
    //
    printf("\n8.2.1.2 Origin of the base deterministic network\n");
    printf("************************************************\n");
    //
    printf("\nThe base deterministic network was originally produced by the determistic\n");
    printf("generator as follows.\n");
    //
    printf("\nLinux command line instruction:\n");
	//
    printf("\n./determ_gen +v +F determ_config_eqSrcsEqSnks_10oct2025-18h10.par\n");
    //
    printf("\nConfiguration file:\n");
    printf("determ_config_eqSrcsEqSnks_10oct2025-18h10.par\n");
    //
    printf("\nOutputs in DOW, STD, LP, MPS, GRA (basic graph), TEX (code of graphical\n");
    printf("representation), pdf (compiled image) formats:\n");
    printf("determ_eqSrcsEqSnks_10oct2025-18h10.dow\n");
    printf("determ_eqSrcsEqSnks_10oct2025-18h10.std\n");
    printf("determ_eqSrcsEqSnks_10oct2025-18h10.mps\n");
    printf("determ_eqSrcsEqSnks_10oct2025-18h10.lp\n");
    printf("determ_eqSrcsEqSnks_10oct2025-18h10.gra\n");
    printf("determ_eqSrcsEqSnks_10oct2025-18h10.tex\n");
    printf("determ_eqSrcsEqSnks_10oct2025-18h10.pdf\n");
    //
    exit(1);
}
