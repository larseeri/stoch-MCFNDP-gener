// Copyright (c) 2026 Eric Larsen
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


inline bool is_empty_or_whitespace(const string& str) {
    return str.empty() || all_of(str.begin(), str.end(), [](unsigned char c) {
        return isspace(c);
    });
}

inline bool isDigitsOnly(const string& str)
{   
    if (is_empty_or_whitespace(str))
		return false; 
    
    // cast to unsigned char prevents undefined behavior with negative values
    return all_of(str.begin(), str.end(), [](unsigned char c) { 
        return isdigit(c); 
    });
}

inline bool isRealNumber(const string& str)
{
    try
	{
        size_t idx;
        // Attempt to convert string to double
        stod(str, &idx);
        
        // Ensure the entire string was consumed (rejects "3.14abc")
        return (idx == str.length());
    } 
    catch (const std::invalid_argument&)
	{
        return false; // Not a number at all
    } 
    catch (const std::out_of_range&)
	{
        return false; // Number too large for a double
    }
}

// Constructor
Configuration::Configuration(int argc, char **argv, string _configFileName)
    : nbArg(argc), vectArg(argv), verbose(false), configFileName(_configFileName), instanceFormat(DOW),
      specCorrBlocks(false), momentsCharacts(UNIFORM), generMoments(false), alpha(0.25), beta(0.25),
      inputInstFileName("./base-determ-inst.dow"), skipFeasibChecks(false), outputInstFileName("./generated-stoch-inst.txt"),
      momentsFileName("./moments.dat"), correlsFileName("./correlations.dat"), numbScenarios(1000),
    //
    // distCharactMoments = 7 (7=1+2+4), ensures that the target moments
	// supplied to HKW_ScenGen() in HKW_sg.c are properly interpreted.
    // 1 -> HKW says in its documentation that this calls for "population
	// estimators (as in spreadsheets)", (In fact, this is a misnomer. 1 calls
	// for sums of terms that are divided by N-1 where N is the sample size. 
	// this is called in the statistical literature sample estimators.
	// "Population estimators" would have N as their denominator.)
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
    maxTrial(10), maxIter(50), randomSeed(7654321), randomStream(1000), outScenariosFileName("raw_scens.dat"),
    readProbs(false), probsFileName(""), readStartScenarios(false), startScenariosFileName(""), paramX(false),
    paramT(false)
{
    correlCoeffsMatrix.resize(MAX_TYPE, VectDbl(MAX_TYPE, 0.0));

    parseConfigFileNameFromCmdLine();
    processConfigFile();
    parseCmdLine();
    verify();
}

// Parses command line.
void Configuration::parseCmdLine()
{	
	if (verbose)
		cout << "****** Begin parsing configuration options from instruction line.\n" << endl;

    int elemStoch = DEMAND;

    bool reinitCorrMatrix = true;

    int a = 1;
    const char *shortKey;
	const char *longKey;
    const char *value;

    while (a < nbArg)
    {
        if ((vectArg[a][0] == PARAM_PREFIX))
        {

            shortKey = &vectArg[a][1];
			longKey = vectArg[a];

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
                cout << "[" << shortKey << "] [" << value << "]" << endl;

            switch (*shortKey)
            {
            case 'F':
				checkKeyLengthFromCmd(longKey, 2);			
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
				else
				{				
					cerr << "\nERROR:\n"
						<< "Option key F in instruction line can only be followed\n"
						<< "by one of the three following values: D, G, R.\n"
						<< "However, value [" << *value <<  "] was received.\n"
						<< "EXECUTION ABORTED" << endl;
					exit(1);                
					break;    
				}	
            case 'T':
				checkKeyLengthFromCmd(longKey, 2);
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
				else
				{				
					cerr << "\nERROR:\n"
						<< "Option key T in instruction line can only be followed\n"
						<< "by one of the two following values: U, T.\n"
						<< "However, value [" << *value <<  "] was received.\n"
						<< "EXECUTION ABORTED" << endl;
					exit(1);                
					break;    
				}
            case 'S':			
				checkKeyLengthFromCmd(longKey, 2);
				// only the last value of option S is retained 
				stochTypes.clear();
				nonStochTypes.clear();
				
				if (isDigitsOnly(value))				
					elemStoch = atoi(value);
				else
					displayPosIntError('S', value);

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
				checkKeyLengthFromCmd(longKey, 4);
                // whenever any information about correlations is supplied in
                // instruction line, any other information previously supplied
                // in configuration file is disregarded
                if (reinitCorrMatrix)
                {
                    correlCoeffsMatrix.resize(MAX_TYPE, VectDbl(MAX_TYPE, 0.0));
                    reinitCorrMatrix = false;
                }
                STOCH_ELEMS_TYPE elem[2];
                if (strlen(longKey) == 4)
                {
                    for (int j = 0; j < 2; j++)
                        switch (longKey[j + 2])
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
							cerr << "\nERROR:\n"
								<< "Option X must be specified along with a pair of characters in {D, A, C, F, V}.\n"
								<< "However, value [" << *value <<  "] was received.\n"
								<< "EXECUTION ABORTED" << endl;
							exit(1);
                            break;
                        }
                    paramX = true;
					if (isRealNumber(value))
					{
						correlCoeffsMatrix[elem[0]][elem[1]] = atof(value);
						correlCoeffsMatrix[elem[1]][elem[0]] = atof(value);
					}
					else
						displayRealError('X', value);				
					break;					
                }
                else
                {
                    cerr << "\nERROR:\n"
						<< "Option X must be specified along with a pair of characters in {D, A, C, F, V}.\n"
						<< "EXECUTION ABORTED" << endl;
					exit(1);
                    break;
                }
            case 'G':
				checkKeyLengthFromCmd(longKey, 2);
                generMoments = true;
                break;
            case 'M':
				checkKeyLengthFromCmd(longKey, 2);
                momentsFileName = value;
                break;
            case 'K':
				checkKeyLengthFromCmd(longKey, 2);			
                specCorrBlocks = true;
                break;
            case 'C':
				checkKeyLengthFromCmd(longKey, 2);
                correlsFileName = value;
                break;
            case 'A':
				checkKeyLengthFromCmd(longKey, 2);
				if (isRealNumber(value))
					alpha = atof(value);				
				else
					displayRealError('A', value);				
                break;
            case 'B':             
				checkKeyLengthFromCmd(longKey, 2);
				if (isRealNumber(value))
					beta = atof(value);				
				else
					displayRealError('B', value);				
                break;
            case 'I':
				checkKeyLengthFromCmd(longKey, 2);
                inputInstFileName = value;
                break;
			case 'J':
				checkKeyLengthFromCmd(longKey, 2);
			    skipFeasibChecks = true;
			    break;				
            case 'O':
				checkKeyLengthFromCmd(longKey, 2);
                outputInstFileName = value;
                break;

            // the following arguments are supplied to the
            // Hoyland-Kaut-Wallace algorithm
            case 'n':
				checkKeyLengthFromCmd(longKey, 2);
				if (isDigitsOnly(value))			
					numbScenarios = atoi(value);
				else
					displayPosIntError('n', value);
                break;
            case 't':
				checkKeyLengthFromCmd(longKey, 2);
				if (isDigitsOnly(value))			
					maxTrial = atoi(value);
				else
					displayPosIntError('t', value);                
                break;
            case 'i':
				checkKeyLengthFromCmd(longKey, 2);
				if (isDigitsOnly(value))			
					 maxIter = atoi(value);
				else
					displayPosIntError('i', value);
                break;
            case 'm': 
				checkKeyLengthFromCmd(longKey, 2);
				if (isRealNumber(value))
					maxErrorMoment = atof(value);				
				else
					displayRealError('m', value);				
                break;
            case 'c':          
				checkKeyLengthFromCmd(longKey, 2);
				if (isRealNumber(value))
					maxErrorCorrel = atof(value);			
				else
					displayRealError('c', value);
                break;
            case 'l':       
				checkKeyLengthFromCmd(longKey, 2);
				if (isDigitsOnly(value))			
					displayLevel = atoi(value);
				else
					displayPosIntError('l', value);				
                break;
            case 'r':  
				checkKeyLengthFromCmd(longKey, 2);
				if (isDigitsOnly(value))			
					randomSeed = atoi(value);
				else
					displayPosIntError('r', value);
                break;
            case 'R':         
				checkKeyLengthFromCmd(longKey, 2);
				if (isDigitsOnly(value))				
					randomStream = atoi(value);
				else
					displayPosIntError('R', value);				
                break;
            case 'o':
				checkKeyLengthFromCmd(longKey, 2);
                outScenariosFileName = value;
                break;
            case 'p':
				checkKeyLengthFromCmd(longKey, 2);
                readProbs = true;
                probsFileName = value;
                break;
            case 's':
				checkKeyLengthFromCmd(longKey, 2);
                readStartScenarios = true;
                startScenariosFileName = value;
                break;
			default:
				cerr << "\nERROR:\n"
					<< "Only the following option keys are admissible in\n"
					<< "the instruction line: F, T, S, X, G, M, K, C, A,\n"
					<< "B, I, J, O, n, t, i, m, c, l, r, R, o, p, s.\n"
					<< "However, the key " << *shortKey <<  " was received.\n"
					<< "EXECUTION ABORTED" << endl;
				exit (1);                
				break;				
            }
        }

        a++;
    }

	if (verbose)
		cout << "\n****** End parsing configuration options from instruction line." << endl;
}

// Parses key or (key, value) pair found on a line of the configuration file.
inline void Configuration::parseKeyValue(string key, string value)
{	
    int elemStoch = DEMAND;

    switch (key[0])
    {
    case 'F':
		checkKeyLengthFromFile(key, 2);
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
        else
		{				
			cerr << "\nERROR:\n"
				<< "Option key F in configuration file can only be followed\n"
				<< "by one of the three following values: D, G, R.\n"
				<< "However, value [" << value <<  "] was received.\n"
				<< "EXECUTION ABORTED" << endl;
			exit(1);                
			break;    
		}	
    case 'T':
		checkKeyLengthFromFile(key, 2);
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
        else
		{				
			cerr << "\nERROR:\n"
				<< "Option key T in configuration file can only be followed\n"
				<< "by one of the two following values: U, T.\n"
				<< "However, value [" << value <<  "] was received.\n"
				<< "EXECUTION ABORTED" << endl;
			exit(1);                
			break;    
		}
    case 'S':	
		checkKeyLengthFromFile(key, 2);
		// only the last value of option S is retained 
		stochTypes.clear();
		nonStochTypes.clear();
		
		if (isDigitsOnly(value))				
			elemStoch = stoi(value);
		else
			displayPosIntError('S', value);	
	
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
		checkKeyLengthFromFile(key, 4);
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
					cerr << "\nERROR:\n"
						<< "Option X must be specified along with a pair of characters in {D, A, C, F, V}.\n"
						<< "However, value [" << value <<  "] was received.\n"
						<< "EXECUTION ABORTED" << endl;
					exit(1);
				break;
                }
            paramX = true;
			
			if (isRealNumber(value))
			{
				correlCoeffsMatrix[elem[0]][elem[1]] = stof(value);
				correlCoeffsMatrix[elem[1]][elem[0]] = stof(value);
			}
			else
				displayRealError('X', value);
            break;
        }
        else
        {
			cerr << "\nERROR:\n"
				<< "Option X must be specified along with a pair of characters in {D, A, C, F, V}.\n"
				<< "EXECUTION ABORTED" << endl;
			exit(1);
			break;
		}
    case 'G':
		checkKeyLengthFromFile(key, 2);
        generMoments = true;
        break;
    case 'M':
		checkKeyLengthFromFile(key, 2);
        momentsFileName = value;
        break;
    case 'K':
		checkKeyLengthFromFile(key, 2);
        specCorrBlocks = true;
        break;
    case 'C':
		checkKeyLengthFromFile(key, 2);
        correlsFileName = value;
        break;
    case 'A':
		checkKeyLengthFromFile(key, 2);
		if (isRealNumber(value))
			alpha = stof(value);				
		else
			displayRealError('A', value);
        break;
    case 'B':
		checkKeyLengthFromFile(key, 2);
		if (isRealNumber(value))
			beta = stof(value);				
		else
			displayRealError('B', value);	
        break;
    case 'I':
		checkKeyLengthFromFile(key, 2);
        inputInstFileName = value;
        break;
	case 'J':
		checkKeyLengthFromFile(key, 2);
		skipFeasibChecks = true;
		break;
    case 'O':
		checkKeyLengthFromFile(key, 2);
        outputInstFileName = value;
        break;

    // the following arguments are supplied to the
    // Hoyland-Kaut-Wallace algorithm
    case 'n':
		checkKeyLengthFromFile(key, 2);		
		if (isDigitsOnly(value))			
			numbScenarios = stoi(value);
		else
			displayPosIntError('n', value);		
        break;
    case 't':
		checkKeyLengthFromFile(key, 2);
		if (isDigitsOnly(value))			
			maxTrial = stoi(value);
		else
			displayPosIntError('t', value); 
        break;
    case 'i':
		checkKeyLengthFromFile(key, 2);
		if (isDigitsOnly(value))			
			maxIter = stoi(value);
		else
			displayPosIntError('i', value);		
        break;
    case 'm':
		checkKeyLengthFromFile(key, 2);
		if (isRealNumber(value))
			maxErrorMoment = stof(value);				
		else
			displayRealError('m', value);		
        break;
    case 'c':
		checkKeyLengthFromFile(key, 2);
		if (isRealNumber(value))
			maxErrorCorrel = stof(value);			
		else
			displayRealError('c', value);		
        break;
    case 'l':
		checkKeyLengthFromFile(key, 2);
		if (isDigitsOnly(value))			
			displayLevel = stoi(value);
		else
			displayPosIntError('l', value);
        break;
    case 'r':
		checkKeyLengthFromFile(key, 2);
		if (isDigitsOnly(value))			
			randomSeed = stoi(value);
		else
			displayPosIntError('r', value);		
        break;
    case 'R':
		checkKeyLengthFromFile(key, 2);		
		if (isDigitsOnly(value))				
			randomStream = stoi(value);
		else
			displayPosIntError('R', value);		
        break;
    case 'o':
		checkKeyLengthFromFile(key, 2);
        outScenariosFileName = value;
        break;
    case 'p':
		checkKeyLengthFromFile(key, 2);
        readProbs = true;
        probsFileName = value;
        break;
    case 's':
		checkKeyLengthFromFile(key, 2);
        readStartScenarios = true;
        startScenariosFileName = value;
        break;
    default:
		cerr << "\nERROR:\n"
			<< "Only the following option keys are admissible in\n"
			<< "the configuration file: F, T, S, X, G, M, K, C, A,\n"
			<< "B, I, J, O, n, t, i, m, c, l, r, R, o, p, s.\n"
			<< "However, the key " << key <<  " was received.\n"
			<< "EXECUTION ABORTED" << endl;
		exit (1);                
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
        cerr << "\nERROR:\n"
			<< "No deterministic base instance file was specified.\n"
			<< "Use option -I.\n"
			<< "EXECUTION ABORTED" << endl;
		exit(1);
    }
    else if (!filesystem::exists(inputInstFileName))
    {
        cerr << "\nERROR:\n"
			<< "inputInstFileName: " << inputInstFileName << ".\n"
			<< "The specified deterministic base instance file could not be found.\n"
			<< "Use option -I.\n"
			<< "EXECUTION ABORTED" << endl;
		exit(1);
    }
    if (stochTypes.size() == 0)
    {
        cerr << "\nERROR:\n"
			<< "No stochastic type was specified.\n"
			<< "Use option -S.\n"
			<< "EXECUTION ABORTED" << endl;
		exit(1);
    }

    if (alpha < 0.0 || alpha > 0.95)
    {
        cerr << "\nERROR:\n"
			<< "alpha must fall in [0, 0.95].\n"
			<< "EXECUTION ABORTED" << endl;
		exit(1);
    }

    if (beta < 0.0)
    {
        cerr << "\nERROR:\n"
			<< "beta must be non-negative.\n"
			<< "EXECUTION ABORTED" << endl;
		exit(1);
    }

    if (!specCorrBlocks)
    {		
		if (paramX)
		{		
			cerr << "\nERROR:\n"
				<< "Option(s) X?? cannot be used without option K.\n"
				<< "EXECUTION ABORTED" << endl;
			exit(1);
		}
	}	
	else
		verifyOutputPath(correlsFileName);

    if (!generMoments)
    {		
		if (paramT)
		{		
			cerr << "\nERROR:\n"
				<< "Option T cannot be used without option G.\n"
				<< "EXECUTION ABORTED" << endl;
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
                cerr << "\nERROR:\n"
					<< "Correlations must fall in [-1, 1].\n"
					<< "EXECUTION ABORTED" << endl;
                exit(1);
            }
        }
    }	
	
	if (((find(stochTypes.begin(), stochTypes.end(), CAP_COMMOD) != stochTypes.end()) ||	
		(std::find(stochTypes.begin(), stochTypes.end(), VAR_COST) != stochTypes.end())) &&
		(instanceFormat == DOW))
	{
		cerr << "\nERROR:\n"
			<< "Randomization of variable costs and/or commodity-specific\n"
			<< "capacities is not allowed when the basic instance is received in\n"
			<< "the DOW format.\n"
			<< "EXECUTION ABORTED" << endl;	
		exit(1);
	}
	
	if (readStartScenarios)
	{
		if (!filesystem::exists(startScenariosFileName))
		{
			cerr << "\nERROR:\n"
				<< "The specified starting scenarios file " << startScenariosFileName << "\n"
				<< "could not be found.\n"
				<< "EXECUTION ABORTED" << endl;
			exit(1);
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
			cerr << "\nERROR:\n"    
				<< "Output file path " << fullPath << " was specified.\n" 
				<< "However, folder " << folderPath << " is unreachable.\n"
				<< "Make sure to create individual folders on output file path\n"
				<< "and to grant sufficient access.\n"
				<< "EXECUTION ABORTED\n" << endl;
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
			else if (strcmp(vectArg[a], CONFIG_HELP_MODE) == 0)
            {
                printHelpAndExit();
            }
            else if (isalpha(vectArg[a][1]))
 				cout << "\nWARNING:\n"
					<< "configuration prefix [" << vectArg[a] << "] unknown." << endl;
        }
    }
    if (verbose)
    {
        cout << "\n********* VERBOSE MODE *********\n" << endl;  

        if (foundFile)
            cout << "****** Configuration file name: " << configFileName << endl;
    }
}

// Processes keys or (key, value) pairs specified on a non-comment lines of
// configuration file.
bool Configuration::processConfigFile()
{
    string line;
    ifstream inStr;

    // checking if configuration file exists
    if (!filesystem::is_regular_file(configFileName))
    {
        cerr << "\nERROR:\n"
			<< "The specified configuration file " << configFileName << " could not be found.\n"
			<< "EXECUTION ABORTED" << endl;
        exit(1);
    }
	
	// attempting to open configuration file
    inStr.open(configFileName);

    if (!inStr.is_open())
    {
        cerr << "\nERROR:\n"
			<< "The specified configuration file " << configFileName << " could not be opened.\n"
			<< "EXECUTION ABORTED" << endl;
        exit(1);
    }
	
    bool found = 0;

	if (verbose)
		cout << "****** Begin parsing configuration file." <<  endl;	

    getline(inStr, line);
    while (inStr)
    {
        istringstream configLine(line);

        string key;
        configLine >> key;
        string rawValue = removeQuotes(configLine);

        int pos = (key[0] == PARAM_PREFIX) ? 1 : 0;

        // ensures that it is not a comment
        if (isalpha(key[pos]))
        {
            found = 1;
			
			// Find the first character of rawValue that is equal to '#'
			auto it = find_if(rawValue.begin(), rawValue.end(), [](unsigned char c) {
			return c=='#';
			});
    
			// Extract everything up to that character
			string cleanValue(rawValue.data(), distance(rawValue.begin(), it));
					            
            if (verbose)
                cout << "[" << key.substr(pos) << "] = [" << cleanValue << "]" << endl;
			
            parseKeyValue(key.substr(pos), cleanValue);
        }
        getline(inStr, line);
    }
	
 	if (verbose)
		cout << "****** End parsing configuration file." <<  endl;	
	
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

void Configuration::displayPosIntError(char key, string value)
{
    cerr << "\nERROR:\n"
		<< "The value [" << value << "] cannot be interpreted as a positive integer\n"
		<< "admissible for key [" << key << "]\n"
		<< "EXECUTION ABORTED\n" << endl;	
	exit(1);
}

void Configuration::displayRealError(char key, string value)
{
    cerr << "\nERROR:\n"
		<< "The value [" << value << "] cannot be interpreted as a real number\n"
		<< "admissible for key [" << key << "]\n"
		<< "EXECUTION ABORTED\n" << endl;	
	exit(1);
}

void Configuration::checkKeyLengthFromCmd(const char* key, int expLen)
{	
	if (expLen != static_cast<int> (strlen(key)))
	{
		cerr << "\nERROR:\n"
			<< "The key beginning with characters [" << string_view(key).substr(0, expLen)  << "] has expected total length equal to " << expLen << ".\n"
			<< "However, its actual length is " << strlen(key) << ".\n"
			<< "EXECUTION ABORTED\n" << endl;	
		exit(1);		
	}
}

void Configuration::checkKeyLengthFromFile(string key, int expSize)
{
	int length = expSize;
	
	if (key[0] != '-')
		length-=1;
	
	if (length != static_cast<int> (key.size()))
	{
		cerr << "\nERROR:\n"
			<< "The key beginning with characters [" << key[0] << "] has expected length equal to " << length << ".\n"
			<< "However, its actual length is " << key.size() << ".\n"
			<< "EXECUTION ABORTED\n" << endl;	
		exit(1);		
	}
}

// Unused.
void Configuration::displayError(string error)
{
    cerr << "\nERROR: " + error << endl;	
	cerr << "EXECUTION ABORTED\n" << endl;
	
	cerr << "For details about the sources of this error and its remedies\n"
			<< "either launch the program with the \"help\" option to display the user guide\n"
			<< "or examine the document USER_GUIDE.md in the online repository." << endl;		
	exit(1);
}

void Configuration::printHelpAndExit()
{    
	printf("\nHELP REQUESTED: DISPLAYING USER GUIDE\n");
	printf("\n");
	printf("GENERATOR OF INSTANCES OF THE STOCHASTIC LINEAR, TWO-STAGE, MULTI-COMMODITY,\n");
	printf("CAPACITATED, FIXED CHARGE NETWORK DESIGN PROBLEM WITH RECOURSE\n");
	printf("\n");
	printf("This application generates finite sets of scenarios matching target moments\n");
	printf("and correlations for specified subsets of parameters of a given instance of\n");
	printf("the deterministic multi-commodity, capacitated, fixed charge network design\n");
	printf("problem (MCFNDP). Details of this process are governed by user-specified\n");
	printf("configuration settings. The code implementing the algorithm presented in\n");
	printf("Hoyland, Kaut and Wallace (2003) (HKW) and built by Michal Kaut and Diego\n");
	printf("Mathieu is used to compute the scenarios.\n");
	printf("\n");
	printf("For detailed background, see\n");
	printf("Larsen, E., Bisaillon, S., Cordeau, J.F. and Frejinger, E.\n");
	printf("Pseudo-random Instance Generators in C++ for Deterministic and Stochastic\n");
	printf("Multi-commodity Network Design Problems.\n");
	printf("\n");
	printf("The two following paragraphs are partly adapted from\n");
	printf("../HOYLAND_KAUT_WALLACE/sg_HKW.c:\n");
	printf("\n");
	printf("Scenario generation code is based on 'A Heuristic for Moment-matching\n");
	printf("Scenario Generation' by K. Høyland, M. Kaut & S.W. Wallace; Computational\n");
	printf("Optimization and Applications, 24 (2-3), pp. 169–185, 2003;\n");
	printf("doi:10.1023/A:1021853807313. Code by Michal Kaut (michal.kaut@iot.ntnu.no)\n");
	printf("& Diego Mathieu.\n");
	printf("\n");
	printf("The HKW code generates scenarios for multivariate random variables by\n");
	printf("matching specified first four moments and correlations. The algorithm\n");
	printf("is iterative and gradually closes the gap from target values. Although\n");
	printf("convergence has not been demonstrated theoretically, extensive empirical\n");
	printf("testing and application have shown that convergence occurs provided the\n");
	printf("requested number of scenarios is sufficiently large. Several trials may be\n");
	printf("required, each one starting from new random initializations automatically.\n");
	printf("(See Section 6 below.)\n");
	printf("\n");
	printf("This user guide comprises the following sections:\n");
	printf("1 Building and running the program directly from the Linux command line\n");
	printf("2 Running the program through Docker instructions\n");
	printf("3 Randomized elements, target moments and target correlations\n");
	printf("4 Output file\n");
	printf("5 Formats of data files\n");
	printf("6 Convergence, target moments, target correlations and number of scenarios\n");
	printf("7 Usage of options\n");
	printf("8 Detailed examples\n");
	printf("\n");
	printf("1 Building and running the program directly from the Linux command line\n");
	printf("***********************************************************************\n");
	printf("\n");
	printf("1.1 Building the program from the Linux command line\n");
	printf("****************************************************\n");
	printf("\n");
	printf("Under Linux, follow these steps to build the program:\n");
	printf("\n");
	printf("a- Install the CPLEX solver and/or the HIGHS solver (if not already done)\n");
	printf("and note the installation folder(s). According to our verifications, CPLEX\n");
	printf("and HIGHS work equally well in our application. HIGHS is freely available\n");
	printf("under the MIT license at https://highs.dev/. See Huangfu, Q., Hall, J.A.J.\n");
	printf("Parallelizing the dual revised simplex method. Math. Prog. Comp. 10, 119-142\n");
	printf("(2018).\n");
	printf("\n");
	printf("b- Either (i) use Git to clone locally the remote repository\n");
	printf("https://github.com/larseeri/stoch-MCFNDP-gener, or (ii) generate, download\n");
	printf("locally to an appropriate folder and uncompress a .zip archive containing\n");
	printf("the remote repository file tree.\n");
	printf("\n");
	printf("c- Edit the file named \"makefile\" to set the variables indicating the\n");
	printf("locations relevant to the CPLEX solver and/or the HIGHS solver, as noted in\n");
	printf("step a-. These locations are chiefly the absolute paths of CPLEX_HOME and\n");
	printf("HIGHS_HOME. Save the modifications.\n");
	printf("\n");
	printf("d- In a Linux terminal, go to the /program folder just created.\n");
	printf("\n");
	printf("e- Make sure that the gcc C++ compiler and the make C++ build tool are\n");
	printf("available.\n");
	printf("\n");
	printf("Remark:\n");
	printf("If the gcc compiler and/or the make build tool are unavailable, consider\n");
	printf("running the program through Docker instructions (see Section 2 below).\n");
	printf("\n");
	printf("f- In the terminal, run ONLY ONE of these commands to build the program\n");
	printf("equipped either with the CPLEX solver or with the HIGHS solver:\n");
	printf("\n");
	printf("make CPLEX_SOLVER=1\n");
	printf("\n");
	printf("make HIGHS_SOLVER=1\n");
	printf("\n");
	printf("An executable file named \"stoch_gen\" will be created inside the /program\n");
	printf("folder. Building with versions 11.5.0, 12.2.0, 13.2.0, 14.2.0 and 15.1.0 of\n");
	printf("the gcc compiler has been verified.\n");
	printf("\n");
	printf("Remarks:\n");
	printf("- To build once again with HIGHS if the program has already been built with\n");
	printf("CPLEX, run the following:\n");
	printf("\n");
	printf("make clean CPLEX_SOLVER=1\n");
	printf("make HIGHS_SOLVER=1\n");
	printf("\n");
	printf("- To build once again with CPLEX if the program has already been built with\n");
	printf("HIGHS, run the following:\n");
	printf("\n");
	printf("make clean HIGHS_SOLVER=1\n");
	printf("make CPLEX_SOLVER=1\n");
	printf("\n");
	printf("1.2 Running the program with Linux command line instructions\n");
	printf("************************************************************\n");
	printf("\n");
	printf("The program can be run as follows from the /program folder just created\n");
	printf("along with default option parameter values:\n");
	printf("\n");
	printf("./stoch_gen\n");
	printf("\n");
	printf("Options for running the program can be specified directly in the command\n");
	printf("line or in a configuration file. By default, the program will attempt to\n");
	printf("read option values from file stoch_config.par. (This default setting can be\n");
	printf("updated by modifying line 1238 of file main.cpp.) The default file can\n");
	printf("be overridden with option +F:\n");
	printf("\n");
	printf("./stoch_gen +F newParamfile.par\n");
	printf("\n");
	printf("This assumes that stoch_config.par and newParamfile.par are located at the\n");
	printf("root of the /program folder. Otherwise, their absolute paths or their paths\n");
	printf("relative to the /program folder must be explicited.\n");
	printf("\n");
	printf("The default configuration file stoch_config.par has initially been filled\n");
	printf("with illustrative values and should be edited appropriately if no\n");
	printf("replacement configuration file is indicated on the command line with option\n");
	printf("key +F.\n");
	printf("\n");
	printf("Verbosity option +v should be used to obtain more information about all\n");
	printf("options that have been set in the parameter file or in the command line:\n");
	printf("\n");
	printf("./stoch_gen +v +F newParamfile.par\n");
	printf("\n");
	printf("The configuration file should have this format:\n");
	printf("\n");
	printf("optionKey1   value1\n");
	printf("optionKey2   value2\n");
	printf(".\n");
	printf(".\n");
	printf("optionKeyN   valueN\n");
	printf("\n");
	printf("It may specify at most one option per line. Option keys may be prefixed or\n");
	printf("not with '-'. Comments may appear either at the beginning of a line, in\n");
	printf("which case they must start with a non-alphabetic character (but not with\n");
	printf("'-'), or at the end of a line.\n");
	printf("\n");
	printf("The general command line for running the program is as follows:\n");
	printf("\n");
	printf("./stoch_gen [+v] [+F nameOfConfigFile] [main and HKW options]\n");
	printf("\n");
	printf("Once the configuration file has been processed, options appearing in the\n");
	printf("command line (other than +F and +v) are processed. Option values specified\n");
	printf("in the command line supersede those specified in the configuration file.\n");
	printf("That is, if included, they replace, modify and supplement the latter.\n");
	printf("Option keys appearing in the command line must be prefixed with '-' as\n");
	printf("follow:\n");
	printf("\n");
	printf("-optionkeyA valueA -optionKeyB valueB...\n");
	printf("\n");
	printf("By default, boolean values are set to false and presence of option key is\n");
	printf("sufficient to set the corresponding value to true. For options expecting\n");
	printf("other types, values should be specified.\n");
	printf("\n");
	printf("To display the complete documentation about the application and the\n");
	printf("available options, run\n");
	printf("\n");
	printf("./stoch_gen -H\n");
	printf("or\n");
	printf("./stoch_gen -h\n");
	printf("\n");
	printf("Summing up, option values are successively updated from these locations in\n");
	printf("this order:\n");
	printf("(i) hard-coded initialization values in lines 27-51 of configuration.cpp\n");
	printf("(ii) values appearing in configuration file whose name is specified in\n");
	printf("command line with +F (if it is specified and if it exists); if no valid\n");
	printf("configuration file is specified with +F in command line, values will be\n");
	printf("taken from default configuration file whose name is specified in line 1238\n");
	printf("of main.cpp; at this point, the latter should exist\n");
	printf("(iii) values specified in the command line\n");
	printf("Notice that location (i) is insufficiently transparent to the user and\n");
	printf("should not be relied upon. Parameter values should be set through (ii) and \n");
	printf("(iii).\n");
	printf("\n");
	printf("Remark:\n");
	printf("See Section 8 below for examples of configuration files and instructions\n");
	printf("when the program is run with command line instructions.\n");
	printf("\n");
	printf("2 Running the program through Docker instructions\n");
	printf("*************************************************\n");
	printf("\n");
	printf("An advantage of running the executable stoch_gen inside a Docker container\n");
	printf("is the possibility to avoid selection and installation of suitable C++\n");
	printf("compiler, C++ make build tool and required C++ libraries (e.g., on Apple\n");
	printf("computers). It suffices that the Docker engine be installed on the host\n");
	printf("where the computation will be run. For this, see\n");
	printf("https://docs.docker.com/engine/install/.\n");
	printf("\n");
	printf("Once the Docker engine is installed, follow these steps:\n");
	printf("\n");
	printf("a- Open a terminal and retrieve from the GitHub container registry the\n");
	printf("Docker image ghcr.io/larseeri/stoch_gen_image:latest as follows:\n");
	printf("\n");
	printf("docker pull ghcr.io/larseeri/stoch_gen_image:latest\n");
	printf("\n");
	printf("(From this image, Docker containers running the executable 'stoch_gen' will\n");
	printf("be generated.)\n");
	printf("\n");
	printf("b- For simplicity, define the alias 'stoch_gen_image' through\n");
	printf("\n");
	printf("docker tag ghcr.io/larseeri/stoch_gen_image:latest stoch_gen_image\n");
	printf("\n");
	printf("Remark:\n");
	printf("If the Docker image stoch_gen_image were instead available as the tar ball\n");
	printf("stoch_gen_image.tar, then it could be loaded through\n");
	printf("\n");
	printf("docker load --input stoch_gen_image.tar\n");
	printf("\n");
	printf("c- Either (i) use Git to clone locally the remote repository\n");
	printf("https://github.com/larseeri/stoch-MCFNDP-gener, or (ii) generate, download\n");
	printf("locally to an appropriate folder and uncompress a .zip archive containing\n");
	printf("the remote repository file tree. The /docker folder contains (i) the\n");
	printf("primitive files named 'Dockerfile' and 'makefile' based on which the\n");
	printf("Docker image was originally generated, (ii) a subfolder /inout used\n");
	printf("for exchanging files between the program running inside a Docker container\n");
	printf("and the directory structure of the host, (iii) some auxiliary files, and\n");
	printf("(iv) the subfolder /examples containing a number of configuration files and\n");
	printf("their resulting outputs (examined in Section 5 below).\n");
	printf("\n");
	printf("d- In a terminal, go to the /docker folder just created.\n");
	printf("\n");
	printf("Then, instructions similar to those appearing in Section 1.2 can be run from\n");
	printf("the terminal on AMD and Intel CPUs (other CPUs are addressed further below).\n");
	printf("For example, the default Docker instruction\n");
	printf("\n");
	printf("docker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image /stoch_gen\n");
	printf("\n");
	printf("is equivalent to the default Linux command line instruction\n");
	printf("\n");
	printf("./stoch_gen\n");
	printf("\n");
	printf("In general, Docker instructions can be formulated from those appearing in\n");
	printf("Section 1.2 by (i) adding the prefix \n");
	printf("\n");
	printf("docker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image\n");
	printf("\n");
	printf("and by (ii) replacing the executable location ./stoch_gen with /stoch_gen\n");
	printf("(i.e. removing the initial dot), since the executable stoch_gen is located\n");
	printf("at the root of the Docker container being generated and since locations\n");
	printf("inside a Docker container must be accessed through absolute paths.\n");
	printf("\n");
	printf("Hence, the following command line instructions:\n");
	printf("\n");
	printf("./stoch_gen +v +F newParamfile.par\n");
	printf("\n");
	printf("./stoch_gen -h\n");
	printf("\n");
	printf("./stoch_gen +v +F newParamFile.par -r 1234 -R 4567\n");
	printf("\n");
	printf("respectively correspond to the following Docker instructions:\n");
	printf("\n");
	printf("docker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image \\n");
	printf("/stoch_gen +v +F /inout/newParamfile.par\n");
	printf("\n");
	printf("docker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image \\n");
	printf("/stoch_gen -h\n");
	printf("\n");
	printf("docker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image \\n");
	printf("/stoch_gen +v +F /inout/newParamFile.par -r 1234 -R 4567\n");
	printf("\n");
	printf("Important:\n");
	printf("When running the program through Docker instructions, all input and output\n");
	printf("files must be located inside the docker/inout subfolder. Within the Docker\n");
	printf("instructions and the configuration files, this subfolder must be denoted\n");
	printf("/inout since it is located at the root of the Docker container. Hence, a\n");
	printf("replacement configuration file located at docker/inout/newParamFile.par\n");
	printf("will be specified in the Docker instruction with +F /inout/newParamFile.par.\n");
	printf("Similarly, a base network file at docker/inout/baseNetworkFile.std will be\n");
	printf("specified in the Docker instruction or inside the configuration file with\n");
	printf("-I /inout/baseNetworkFile.std.\n");
	printf("\n");
	printf("The default configuration file is located at docker/inout/stoch_config.par.\n");
	printf("It has initially been filled with illustrative values and should be edited\n");
	printf("appropriately if no replacement configuration file is indicated on the\n");
	printf("Docker instruction line with option key +F.\n");
	printf("\n");
	printf("Option parameter values are successively updated from these locations in\n");
	printf("this order:\n");
	printf("(i) initialization values hard-coded in lines 27-51 of configuration.cpp,\n");
	printf("(ii) either replacement configuration file specified in the Docker\n");
	printf("instruction line with option key +F, or default configuration file\n");
	printf("./inout/stoch_config.par if no replacement configuration file is\n");
	printf("specified in the Docker instruction line with option key +F,\n");
	printf("(iii) values specified in the Docker instruction\n");
	printf("Notice that location (i) is insufficiently transparent to the user and\n");
	printf("should not be relied upon. Parameter values should be set through (ii) and\n");
	printf("(iii).\n");
	printf("\n");
	printf("Important:\n");
	printf("The supplied Docker image conforms to the amd64 CPU instructions standard\n");
	printf("and will run natively on AMD and Intel CPUs. On CPUs conforming instead to\n");
	printf("the arm64 standard, e.g., on recent Apple computers, the supplied Docker\n");
	printf("image can be run through emulation. Essentially, this may be achieved as\n");
	printf("follows:\n");
	printf("\n");
	printf("a- Running the following preliminary instruction once. (For details, see\n");
	printf("https://hub.docker.com/r/tonistiigi/binfmt.)\n");
	printf("\n");
	printf("docker run --privileged --rm tonistiigi/binfmt --install all\n");
	printf("\n");
	printf("b- Adding the option --platform linux/amd64 to the Docker instruction line\n");
	printf("prefix as follows:\n");
	printf("\n");
	printf("docker run --platform linux/amd64 --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image\n");
	printf("\n");
	printf("The following characteristics are specific to the operation of the program\n");
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
	printf("\n");
	printf("The following Docker instructions may be useful:\n");
	printf("\n");
	printf("docker image ls (list all currently existing docker images)\n");
	printf("\n");
	printf("docker image rm --force imageName (forcefully remove imageName)\n");
	printf("\n");
	printf("docker container ls (list all running containers)\n");
	printf("\n");
	printf("docker container ls --all (list all containers, running or not)\n");
	printf("\n");
	printf("See Section 8 below for examples of configuration files and instructions\n");
	printf("when the program is run through Docker.\n");
	printf("\n");
	printf("3 Randomized elements, target moments and target correlations\n");
	printf("*************************************************************\n");
	printf("\n");
	printf("The generator rests on three primary components: a subset of randomized\n");
	printf("parameters, target moments and target correlations. We proceed to describe\n");
	printf("how each one is specified.\n");
	printf("\n");
	printf("Remarks:\n");
	printf("- Section 7 below describes detailed usage of every option. Section 8\n");
	printf("presents systematically some examples of configuration options specified in\n");
	printf("configuration files.\n");
	printf("- For the purposes of the current discussion, <...> indicates that some\n");
	printf("value or values of the type specified inside the brackets is expected. For\n");
	printf("instance, <string> means that a string is expected whereas <> means that no\n");
	printf("value is expected.\n");
	printf("\n");
	printf("3.1 Base deterministic MCFNDP instance and randomized elements\n");
	printf("**************************************************************\n");
	printf("\n");
	printf("Option S <int> must always be specified, either in the instruction line or\n");
	printf("in the configuration file. It indicates which subset of parameters of a\n");
	printf("given base deterministic MCFNDP instance will be assigned a joint probability\n");
	printf("distribution and will thus vary between scenarios. We call these parameters\n");
	printf("randomized elements. Option S expects an integer from 1 to 31 that is the\n");
	printf("sum of the following, where symbol -> is a shorthand notation for \"calls\n");
	printf("for randomization of\":\n");
	printf("1 -> demand volumes\n");
	printf("2 -> total capacities of arcs\n");
	printf("4 -> commodity-specific capacities of arcs\n");
	printf("8 -> fixed costs\n");
	printf("16 -> variable costs (either commodity-specific or not)\n");
	printf("Hence, specifying -S 3 in the instruction line indicates that total\n");
	printf("capacities of arcs (2) and demand volumes (1) must be randomized\n");
	printf("(2 + 1 = 3).\n");
	printf("\n");
	printf("By default, the name of the file describing the base deterministic MCFNDP\n");
	printf("instance is './base-determ-network.dow' (Linux instructions),\n");
	printf("'./inout/base-determ-network.dow' (Docker instructions) and the default input\n");
	printf("format is DOW. This name may be changed with option I <string> and the input\n");
	printf("format can be indicated accordingly with option F <char>. For example,\n");
	printf("specifying -S 3 -I instB.std -F G in the instruction line indicates that\n");
	printf("demands and total capacities of arcs will be randomized based on the\n");
	printf("deterministic instance stored in file instB.std using the generic STD format.\n");
	printf("(See Section 5.1 below for details.)\n");
	printf("\n");
	printf("Important:\n");
	printf("It is critical that the content of the file describing the base\n");
	printf("deterministic MCFNDP instance and the indicated format of its content match\n");
	printf("together.\n");
	printf("\n");
	printf("Remark:\n");
	printf("While the generator is sufficiently general to randomize fixed costs, the\n");
	printf("latter are usually viewed as being non-stochastic in the context of two-\n");
	printf("stage stochastic programming. Hence, in standard applications, randomization\n");
	printf("of the MCFNDP will be limited to a subset of the following: demand volumes,\n");
	printf("total capacities of arcs, commodity-specific capacities of arcs and variable\n");
	printf("costs (either commodity-specific or not).\n");
	printf("\n");
	printf("3.2 Target moments\n");
	printf("******************\n");
	printf("\n");
	printf("Target moments can be set and saved as follows:\n");
	printf("\n");
	printf("a- If option key G is included, target moments will be generated according\n");
	printf("to a distributional assumption and will be written individually in a file\n");
	printf("whose default name is './moments.dat' (Linux instructions) or\n");
	printf("'./inout/moments.dat' (Docker instructions). This name can be specified\n");
	printf("otherwise with option key M. The distributional assumption can be specified\n");
	printf("with option key T <char> where <char> stands for either U (uniform) or D\n");
	printf("(triangular) and with two parameters, alpha and beta, that can be specified\n");
	printf("using options A <double> and B <double>.\n");
	printf("\n");
	printf("b- If option key G is not included, individual target moments will be read\n");
	printf("from the file whose default name is is './moments.dat' (Linux\n");
	printf("instructions) or './inout/moments.dat' (Docker instructions). The latter\n");
	printf("can be specified otherwise with option key M. Format of the target moments\n");
	printf("file is specified in Section 5.2.1 below. Without option G, option T will be\n");
	printf("rejected and options A and B will be ignored.\n");
	printf("\n");
	printf("Important:\n");
	printf("In case b-, that is when option G is omitted and target moments are read\n");
	printf("from a file, it is crucial that they are compatible with a proper\n");
	printf("probability distribution. If this condition fails, an attempt at\n");
	printf("randomization is meaningless and convergence of the HKW algorithm is\n");
	printf("expected to fail. In case a-, that is when option G is included, the target\n");
	printf("moments that are generated are by construction compatible with a proper\n");
	printf("probability distribution.\n");
	printf("\n");
	printf("3.3 Target correlations\n");
	printf("***********************\n");
	printf("\n");
	printf("Target correlations can be set and saved as follows.\n");
	printf("\n");
	printf("a- If option key K is included, target correlations can be directly\n");
	printf("specified in blocks either (i) in the instruction line, or (ii) in the\n");
	printf("configuration file. A common target correlation coefficient within subset of\n");
	printf("parameters U or between subsets of parameters U and V is specified with\n");
	printf("X<UU> <double> and X<UV> <double> respectively, where <double> is the target\n");
	printf("correlation coefficient and U and V, can take these values, where symbol ->\n");
	printf("is a shorthand notation for \"stands for\":\n");
	printf("D -> demands\n");
	printf("A -> total capacities of arcs\n");
	printf("C -> commodity-specific capacities of arcs\n");
	printf("F -> fixed costs\n");
	printf("V -> variable costs\n");
	printf("Target correlations that are left unspecified are assumed by default to be\n");
	printf("zero. The resulting individual correlations will be written to a file whose\n");
	printf("default name is './correlations.dat' (Linux instructions) or\n");
	printf("'./inout/correlations.dat' (Docker instructions). This name can be\n");
	printf("specified otherwise with option key C.\n");
	printf("\n");
	printf("Important:\n");
	printf("Cases (i) and (ii) are not coexistent: if any correlation information is\n");
	printf("supplied in the instruction line, any correlation information supplied in\n");
	printf("the configuration file is ignored.\n");
	printf("\n");
	printf("b- If option key K is not included, individual correlations will be read\n");
	printf("from a file whose default name is './correlations.dat' (Linux instructions)\n");
	printf("or './inout/correlations.dat' (Docker instructions). The latter can be\n");
	printf("specified otherwise with option key C. Format of the target correlations\n");
	printf("file is specified in Section 5.2.2 below. Without option K, options X?? will\n");
	printf("be rejected.\n");
	printf("\n");
	printf("Important:\n");
	printf("In case a-, when option key K is included, it is crucial that the block\n");
	printf("correlations that are specified result in a positive definite correlation\n");
	printf("matrix. Similarly, in case b-, when option key K is omitted, the supplied\n");
	printf("matrix must define a positive definite correlation matrix. Whenever the\n");
	printf("alleged correlation matrix, either constructed (case a-) or read from file\n");
	printf("(case b-), fails positive definitess, its Cholesky decomposition will fail\n");
	printf("in the HKW algorithm. A suggestion to verify and correct the alleged\n");
	printf("correlation matrix will be displayed and the program will be interrupted.\n");
	printf("\n");
	printf("3.4 Example\n");
	printf("***********\n");
	printf("\n");
	printf("Consider this hypothetical Linux command line instruction:\n");
	printf("\n");
	printf("./stoch_gen -S 3 -I instB.std -F G -K -XDD 0.5 -XDA -0.3 -XAA 0.7 -G -T U -A 0.25 -B 0.25\n");
	printf("\n");
	printf("It is interpreted as follows:\n");
	printf("-S 3: demands and total capacities of arcs will be randomized\n");
	printf("-I instB.std: base deterministic network is to be read from file instB.std\n");
	printf("-F G: format of the latter is generic STD\n");
	printf("-K: target correlations are specified in blocks (based on explicit and\n");
	printf("default values)\n");
	printf("-XDD 0.5: individual target cross-correlations among demands are all equal\n");
	printf("to 0.5 and individual target self-correlations are all equal 1.0\n");
	printf("-XDA -0.3: individual target cross-correlations between demands and total\n");
	printf("capacities of arcs are all equal -0.3\n");
	printf("-XAA 0.7: individual target cross-correlations among total capacities of\n");
	printf("arcs are all equal to 0.7 and individual target self-correlations are all\n");
	printf("equal to 1.0\n");
	printf("-G: target moments are generated (based on explicit and default\n");
	printf("specifications)\n");
	printf("-T U: target moments are generated under the assumption that\n");
	printf("randomizations each obey a uniform distribution\n");
	printf("-A 0.25: parameter alpha defining the latter distribution is 0.25\n");
	printf("-B 0.25: parameter beta defining the latter distribution is 0.25\n");
	printf("\n");
	printf("See Section 7 for additional details about each option key.\n");
	printf("\n");
	printf("4 Output file\n");
	printf("*************\n");
	printf("\n");
	printf("By default, the name of the file containing the MCFNDP instances generated\n");
	printf("by the application is './generated-stoch-network.txt' (Linux instructions) or\n");
	printf("'./inout/generated-stoch-network.txt' (Docker instructions). This can be\n");
	printf("changed with option O. This text file superposes representations of complete\n");
	printf("deterministic MCFNDP instances, one for each scenario realization. Each\n");
	printf("representation is preceded by a separator as follows:\n");
	printf("\n");
	printf("SCENARIO   <int>\n");
	printf("\n");
	printf("Immediately after the separator, the generated MCFNDP instance corresponding\n");
	printf("to the current scenario is described in its entirety with a format (DOW,\n");
	printf("restricted STD or generic STD; see Section 5.1 below), matching that used to\n");
	printf("supply the base deterministic MCFNDP instance. Stacking and separating full\n");
	printf("descriptions of the MCFNDP realizations in this fashion makes it easy to\n");
	printf("pinpoint and read the first and second stage information relevant to every\n");
	printf("individual scenario.\n");
	printf("\n");
	printf("Remark:\n");printf("If specified, feasibility of the second stage of each generated MCFNDP\n");
	printf("instance is verified with the selected solver (CPLEX or HIGHS) after opening\n");
	printf("all arcs. Only the generated networks that are demonstrated feasible are\n");
	printf("added to the output file. A report indicating how many scenarios were\n");
	printf("generated and how many among the latter are feasible is produced. A warning\n");
	printf("is issued if one or more scenarios are infeasible.\n");
	printf("\n");
	printf("5 Formats of data files\n");
	printf("***********************\n");
	printf("\n");
	printf("5.1 Base deterministic MCFNDP instance and generated stochastic MCFNDP instance\n");
	printf("*******************************************************************************\n");
	printf("\n");
	printf("DOW and STD formats are commonly used to store deterministic MCFNDP\n");
	printf("instances in text files. These formats are described in files\n");
	printf("DOW-format-desc.txt and STD-format-desc.txt. The application can read from\n");
	printf("these two formats and, in addition, from a restricted version of STD. Thus,\n");
	printf("three formats are available for supplying the base deterministic MCFNDP\n");
	printf("instance: DOW, STD (calling it generic STD, or GSTD) and a restricted\n");
	printf("version of STD (calling it restricted STD, or RSTD). One of these three\n");
	printf("formats may be selected with option F <char>.\n");
	printf("a- DOW presupposes that (i) for each arc, capacities and variable\n");
	printf("costs are identical for all commodities and (ii) each commodity has a\n");
	printf("single source node (at which volume > 0) and a single sink node (at which\n");
	printf("volume < 0).\n");
	printf("b- Generic STD imposes neither (i), nor (ii).\n");
	printf("c- Restricted STD imposes only (ii).\n");
	printf("\n");
	printf("Important:\n");
	printf("When the base deterministic MCFNDP instance is received in the generic STD\n");
	printf("format, there are potentially multiple source and/or sink nodes for some or\n");
	printf("all commodities. In this case, whenever randomization of commodity volumes\n");
	printf("is requested through option S, it is handled through these steps:\n");
	printf("a- For each commodity, the total volume, the flows out of each source node\n");
	printf("and the flows into each sink node occurring in the base deterministic MCFNDP\n");
	printf("instance are saved for ulterior computations while reading (see\n");
	printf("InstanceGStd::read() in generate.cpp).\n");
	printf("b- Total commodity volumes are randomized (along with other stochastic\n");
	printf("elements if there are any) as they are when there are single source and sink\n");
	printf("nodes for each commodity.\n");
	printf("c- For each commodity, the randomized flows out of each source node and\n");
	printf("randomized flows into each sink node are calculated by distributing the\n");
	printf("randomized total commodity volume using the original ratios of inflows and\n");
	printf("outflows over total volume obtained in step a-.\n");
	printf("\n");
	printf("Remark:\n");
	printf("The MCFNDP instances generated by the application (one for each scenario) are\n");
	printf("expressed in the same format as that used for supplying the base\n");
	printf("deterministic MCFNDP instance. Their descriptions are stacked and returned in\n");
	printf("the output file whose default name './generated-stoch-network.txt' (Linux\n");
	printf("instructions), './inout/generated-stoch-network.txt' (Docker instructions) can be\n");
	printf("modified with option O. (See Section 4 above.)\n");
	printf("\n");
	printf("5.2 Input files supplied to HKW algorithm\n");
	printf("*****************************************\n");
	printf("(partly adapted from ../HOYLAND_KAUT_WALLACE/sg_HKW.c)\n");
	printf("\n");
	printf("5.2.1 Target moments\n");
	printf("********************\n");
	printf("\n");
	printf("When option G is omitted from instruction line and from configuration file,\n");
	printf("target moments are read from a file named './moments.dat' (Linux\n");
	printf("instructions) or './inout/moments.dat' (Docker instructions) unless this\n");
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
	printf("\n");
	printf("Important:\n");
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
	printf("\n");
	printf("5.2.2 Target correlations\n");
	printf("*************************\n");
	printf("\n");
	printf("When option K is omitted from instruction line and from configuration file,\n");
	printf("target correlations are read from a file named './correlations.dat' (Linux\n");
	printf("instructions) or './inout/correlations.dat' (Docker instructions), unless\n");
	printf("this is overridden with option C. The file containing the target\n");
	printf("correlations matrix is formatted as follows (see Section 8 for examples):\n");
	printf("a- 1st and 2nd lines are identical and show the total number of individual\n");
	printf("randomized elements.\n");
	printf("b- Remaining lines of the file contain the target correlations matrix, one\n");
	printf("row of the matrix per line of the file. Individual correlations are\n");
	printf("separated by spaces.\n");
	printf("c- The correlation pairs appearing in the matrix are ordered by row and by\n");
	printf("column. Order in a row and order in a column is similar to that in a row\n");
	printf("of the matrix of target moments (see Section 4.1.1).\n");
	printf("d- Any additional line in the file is ignored.\n");
	printf("\n");
	printf("5.2.3 Probabilities of scenarios\n");
	printf("********************************\n");
	printf("\n");
	printf("If option p is omitted, then scenarios are assumed to be equiprobable.\n");
	printf("Option p can be used to supply the name of a file containing probabilities\n");
	printf("of all scenarios. The following format must be satisfied (see supplied\n");
	printf("example probs_100_eq_example.dat):\n");
	printf("a- 1st line states the number of scenarios\n");
	printf("b- 2nd line states all probabilities, separated by spaces\n");
	printf("c- Remaining lines of the file are ignored\n");
	printf("This file would be supplied as './probs_100_eq_example.dat' (Linux\n");
	printf("instructions) or './inout/probs_100_eq_example.dat' (Docker instructions).\n");
	printf("\n");
	printf("5.2.4 Raw scenarios in HKW format\n");
	printf("*********************************\n");
	printf("\n");
	printf("Raw scenarios are saved in HKW format to './raw_scens.dat' (Linux\n");
	printf("instructions) or './inout/raw_scens.dat' (Docker instructions), unless this\n");
	printf("is overriden with option o. Raw scenarios saved in previous run can be used\n");
	printf("as starting values in a following run (using option s). These will only be\n");
	printf("used in the first trial of the current run of the HKW algorithm. If this\n");
	printf("trial does not converge and the maximum number of trials specified by option\n");
	printf("t is larger than 1, the next trial starts from random guesses.\n");
	printf("\n");
	printf("Important:\n");
	printf("Raw scenarios pertain to the individual stochastic elements. They are\n");
	printf("distinct from the MCFNDP instance scenarios produced by the application and\n");
	printf("described in Section 4 above.\n");
	printf("\n");
	printf("6 Convergence, target moments, target correlations and number of scenarios\n");
	printf("**************************************************************************\n");
	printf("\n");
	printf("Remark:\n");
	printf("The level of on-screen reporting by the HKW algorithm should be set at\n");
	printf("no less than the default value of 3 (option l) to adequately follow the\n");
	printf("convergence process.\n");
	printf("\n");
	printf("The HKW algorithm might fail to satisfy the specified target moments and\n");
	printf("correlations within the default tolerances with the requested number of\n");
	printf("scenarios within the default number of trials and iterations. If\n");
	printf("convergence fails, the following course of action is recommended:\n");
	printf("\n");
	printf("1- Verify that target moments and correlations are proper.\n");
	printf("\n");
	printf("If option G has not been selected, target moments are read from a file. It\n");
	printf("should be ascertained that they are compatible with a probability\n");
	printf("distribution. If this condition fails, convergence of the HKW algorithm is\n");
	printf("expected to fail. By construction, target moments that are generated when\n");
	printf("option G is selected are correctly specified.\n");
	printf("\n");
	printf("If option K has not been selected, a target correlation matrix is read from\n");
	printf("a file. it should be ascertained that this matrix is well-specified\n");
	printf("(square, symmetric, all elements in [-1,1], diagonal elements equal to 1,\n");
	printf("positive definite). If this condition fails, the HKW algorithm is expected\n");
	printf("to fail.\n");
	printf("\n");
	printf("If option K has been selected, a target correlation matrix is built from\n");
	printf("block correlations specified with options X??. it should be ascertained\n");
	printf("that the resulting alleged target correlation matrix is well-specified.\n");
	printf("\n");
	printf("Remarks:\n");
	printf("- If the alleged target correlation matrix fails positive definiteness, the\n");
	printf("HKW algorithm initially (i.e. before any iteration) reports a failure of\n");
	printf("its Cholesky decomposition and displays an explanation.\n");
	printf("- However, failure of the Cholesky decomposition during the iterations of\n");
	printf("the HKW algorithm is a distinct issue that might be due to an insufficient\n");
	printf("requested number of scenarios (see Paragraph 3 below).\n");
	printf("\n");
	printf("2- Verify that slow convergence is not responsible.\n");
	printf("\n");
	printf("If both target moments and target correlations are known to be properly\n");
	printf("specified, consider launching a new run after increasing the maximum numbers\n");
	printf("of trials and iterations, respectively options t and i. Doubling their\n");
	printf("values from default is suggested.\n");
	printf("\n");
	printf("Remark:\n");
	printf("Default values for the numbers of trials and iterations are usually\n");
	printf("sufficient to reach convergence under the default error tolerances. Failure\n");
	printf("of convergence is most frequently due to the causes examined in Paragraphs 1\n");
	printf("and 3 above and below.\n");
	printf("\n");
	printf("3- Adjust the number of scenarios.\n");
	printf("\n");
	printf("If both target moments and target correlations are properly specified and\n");
	printf("slow convergence is excluded, then the requested number of scenarios is too\n");
	printf("small to support the specified target moments and target correlations.\n");
	printf("Hence, the requested number of scenarios, option n, may be progressively\n");
	printf("increased by a moderate factor (e.g., between 1.5 and 5) until the algorithm\n");
	printf("converges. Once it has converged, the number of requested scenarios may be\n");
	printf("progressively adjusted downward, if desired, while the algorithm is still\n");
	printf("converging.\n");
	printf("\n");
	printf("Remarks:\n");
	printf("- The lower bound on the number of requested scenarios leading to\n");
	printf("convergence is positively related to the following: number of random\n");
	printf("elements, magnitude of 3rd and 4th target moments, heterogeneity of target\n");
	printf("correlations, heterogeneity of probabilities.\n");
	printf("- A rule of thumb is that the number of scenarios should be at least as\n");
	printf("large as the number of randomized elements.\n");
	printf("\n");
	printf("7 Usage of options\n");
	printf("******************\n");
	printf("\n");
	printf("(Some parts are adapted from ../HOYLAND_KAUT_WALLACE/sg_HKW.c.)\n");
	printf("\n");
	printf("Except for option key S (and leaving aside option key F when it is prefixed\n");
	printf("by '+' in the instruction line), non-boolean user-selected options that are\n");
	printf("specified neither in the instruction line nor the configuration file are set\n");
	printf("to the values that are hard-coded in the constructor in lines 27-51 of\n");
	printf("configuration.cpp. Boolean user-selected options that are specified neither\n");
	printf("in the instruction line nor in the configuration file are set to false by\n");
	printf("default.\n");
	printf("\n");
	printf("Notation:\n");
	printf("For the purposes of this presentation, <...> indicates that some value or\n");
	printf("values of the type specified inside the brackets is expected. For instance,\n");
	printf("<string> means that a string is expected whereas <> means that no value is\n");
	printf("expected.\n");
	printf("\n");
	printf("7.1 Options available in instruction line only (Linux or Docker)\n");
	printf("****************************************************************\n");
	printf("\n");
	printf("+F <string>; name of configuration file supplying option parameter values;\n");
	printf("this replaces the default file whose name is hard-coded in line 1238 of\n");
	printf("main.cpp; any option value supplied in instruction line replaces value of\n");
	printf("same option supplied in configuration file\n");
	printf("\n");
	printf("+v <> (no value supplied); requests verbosity while parsing options in\n");
	printf("instruction line or configuration file\n");
	printf("\n");
	printf("Remark:\n");
	printf("Unknown option keys prefixed with '+' in the instruction line are\n");
	printf("disregarded.\n");
	printf("\n");
	printf("7.2 Options available in instruction line (Linux or Docker) and in configuration file\n");
	printf("*************************************************************************************\n");
	printf("\n");
	printf("The following options must be prefixed by '-' when they are called from the\n");
	printf("instruction line. This prefix can be omitted when the options are called\n");
	printf("from a configuration file. To supply an option, a line of a configuration\n");
	printf("file must begin with a dash '-' or a letter. A line beginning with any other\n");
	printf("character will be considered a comment. The remainder of a line beyond the\n");
	printf("option key and its value (if any is expected) is disregarded and may\n");
	printf("therefore be used for adding comments.\n");
	printf("\n");
	printf("Remark:\n");
	printf("Unknown option keys following '-' in the instruction line, stray characters\n");
	printf("or strings appearing without the prefix '-' in the instruction line and\n");
	printf("unknown option keys (following '-' or not) at beginning of line in a\n");
	printf("configuration file are disregarded.\n");
	printf("\n");
	printf("7.2.1 Main generator options\n");
	printf("****************************\n");
	printf("\n");
	printf("H <> (no value supplied); display complete documentation about the\n");
	printf("application and the available options\n");
	printf("\n");
	printf("h <> (no value supplied); display complete documentation about the\n");
	printf("application and the available options\n");
	printf("\n");
	printf("F <char>; format of file where base deterministic MCFNDP instance is to be\n");
	printf("read; D: DOW, G: generic STD, R: restricted STD; default: D\n");
	printf("\n");
	printf("I <string>; name of file where base deterministic MCFNDP instance is to be\n");
	printf("read; default value is './base-determ-network.dow' (Linux instructions) or\n");
	printf("'./inout/base-determ-network.dow' (Docker instructions)\n");
	printf("\n");
	printf("O <string>; name of file where MCFNDP instance resulting from\n");
	printf("scenario generation is to be written; default value is\n");
	printf("'./generated-stoch-network.txt' (Linux instructions) or\n");
	printf("'./inout/generated-stoch-network.txt' (Docker instructions)\n");
	printf("\n");
	printf("J <> (no value supplied); write the MCFNDP instance resulting from\n");
	printf("scenario generation without first checking for feasibility of each\n");
	printf("scenario and reporting on results.\n");
	printf("\n");
	printf("S <int>; identifies which subsets of parameters should vary between\n");
	printf("scenarios; expects a number (from 1 to 31) which is a sum of the following:\n");
	printf("1 -> demands\n");
	printf("2 -> total capacities of arcs\n");
	printf("4 -> commodity-specific capacities of arcs\n");
	printf("8 -> fixed costs\n");
	printf("16 -> variable costs\n");
	printf("\n");
	printf("K <> (no value supplied); if specified, target correlations are directly\n");
	printf("specified in blocks using options X and written to file\n");
	printf("'./correlations.dat' (Linux instructions) or './inout/correlations.dat'\n");
	printf("(Docker instructions), unless this default name is changed with option C;\n");
	printf("if omitted, target correlations are read from file './correlations.dat'\n");
	printf("(Linux instructions) or './inout/correlations.dat' (Docker instructions),\n");
	printf("unless this default name is changed with option C\n");
	printf("\n");
	printf("C <string>; name of file where target correlations are read if K is omitted\n");
	printf("and where generated target correlations are written if K is specified;\n");
	printf("default: './correlations.dat' (Linux instructions) or\n");
	printf("'./inout/correlations.dat' (Docker instructions)\n");
	printf("\n");
	printf("X<char><char> <double>; identifies target correlations within or between\n");
	printf("subsets of parameters that are randomized (i.e, that vary between\n");
	printf("scenarios); state once for each pair of subsets; target self-correlations\n");
	printf("are equal to 1.0; any use of X without G will be rejected;\n");
	printf("default: 0.0;\n");
	printf("D -> demands\n");
	printf("A -> total capacities of arcs\n");
	printf("C -> commodity-specific capacities of arcs\n");
	printf("F -> fixed costs\n");
	printf("V -> variable costs\n");
	printf("Example: in the instruction line, -XDD 0.9 -XFD 0.1 -XFF 0.01 specifies\n");
	printf("target individual cross-correlations between distinct demand parameters\n");
	printf("all equal to 0.9 (all individual self-correlations being equal to 1.0),\n");
	printf("between individual fixed cost parameters and demand parameters all equal\n");
	printf("to 0.1, between distinct fixed cost parameters all equal to 0.01 (all\n");
	printf("individual self-correlations being equal to 1.0)\n");	
	printf("\n");
	printf("Remark:\n");
	printf("Options -XST and -XTS are equivalent.\n");
	printf("\n");
	printf("G <> (no value supplied); if specified, target moments are generated based on\n");
	printf("option T and written to file './moments.dat' (Linux instructions) or\n");
	printf("'./inout/moments.dat' (Docker instructions), unless this default name is\n");
	printf("changed with option M; if omitted, target moments are read from file\n");
	printf("'./moments.dat' (Linux instructions) or './inout/moments.dat' (Docker\n");
	printf("instructions) unless this default name is changed with option M\n");
	printf("\n");
	printf("M <string>; name of file where target moments are read if option G is\n");
	printf("omitted or where generated target moments are written if option G is\n");
	printf("specified; default: './moments.dat' (Linux instructions) or\n");
	printf("'./inout/moments.dat' (Docker instructions)\n");
	printf("\n");
	printf("T <char>; distributional characteristics of generated target moments\n");
	printf("(U : UNIFORM , T : TRIANGULAR); any use of option T without option G\n");
	printf("will be rejected; default: U\n");
	printf("\n");
	printf("A <double>; parameter alpha used in generating target moments (see\n");
	printf("discussion immediately below); without option G, any use of option A will be\n");
	printf("ignored; default: 0.25\n");
	printf("\n");
	printf("B <double>; parameter beta used in generating target moments (see\n");
	printf("discussion immediately below); without option G, any use of option B will be\n");
	printf("ignored; default: 0.25\n");
	printf("\n");
	printf("When option T is set to U, then the stochastic elements of the network that\n");
	printf("are changed according to scenarios are each assumed to follow a\n");
	printf("uniform(a, b) distribution, where\n");
	printf("a = D - (alpha * D),\n");
	printf("b = D + (beta * D),\n");
	printf("D is the base value of the stochastic element taken from the base\n");
	printf("deterministic network,\n");
	printf("alpha must be in [0, 0.99],\n");
	printf("beta >= 0.\n");
	printf("The distribution is symmetric around D when alpha = beta.\n");
	printf("\n");
	printf("When option T is set to T, then the stochastic elements of the network that\n");
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
	printf("\n");
	printf("7.2.2 Options supplied to HKW algorithm\n");
	printf("***************************************\n");
	printf("\n");
	printf("The scenario generation code is based on article 'A Heuristic for Moment-\n");
	printf("matching Scenario Generation' by K. Høyland, M. Kaut & S.W. Wallace;\n");
	printf("Computational Optimization and Applications, 24 (2-3), pp. 169–185, 2003;\n");
	printf("doi:10.1023/A:1021853807313. Code by Michal Kaut (michal.kaut@iot.ntnu.no)\n");
	printf("& Diego Mathieu. Additional details about HKW options available in\n");
	printf("../HOYLAND_KAUT_WALLACE/sg_HKW.c.\n");
	printf("\n");
	printf("n <int>; number of scenarios to generate; default: 1000\n");
	printf("\n");
	printf("t <int>; maximum number of trials (attempts to generate scenarios using\n");
	printf("alternative random starting values); default: 10\n");
	printf("\n");
	printf("i <int>; maximum number of iterations in a trial; default: 50\n");
	printf("\n");
	printf("m <double>; maximum RMS error in matching moments (scaled to var=1);\n");
	printf("default: 0.001\n");
	printf("\n");
	printf("c <double>; maximum RMS error in matching correlations;\n");
	printf("default: 0.001\n");
	printf("\n");
	printf("l <int>; level of on-screen reporting by HKW algorithm (between 0 and 11);\n");
	printf("default: 3\n");
	printf("\n");
	printf("o <string>; name of file where resulting matrix of raw scenarios in HKW\n");
	printf("format is written; default: './raw_scens.dat' (Linux instructions) or\n");
	printf("'./inout/raw_scens.dat' (Docker instructions)\n");
	printf("\n");
	printf("s <string>; if specified, name of file where a matrix of starting scenarios\n");
	printf("saved in HKW format in a previous run (possibly using option o) is read and\n");
	printf("used for starting values; no default\n");
	printf("\n");
	printf("p <string>; if specified, read probabilities from this file; otherwise,\n");
	printf("scenarios are assumed equiprobable; no default\n");
	printf("\n");
	printf("r <int>; seed of the PCG pseudo-random generator; default: 7654321\n");
	printf("\n");
	printf("R <int>; stream of the PCG pseudo-random generator; default: 1000\n");
	printf("\n");
	printf("8 Detailed examples\n");
	printf("*******************\n");
	printf("\n");
	printf("The program/examples subfolder holds examples of configuration files and \n");
	printf("instruction lines occurring when the executable of the program is run\n");
	printf("directly from the Linux instruction line. It also contains the output files\n");
	printf("resulting from each pair of configuration file and instruction line.\n");
	printf("\n");
	printf("Similarly, the docker/examples subfolder holds examples of configuration\n");
	printf("files and instruction lines occurring when the executable of the program is\n");
	printf("run through Docker instructions. These examples are identified by the\n");
	printf("addition of the prefix 'docker_'. The subfolder also contains the output\n");
	printf("files resulting from each pair of configuration file and instruction line.\n");
	printf("The former are identical to the files generated when the executable is run\n");
	printf("from the Linux instruction line.\n");
	printf("\n");
	printf("8.1 One source and one sink for each commodity\n");
	printf("**********************************************\n");
	printf("\n");
	printf("8.1.1 Folder /examples/oneSrcOneSnk/DOW/\n");
	printf("*****************************************\n");
	printf("\n");
	printf("In the base network, there is one source and one sink for each commodity and\n");
	printf("no differentiation of variable costs and capacities according to\n");
	printf("commodities. The base network is received in the DOW format and the\n");
	printf("generated networks are saved in the DOW format.\n");
	printf("\n");
	printf("8.1.1.1 Operation of the stochastic generator\n");
	printf("*********************************************\n");
	printf("\n");
	printf("Linux command line instruction:\n");
	printf("\n");
	printf("./stoch_gen +v +F stoch_config_oneSrcOneSnk_12oct2026-17h09_dow.par\n");
	printf("\n");
	printf("Configuration file (input) used with Linux command line:\n");
	printf("stoch_config_oneSrcOneSnk_12oct2026-17h09_dow.par\n");
	printf("\n");
	printf("Docker instruction:\n");
	printf("\n");
	printf("docker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image \\n");
	printf("/stoch_gen +v +F /inout/docker_stoch_config_oneSrcOneSnk_12oct2026-17h09_dow.par\n");
	printf("\n");
	printf("Configuration file (input) used with Docker instruction:\n");
	printf("docker_stoch_config_oneSrcOneSnk_12oct2026-17h09_dow.par\n");
	printf("\n");
	printf("Base network file (input):\n");
	printf("determ_oneSrcOneSnk_12oct2026-17h09.dow\n");
	printf("\n");
	printf("Generated moments (output):\n");
	printf("stoch_moms_oneSrcOneSnk_12oct2026-17h09_dow.dat\n");
	printf("\n");
	printf("Generated correlations (output):\n");
	printf("stoch_corrs_oneSrcOneSnk_12oct2026-17h09_dow.dat\n");
	printf("\n");
	printf("Raw HKW scenarios (output):\n");
	printf("stoch_rawscen_oneSrcOneSnk_12oct2026-17h09_dow.dat\n");
	printf("\n");
	printf("Stacked scenarios of the MCFNDP instance in DOW fomat (output):\n");
	printf("stoch_oneSrcOneSnk_12oct2026-17h09_dow.txt\n");
	printf("\n");
	printf("8.1.1.2 Origin of the base deterministic network\n");
	printf("************************************************\n");
	printf("\n");
	printf("The base deterministic network was originally produced by the determistic\n");
	printf("generator as follows.\n");
	printf("\n");
	printf("Linux command line instruction:\n");
	printf("\n");
	printf("./determ_gen +v +F determ_config_oneSrcOneSnk_12oct2026-17h09.par\n");
	printf("\n");
	printf("Configuration file:\n");
	printf("determ_config_oneSrcOneSnk_12oct2026-17h09.par\n");
	printf("\n");
	printf("Outputs in DOW, STD, LP, MPS, GRA (basic graph), TEX (code of graphical\n");
	printf("representation), pdf (compiled image) formats:\n");
	printf("determ_oneSrcOneSnk_12oct2026-17h09.dow\n");
	printf("determ_oneSrcOneSnk_12oct2026-17h09.std\n");
	printf("determ_oneSrcOneSnk_12oct2026-17h09.mps\n");
	printf("determ_oneSrcOneSnk_12oct2026-17h09.lp\n");
	printf("determ_oneSrcOneSnk_12oct2026-17h09.gra\n");
	printf("determ_oneSrcOneSnk_12oct2026-17h09.tex\n");
	printf("determ_oneSrcOneSnk_12oct2026-17h09.pdf\n");
	printf("\n");
	printf("8.1.2 Folder /examples/oneSrcOneSnk/RSTD/\n");
	printf("******************************************\n");
	printf("\n");
	printf("Similar to 8.1.1, except that the base network is received in the STD format\n");
	printf("and the generated networks are saved in the RSTD format.\n");
	printf("\n");
	printf("8.1.2.1 Operation of the stochastic generator\n");
	printf("*********************************************\n");
	printf("\n");
	printf("Linux command line instruction:\n");
	printf("\n");
	printf("./stoch_gen +v +F stoch_config_oneSrcOneSnk_12oct2026-17h09_rstd.par\n");
	printf("\n");
	printf("Configuration file (input) used with Linux command line:\n");
	printf("stoch_config_oneSrcOneSnk_12oct2026-17h09_rstd.par\n");
	printf("\n");
	printf("Docker instruction:\n");
	printf("\n");
	printf("docker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image \\n");
	printf("/stoch_gen +v +F /inout/docker_stoch_config_oneSrcOneSnk_12oct2026-17h09_rstd.par\n");
	printf("\n");
	printf("Configuration file (input) used with Docker instruction:\n");
	printf("docker_stoch_config_oneSrcOneSnk_12oct2026-17h09_rstd.par\n");
	printf("\n");
	printf("Base network file (input):\n");
	printf("determ_oneSrcOneSnk_12oct2026-17h09.std\n");
	printf("\n");
	printf("Generated moments (output):\n");
	printf("stoch_moms_oneSrcOneSnk_12oct2026-17h09_rstd.dat\n");
	printf("\n");
	printf("Generated correlations (output):\n");
	printf("stoch_corrs_oneSrcOneSnk_12oct2026-17h09_rstd.dat\n");
	printf("\n");
	printf("Raw HKW scenarios (output):\n");
	printf("stoch_rawscen_oneSrcOneSnk_12oct2026-17h09_rstd.dat\n");
	printf("\n");
	printf("Stacked scenarios of the MCFNDP instance in RSTD fomat (output):\n");
	printf("stoch_oneSrcOneSnk_12oct2026-17h09_rstd.txt\n");
	printf("\n");
	printf("8.1.2.2 Origin the base deterministic network\n");
	printf("*********************************************\n");
	printf("\n");
	printf("Similar to 8.1.1.2.\n");
	printf("\n");
	printf("8.1.3 Folder /examples/oneSrcOneSnk/GSTD/\n");
	printf("******************************************\n");
	printf("\n");
	printf("Similar to 8.1.1, except that the base network is received in the STD format\n");
	printf("and the generated networks are saved in the GSTD format (which is identical\n");
	printf("to RSTD when there are one source and one sink for each commodity).\n");
	printf("\n");
	printf("8.1.3.1 Operation of the stochastic generator\n");
	printf("*********************************************\n");
	printf("\n");
	printf("Linux command line instruction:\n");
	printf("\n");
	printf("./stoch_gen +v +F stoch_config_oneSrcOneSnk_12oct2026-17h09_gstd.par\n");
	printf("\n");
	printf("Configuration file (input):\n");
	printf("stoch_config_oneSrcOneSnk_12oct2026-17h09_gstd.par\n");
	printf("\n");
	printf("Docker instruction:\n");
	printf("\n");
	printf("docker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image \\n");
	printf("/stoch_gen +v +F /inout/docker_stoch_config_oneSrcOneSnk_12oct2026-17h09_gstd.par\n");
	printf("\n");
	printf("Configuration file (input) used with Docker instruction:\n");
	printf("docker_stoch_config_oneSrcOneSnk_12oct2026-17h09_gstd.par\n");
	printf("\n");
	printf("Configuration file (input) used with Linux command line:\n");
	printf("determ_oneSrcOneSnk_12oct2026-17h09.std\n");
	printf("\n");
	printf("Generated moments (output):\n");
	printf("stoch_moms_oneSrcOneSnk_12oct2026-17h09_gstd.dat\n");
	printf("\n");
	printf("Generated correlations (output):\n");
	printf("stoch_corrs_oneSrcOneSnk_12oct2026-17h09_gstd.dat\n");
	printf("\n");
	printf("Raw HKW scenarios (output):\n");
	printf("stoch_rawscen_oneSrcOneSnk_12oct2026-17h09_gstd.dat\n");
	printf("\n");
	printf("Stacked scenarios of the MCFNDP instance in GSTD fomat (output):\n");
	printf("stoch_oneSrcOneSnk_12oct2026-17h09_gstd.txt\n");
	printf("\n");
	printf("8.1.3.2 Origin the base deterministic network\n");
	printf("*********************************************\n");
	printf("\n");
	printf("Similar to 8.1.1.2.\n");
	printf("\n");
	printf("8.2 Multiple but equal sources and multiple but equal sinks\n");
	printf("***********************************************************\n");
	printf("\n");
	printf("8.2.1 Folder /examples/eqSrcsEqSnks/GSTD/\n");
	printf("******************************************\n");
	printf("\n");
	printf("In the base network, there are multiple sources and multiple sinks for each\n");
	printf("commodity. However, sources and sinks are identical across commodities.\n");
	printf("There is differentiation of variable costs and capacities according to\n");
	printf("commodities. The base network must be received in the STD format (as the DOW\n");
	printf("format is unsuitable in this situation) and the generated networks must be\n");
	printf("saved in the GSTD format.\n");
	printf("\n");
	printf("8.2.1.1 Operation of the stochastic generator\n");
	printf("*********************************************\n");
	printf("\n");
	printf("Linux command line instruction:\n");
	printf("\n");
	printf("./stoch_gen +v +F stoch_config_eqSrcsEqSnks_12oct2026-17h09_gstd.par\n");
	printf("\n");
	printf("Configuration file (input) used with Linux command line:\n");
	printf("stoch_config_eqSrcsEqSnks_12oct2026-17h09_gstd.par\n");
	printf("\n");
	printf("Docker instruction:\n");
	printf("\n");
	printf("docker run --rm --mount type=bind,src=absolutePathToDockerFolder/inout,dst=/inout stoch_gen_image \\n");
	printf("/stoch_gen +v +F /inout/docker_stoch_config_eqSrcsEqSnks_12oct2026-17h09_gstd.par\n");
	printf("\n");
	printf("Configuration file (input) used with Docker instruction:\n");
	printf("docker_stoch_config_eqSrcsEqSnks_12oct2026-17h09_gstd.par\n");
	printf("\n");
	printf("Base network file (input):\n");
	printf("determ_eqSrcsEqSnks_12oct2026-17h09.std\n");
	printf("\n");
	printf("Generated moments (output):\n");
	printf("stoch_moms_eqSrcsEqSnks_12oct2026-17h09_gstd.dat\n");
	printf("\n");
	printf("Generated correlations (output):\n");
	printf("stoch_corrs_eqSrcsEqSnks_12oct2026-17h09_gstd.dat\n");
	printf("\n");
	printf("Raw HKW scenarios (output):\n");
	printf("stoch_rawscen_eqSrcsEqSnks_12oct2026-17h09_gstd.dat\n");
	printf("\n");
	printf("Stacked scenarios of the MCFNDP instance in GSTD fomat (output):\n");
	printf("stoch_eqSrcsEqSnks_12oct2026-17h09_gstd.txt\n");
	printf("\n");
	printf("8.2.1.2 Origin of the base deterministic network\n");
	printf("************************************************\n");
	printf("\n");
	printf("The base deterministic network was originally produced by the determistic\n");
	printf("generator as follows.\n");
	printf("\n");
	printf("Linux command line instruction:\n");
	printf("\n");
	printf("./determ_gen +v +F determ_config_eqSrcsEqSnks_12oct2026-17h09.par\n");
	printf("\n");
	printf("Configuration file:\n");
	printf("determ_config_eqSrcsEqSnks_12oct2026-17h09.par\n");
	printf("\n");
	printf("Outputs in STD, LP, MPS, GRA (basic graph), TEX (code of graphical\n");
	printf("representation), pdf (compiled image) formats:\n");
	printf("determ_eqSrcsEqSnks_12oct2026-17h09.std\n");
	printf("determ_eqSrcsEqSnks_12oct2026-17h09.mps\n");
	printf("determ_eqSrcsEqSnks_12oct2026-17h09.lp\n");
	printf("determ_eqSrcsEqSnks_12oct2026-17h09.gra\n");
	printf("determ_eqSrcsEqSnks_12oct2026-17h09.tex\n");
	printf("determ_eqSrcsEqSnks_12oct2026-17h09.pdf\n");
	exit(0);
}


