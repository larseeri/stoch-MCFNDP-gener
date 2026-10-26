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
#include "instance.h"

// Assigns sizes to arrays of stochastic and non-stochastic elements.
void Instance::assignSizes()
{
    stochElems[DEMAND].size = numbCommods;
    stochElems[CAP_ARC].size = numbArcs;
    stochElems[CAP_COMMOD].size = numbArcs * numbCommods;
    stochElems[FIXED_COST].size = numbArcs;
    stochElems[VAR_COST].size = numbArcs * numbCommods;
}

void displayTimeStamp()
{
	// Google AI
	auto now = std::chrono::system_clock::now();
    std::time_t time_now = std::chrono::system_clock::to_time_t(now);

    // convert to local time structure safely
    std::tm* local_time = std::localtime(&time_now);

    // output using stream manipulators
    std::cout << std::put_time(local_time, "%Y-%m-%d %H:%M:%S") << std::endl;
}

// Google
double roundToNPlaces(double value, int n)
{
    double multiplier = std::pow(10.0, n);
	
    // std::round(-3.145) becomes -3.15, correctly rounding away from zero
    return std::round(value * multiplier) / multiplier;
}

bool is_empty_or_whitespace(const string& str) {
    return str.empty() || all_of(str.begin(), str.end(), [](unsigned char c) {
        return isspace(c);
    });
}

bool isDigitsOnly(const string& str)
{   
    if (is_empty_or_whitespace(str))
		return false; 
    
    // cast to unsigned char prevents undefined behavior with negative values
    return all_of(str.begin(), str.end(), [](unsigned char c) { 
        return isdigit(c); 
    });
}

bool isRealNumber(const string& str) {
    try
	{
        size_t idx;
        // Attempt to convert string to double
        std::stod(str, &idx);
        
        // Ensure the entire string was consumed (rejects "3.14abc")
        return idx == str.length();
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

// UNUSED
void Instance::readProbs(string probsFileName, vector<double> &vectProbs)
{		
    if (!filesystem::is_regular_file(probsFileName))
    {
        cerr << "\nERROR:\n"
			<< "The probability input file " << probsFileName << " could not be found.\n"
			<< "EXECUTION ABORTED" << endl;
        exit(1);
    }

    ifstream inFile(probsFileName);

    if (!inFile.is_open())
    {
        cerr << "\nERROR:\n"
			<< "The probability input file " << probsFileName << " could not be opened.\n"
			<< "EXECUTION ABORTED" << endl;
        exit(1);
    }
	
	string numbStr;
	int numb;

	// reading the number of probabilities 
	inFile >> numbStr;

	if (isDigitsOnly(numbStr))
	{
		numb = stoi(numbStr);
			
		if ((numb < 1) || (numb > params->numbScenarios))
		{
			cerr << "\nERROR:\n"
				<< "While reading " + probsFileName + ", the token for the number of probabilities\n"
				<< "[" << numb << "]" << " was non-positive or larger than the number of scenarios " << params->numbScenarios << ".\n"
				<< "EXECUTION ABORTED" << endl;
			exit(1);
		}			
	}
	else
	{			
		cerr << "\nERROR:\n"
			<< "While reading " + probsFileName + ", the token for the number of probabilities\n"
			<< "[" << numbStr << "] could not be interpreted as a positive number.\n"
			<< "EXECUTION ABORTED" << endl;
        exit(1);
	}
		
	vectProbs.resize(params->numbScenarios, 0);
	
	string probStr;
	double prob;
	
	int countProbs = 0;
	double sumProbs = 0;
    while (inFile)
    {
       inFile >> probStr;
	   
		if (isRealNumber(probStr))
		{
			prob = stod(probStr);
				
			if ((prob > 1 ) || (prob <= 0))
			{
				cerr << "\nERROR:\n"
					<< "While reading " + probsFileName + ", probability token\n"
					<< "[" << probStr << "]" << " was non-positive or larger than 1.\n"
					<< "EXECUTION ABORTED" << endl;
				exit(1);
			}		
		}
		else
		{			
			cerr << "\nERROR:\n"
				<< "While reading " + probsFileName + ", the probability token\n"
				<< "[" << probStr << "] could not be interpreted as a real number.\n"
				<< "EXECUTION ABORTED" << endl;
			exit(1);
		}		

		countProbs++;
		
		if (countProbs > params->numbScenarios)
		{		
			cerr << "\nERROR:\n"
				<< "While reading " + probsFileName + ", the number of probability tokens\n"
				<< "exceeds the number of scenarios.\n"
				<< "EXECUTION ABORTED" << endl;
			exit(1);				
		}	
	
	sumProbs+=prob;	
	vectProbs.push_back(prob);
	}

	inFile.close();
  
	if (abs(sumProbs -1)  > 0.0001)
	{		
		cerr << "\nERROR:\n"
			<< "While reading " + probsFileName + ", the sum of probabilities\n"
			<< "differs from 1 by at least 0.0001.\n"
			<< "EXECUTION ABORTED" << endl;
        exit(1);
			
	}    
}

// Generates scenarios with HKW algorithm.
int Instance::getScenarios()
{
    printf("\nBegin generation of raw scenarios with HKW algorithm.\n");

	displayTimeStamp();	

    HoylandKautWallace hkw(numb_moments, stochElems[ALL_TYPES].size, params, stochElems);

    // matrix of HKW scenarios, one scenario per column, individual stochastic
    // elements on rows
    scenarios = hkw.generScenarios();
	
	// update vector of probabilities with that of hkw
	vectProbsPtr = hkw.vectProbsPtr;	
	
	displayTimeStamp();

    printf("End generation of raw scenarios with HKW algorithm.\n");

    return scenarios.ncol;
}

// Copies one scenario to array containing all stochastic elements for the
// purpose of verifying feasibility with solver. (Vectors cannot be copied
// directly since the resulting matrix would be improperly structured.)
void Instance::copyScenario(int idx)
{
    for (int i = 0; i < scenarios.nrow; i++) // index 1
        stochElems[ALL_TYPES].ptr[i] = scenarios.val[i][idx];
}

// Allocates arrays.
void Instance::allocate()
{
    int numbStochElem = 0;
    VectInt begin(MAX_TYPE, -1);

    stochElems.resize(ALL_TYPES + 1);

    for (int i = 0; i < MAX_TYPE; i++)
        stochElems[i].type = static_cast<STOCH_ELEMS_TYPE>(i);

    stochElems[ALL_TYPES].type = ALL_TYPES;
    assignSizes();

    // setting the beginning indexes of subarrays in master array of all
    // stochastic elements
    for (uint i = 0; i < params->stochTypes.size(); i++)
    {
        STOCH_ELEMS_TYPE elem = params->stochTypes[i];
        // stochastic elements of all types are stored in a single array
        // will be starting index of corrresponding subarray in master array
        // containing all stochastic elements
        begin[elem] = numbStochElem;
        numbStochElem += getSize(elem);
    }

    // non-stochastic elements are stored in separate arrays and therefore
    // each start from index 0
    for (uint i = 0; i < params->nonStochTypes.size(); i++)
    {
        STOCH_ELEMS_TYPE elem = params->nonStochTypes[i];
        begin[elem] = 0;
    }

    // setting up master array containing subarrays for all stochastic types
    stochElems[ALL_TYPES].ptr = (double *)malloc(numbStochElem * sizeof(double));
    stochElems[ALL_TYPES].size = numbStochElem;

    // setting up individual arrays for non-stochastic types
    for (uint i = 0; i < params->nonStochTypes.size(); i++)
    {
        STOCH_ELEMS_TYPE elem = params->nonStochTypes[i];
        // there is a distinct array for each non-stochastic type
        stochElems[elem].ptr = (double *)malloc(stochElems[elem].size * sizeof(double));
    }

    // locating subarrays of individual stochastic types inside
    // master array of all stochastic types
    for (uint i = 0; i < params->stochTypes.size(); i++)
    {
        STOCH_ELEMS_TYPE elem = params->stochTypes[i];
        // locate beginning of subarrays for elements of stochastic types
        // within master array
        stochElems[elem].ptr = &stochElems[ALL_TYPES][begin[elem]];
    }	
	
	// initializing commodity-specific capacities to zero since
	// zero-valued figures are excluded from STD format and therefore
	// not read by InstanceRStd::read() and InstanceGStd::read(). 
	int countComm = 0;
	for (int a = 0; a < numbArcs; a++)
	{
		for (int k = 0; k < numbCommods; k++)
        {
			stochElems[CAP_COMMOD][countComm] = 0;
            countComm++;
		}		
	}

    // setting correspondences between pointer defined in class and
    // beginning of corresponding subarray (stochastic type)
    // or array (non-stochastic type)
    Arc::setPtrCapac(stochElems[CAP_ARC].ptr);
    Arc::setPtrCost(stochElems[FIXED_COST].ptr);
    Commodity::setPtr(stochElems[DEMAND].ptr);
    CommodOnArc::setPtrCapac(stochElems[CAP_COMMOD].ptr);
    CommodOnArc::setPtrCost(stochElems[VAR_COST].ptr);

    arcs.resize(numbArcs);
    commods.resize(numbCommods);
    nodes.resize(numbNodes);
}

// Adjusts final commodity flows to restore balance between total
// supplied and demanded volumes when there are several sources
// and/or sinks.
void InstanceGStd::adjCommodFlows()
{	
	const double margin = 1;
	
	for (int k = 0; k < numbCommods; k++)
    {		
		double totalOutflow = 0;
		double totalInflow = 0;				
        for (int n = 0; n < numbNodes; n++)
        {			
			// all non-negative
			commods[k].adjOutflows[n] = commods[k].outflows[n] * (*(commods[k].volume) / commods[k].origVolume);	

			commods[k].adjOutflows[n] = roundToNPlaces(commods[k].adjOutflows[n], DECIMAL_PLACES);
				
			totalOutflow += commods[k].adjOutflows[n];
					   
			// all non-positive
            commods[k].adjInflows[n] = commods[k].inflows[n] * (*(commods[k].volume) / commods[k].origVolume);			
			
			commods[k].adjInflows[n] = roundToNPlaces(commods[k].adjInflows[n], DECIMAL_PLACES);
			
			totalInflow += commods[k].adjInflows[n];
		}
		
		double netFlow = totalOutflow + totalInflow;		
		
		if (netFlow > 0)
		{
			int nOut = 0;
			while (commods[k].adjOutflows[nOut] < netFlow + margin)
			{
				nOut++;			
			}			
			commods[k].adjOutflows[nOut] -= netFlow;			
		}
		else if (netFlow < 0)
		{
			int nIn = 0;
			while (commods[k].adjInflows[nIn] < netFlow - margin)
			{
				nIn++;			
			}			
			commods[k].adjInflows[nIn] += netFlow;			
		}		
    }		
}

// Verifies if at least one of the randomized values in a scenario is negative.
bool Instance::hasNegativeValue()
{		
	if (find(params->stochTypes.begin(), params->stochTypes.end(), DEMAND) != params->stochTypes.end())
	{
		for (int k = 0; k < numbCommods; k++)					
			if (stochElems[DEMAND][k] < 0)					
				return true;
	}
	if (find(params->stochTypes.begin(), params->stochTypes.end(), CAP_ARC) != params->stochTypes.end())
	{
		for (int a = 0; a < numbArcs; a++)		
			if (stochElems[CAP_ARC][a] < 0)
				return true;
	}
	if (find(params->stochTypes.begin(), params->stochTypes.end(), CAP_COMMOD) != params->stochTypes.end())
	{		
		int countComm = 0;
		for (int a = 0; a < numbArcs; a++)		
			for (int k = 0; k < numbCommods; k++)			
				if(stochElems[CAP_COMMOD][countComm++] < 0)
					return true;			
	}		
	if (find(params->stochTypes.begin(), params->stochTypes.end(), FIXED_COST) != params->stochTypes.end())
	{
			for (int a = 0; a < numbArcs; a++)		
			if (stochElems[FIXED_COST][a] < 0)
				return true;			
	}		
	if (find(params->stochTypes.begin(), params->stochTypes.end(), VAR_COST) != params->stochTypes.end())
	{
		int countComm = 0;
		for (int a = 0; a < numbArcs; a++)		
			for (int k = 0; k < numbCommods; k++)			
				if(stochElems[VAR_COST][countComm++] < 0)
					return true;			
	}		
	
	return false;
}

// Reads a base deterministic MCFND instance stored in DOW format.
void InstanceDow::read(istream &is)
{
	if (params->verbose)
		cout << "\n****** Begin reading DOW instance from " << params->inputInstFileName << "." << endl;	
	
    string line, dummy;
    double capac, volume, varCost, fixedCost;
    int countArc = 0, countDemand = 0, countComm = 0;

    getline(is, line);
	
	if (line != "MULTIGEN.DAT:")
	{	
		cerr << "\nERROR:\n"
			<< "The file " << params->inputInstFileName << " specified with option key I\n"
			<< "does not conform to the expected DOW format specified with option key F.\n"
			<< "Header is missing.\n"
			<< "EXECUTION ABORTED" << endl;
			exit(1);	
	}

    is >> numbNodes >> numbArcs >> numbCommods;
	
	if (numbNodes < 1 || numbArcs < 1 || numbCommods < 1)
	{
		cerr << "\nERROR:\n"
			<< "At least one among the number of nodes, number of arcs and number of commodities\n"
			<< "stated at beginning of instance file " << params->inputInstFileName << "\n"
			<< "is non-positive.\n"
			<< "EXECUTION ABORTED" << endl;
			exit(1);		
	}

    allocate();
	
    for (int a = 0; a < numbArcs; a++)
    {
        Arc arc(numbCommods);
        int idxArc;

        is >> arc.orig >> arc.dest >> varCost >> capac >> fixedCost >> dummy >> idxArc;			
			
		if (arc.orig > numbNodes || arc.dest > numbNodes || arc.orig < 1 || arc.dest < 1)
		{
			cerr << "\nERROR:\n"
				<< "At least one of the two node indices (" << arc.orig << ", " << arc.dest << ") stated\n"
				<< "in connection with the " << a+1 << "th arc listed (starting from 1) in instance file\n"
				<< params->inputInstFileName << " is non-positive or exceeds the total number of nodes stated\n"
				<< "at beginning of file.\n"
				<< "EXECUTION ABORTED" << endl;
				exit(1);		
		}
		if (varCost < 0 || capac < 0 || fixedCost < 0)
		{
			cerr << "\nERROR:\n"
				<< "At least one of the variable cost " << varCost << ", capacity " << capac << " and fixed cost\n"
				<< fixedCost << " stated in connection with the " << a+1 << "th arc listed (starting from 1)\n"
				<< "in instance file " << params->inputInstFileName << " is negative.\n"
				<< "EXECUTION ABORTED" << endl;
				exit(1);		
		}
		if (idxArc != (a + 1))
		{
			cerr << "\nERROR:\n"
				<< "The arc index (starting from 1) " << idxArc << " stated in connection with the\n"
				<< a+1 << "th arc listed (starting from 1) in instance file " << params->inputInstFileName << "\n"
				<< "does not agree with expected value " << a+1 << ".\n"
				<< "EXECUTION ABORTED" << endl;
				exit(1);		
		}
				
		if (capac < 0)
			capac = ASSIGNED_CAPAC_IF_READ_NEG;
		
        arc.orig--;
        arc.dest--;		
		
        stochElems[CAP_ARC][countArc] = capac;
        stochElems[FIXED_COST][countArc] = fixedCost;
        nodes[arc.orig].arcOut.push_back(countArc);
        nodes[arc.dest].arcIn.push_back(countArc);
        countArc++;
        idxArc--;

        for (int k = 0; k < numbCommods; k++)
        {
            stochElems[CAP_COMMOD][countComm] = capac;
            stochElems[VAR_COST][countComm] = varCost;
            countComm++;
        }

        arcs[idxArc] = arc;
    }
	
	if (countArc != numbArcs)
	{
		cerr << "\nERROR:\n"
			<< "The number of arcs " << countArc << " detected in instance file " << params->inputInstFileName << "\n"
			<< "does not agree with the number of arcs " << numbArcs << " stated at beginning of file.\n"
			<< "EXECUTION ABORTED" << endl;
			exit(1);		
	}

    for (int k = 0; k < numbCommods; k++)
    {
        is >> commods[k].orig >> commods[k].dest >> volume;
		
		if (commods[k].orig > numbNodes || commods[k].orig < 1 || commods[k].dest > numbNodes || commods[k].dest < 1)
		{	
			cerr << "\nERROR:\n"
				<< "At least one of the origin and destination node indices " << commods[k].orig << " and\n"
				<< commods[k].orig << " stated in connection with a commodity volume in instance file\n"
				<< params->inputInstFileName << " is either non-positive or exceeds the total\n"
				<< "number of nodes stated at beginning of file.\n"
				<< "EXECUTION ABORTED" << endl;
				exit(1);
		}

        commods[k].orig--;
        commods[k].dest--;
        stochElems[DEMAND][countDemand++] = volume;
    }
	
	if (params->verbose)
		cout << "\n****** End reading DOW instance from " << params->inputInstFileName << ".\n" << endl;
}

// Writes a generated MCFND instance to file in DOW format.
void InstanceDow::write(ostream &os, string feasibility, string validity) const
{
	//if (params->verbose)
	//	cout << "\n****** Begin writing DOW instance in " << params->inputInstFileName << ".\n" << endl;	
	
    static int scen = 0;
	os << endl;	
    os << "SCENARIO: " << scen << " Feasibility: " << feasibility << " Validity: " << validity << " Prob: " << (*vectProbsPtr)[scen] << endl;
	scen++;
	
    os << "MULTIGEN.DAT:" << endl;
    os << setw(15) << right << numbNodes << "   " << setw(15) << right << numbArcs << "   " << setw(15) << right
       << numbCommods << endl;
    
	for (int a = 0; a < numbArcs; a++)
	{
		os << fixed << setprecision(DECIMAL_PLACES) << setw(15) << right << arcs[a].orig + 1 << "   " << setw(15) << right
			<< arcs[a].dest + 1 << "   " << setw(15) << right << *(arcs[a].comms[0].cost) << "   " << setw(15)
			<< right << *(arcs[a].capac) << "   " << setw(15) << right << *(arcs[a].cost) << "   " << setw(15)
			<< right << 1 << "   " << setw(15) << right << a + 1 << endl;           
	}

    for (int k = 0; k < numbCommods; k++)
    {
        os << fixed << setprecision(DECIMAL_PLACES) << setw(15) << right << commods[k].orig + 1 << "   " << setw(15) << right
           << commods[k].dest + 1 << "   " << setw(15) << right << *(commods[k].volume) << endl;
    }
	
	//if (params->verbose)
	//	cout << "\n****** End writing DOW instance in " << params->inputInstFileName << ".\n" << endl;	
}

// Reads a deterministic MCFNDP instance stored in restricted STD format.
// Specialized to handle cases where there is a single source node
// and a single sink node for each commodity.
void InstanceRStd::read(istream &is)
{
	if (params->verbose)
		cout << "\n****** Begin reading RSTD instance " << params->inputInstFileName << "." << endl;	
		
    string line, dummy;
    int numbComm, idxComm, idxNode;
	double capac, volume, varCost, fixedCost;
    int countArc = 0, countComm = 0;

    is >> numbNodes >> numbArcs >> numbCommods;
	
	if (numbNodes < 1 || numbArcs < 1 || numbCommods < 1)
	{
		cerr << "\nERROR:\n"
			<< "At least one among the number of nodes, number of arcs and number of commodities\n"
			<< "stated at beginning of instance file " << params->inputInstFileName << "\n"
			<< "is non-positive.\n"
			<< "EXECUTION ABORTED" << endl;
			exit(1);		
	}	

    allocate();

    for (int a = 0; a < numbArcs; a++)
    {
        Arc arc(numbCommods);

        is >> arc.orig >> arc.dest >> fixedCost >> capac >> numbComm;
			
		if (arc.orig > numbNodes || arc.dest > numbNodes || arc.orig < 1 || arc.dest < 1)
		{
			cerr << "\nERROR:\n"
				<< "At least one of the two node indices (" << arc.orig << ", " << arc.dest << ") stated\n"
				<< "in connection with an arc listed in instance file " << params->inputInstFileName << "\n"
				<< "is non-positive or exceeds the total number of nodes stated at beginning of file.\n"
				<< "EXECUTION ABORTED" << endl;
				exit(1);		
		}
		if (fixedCost < 0)
		{
			cerr << "\nERROR:\n"
				<< "Fixed cost " << fixedCost << " stated in connection with the " << a+1 << "th arc row listed\n"
				<< "(starting from 1) in instance file " << params->inputInstFileName << " is negative.\n"
				<< "EXECUTION ABORTED" << endl;
				exit(1);		
		}
		if (numbCommods < numbComm)
		{
			cerr << "\nERROR:\n"
				<< "In instance file " << params->inputInstFileName << " the number of commodities " << numbComm << "\n"
				<< "stated in connection with arc (" << arc.orig << ", " << arc.dest << ") for which the\n"
				<< "capacities are non-zero exceeds the number of commodities stated at beginning of file.\n"
				<< "EXECUTION ABORTED" << endl;
				exit(1);		
		}	

        arc.orig--;
        arc.dest--;
		
		if (capac < 0)
			capac = ASSIGNED_CAPAC_IF_READ_NEG;
		
        stochElems[CAP_ARC][countArc] = capac;
        stochElems[FIXED_COST][countArc] = fixedCost;
        nodes[arc.orig].arcOut.push_back(countArc);
        nodes[arc.dest].arcIn.push_back(countArc);
        countArc++;		
		
        for (int k = 0; k < numbComm; k++)
        {
            is >> idxComm;
			
			if (idxComm > numbCommods)
			{
				cerr << "\nERROR:\n"
					<< "The commodity index " << idxComm << " stated while listing the costs and capacities\n"
					<< "in connection with arc ("<< (arc.orig + 1) << ", " << (arc.dest + 1) << ") in instance file\n"
					<< params->inputInstFileName << " exceeds the total number of commodities at beginning of file.\n"
					<< "EXECUTION ABORTED" << endl;
					exit(1);
			}
			
            idxComm--;
            is >> varCost >> capac;
			
			if (varCost < 0)
			{
				cerr << "\nERROR:\n"
					<< "The variable cost " << varCost << " stated in connection with the " << k+1 << "th commodity row\n"
					<< "listed (starting from 1) in instance file " << params->inputInstFileName << " is negative.\n"
					<< "EXECUTION ABORTED" << endl;
					exit(1);		
			}			
				
			if (capac < 0)
				capac = ASSIGNED_CAPAC_IF_READ_NEG;
		
            stochElems[CAP_COMMOD][a * numbCommods + (idxComm)] = capac;
            stochElems[VAR_COST][a * numbCommods + (idxComm)] = varCost;
            countComm++;
        }
		
        arcs[a] = arc;
    }
	
	if (countArc != numbArcs)
	{
		cerr << "\nERROR:\n"
			<< "The number of arc rows detected in instance file " << params->inputInstFileName << "\n"
			<< "does not agree with the number of arcs stated at beginning of file.\n"
			<< "EXECUTION ABORTED" << endl;
			exit(1);		
	}	

    for (int k = 0; k < numbCommods * 2; k++)
    {
        is >> idxComm >> idxNode >> volume;
		
		if (idxComm > numbCommods || idxComm < 1 || idxNode > numbNodes || idxNode < 1)
		{	
			cerr << "\nERROR:\n"
				<< "At least one of the commodity index " << idxComm << " and the node index\n"
				<< idxNode << " stated in connection with a commodity volume in instance file\n"
				<< params->inputInstFileName << " is either non-positive or exceeds the total\n"
				<< "number of commodities or nodes stated at beginning of file.\n"
				<< "EXECUTION ABORTED" << endl;
				exit(1);
		}

        idxComm--;
        idxNode--;

        if (volume < 0)
            commods[idxComm].dest = idxNode;
        else
        {
            commods[idxComm].orig = idxNode;
            stochElems[DEMAND][idxComm] = volume;
        }
    }
	
	if (params->verbose)
		cout << "\n****** End reading RSTD instance " << params->inputInstFileName << ".\n" << endl;
}

// Writes a generated MCFND instance to file in restricted STD format.
//
// This is specialized to handle cases where there is a single source
// node and a single sink node for each commodity.
void InstanceRStd::write(ostream &os, string feasibility, string validity) const
{	
	//if (params->verbose)
	//	cout << "\n****** Begin writing RSTD instance in " << params->inputInstFileName << ".\n" << endl;	
	
    static int scen = 0;
	os << endl;
    os << "SCENARIO: " << scen << " Feasibility: " << feasibility << " Validity: " << validity << " Prob: " << (*vectProbsPtr)[scen] << endl;
    scen++;
	os << setw(15) << right << numbNodes << "   " << setw(15) << right << numbArcs << "   " << setw(15) << right
       << numbCommods << endl;

    for (int n = 0; n < numbNodes; n++)
    {
        for (int a = 0; a < numbArcs; a++)
        {
            if (arcs[a].orig == n)
            {
                int count = 0;

                for (int k = 0; k < numbCommods; k++)
                    if (*(arcs[a].comms[k].capac) != 0)
                        count++;

                os << fixed << setprecision(DECIMAL_PLACES) << setw(15) << right << arcs[a].orig + 1 << "   " << setw(15) << right
                   << arcs[a].dest + 1 << "   " << setw(15) << right << *(arcs[a].cost) << "   " << setw(15) << right
                   << *(arcs[a].capac) << "   " << setw(15) << right << count << endl;

                for (int k = 0; k < numbCommods; k++)
                    if (*(arcs[a].comms[k].capac) != 0)
                        os << fixed << setprecision(DECIMAL_PLACES) << setw(15) << right << k + 1 << "   " << setw(15) << right
                           << *(arcs[a].comms[k].cost) << "   " << setw(15) << right << *(arcs[a].comms[k].capac)
                           << endl;
            }
        }
    }

    for (int k = 0; k < numbCommods; k++)
    {
        os << fixed << setprecision(DECIMAL_PLACES) << setw(15) << right << k + 1 << "   " << setw(15) << right
           << commods[k].orig + 1 << "   " << setw(15) << right << *(commods[k].volume) << endl;
        os << fixed << setprecision(DECIMAL_PLACES) << setw(15) << right << k + 1 << "   " << setw(15) << right
           << commods[k].dest + 1 << "   " << setw(15) << right << -*(commods[k].volume) << endl;
    }
	
	//if (params->verbose)
	//	cout << "\n****** End writing RSTD instance in " << params->inputInstFileName << ".\n" << endl;	
}

// Reads a deterministic MCFNDP instance stored in generic STD format.
void InstanceGStd::read(istream &is)
{
	if (params->verbose)
		cout << "\n****** Begin reading GSTD instance " << params->inputInstFileName << "." << endl;
	
    string line, dummy;
    int numbComm, idxComm, idxNode;
	double capac, volume, varCost, fixedCost;
    int countArc = 0, countComm = 0;

	is >> numbNodes >> numbArcs >> numbCommods;

	if (numbNodes < 1 || numbArcs < 1 || numbCommods < 1)
	{
		cerr << "\nERROR:\n"
			<< "At least one among the number of nodes, number of arcs and number of commodities\n"
			<< "stated at beginning of instance file " << params->inputInstFileName << "\n"
			<< "is non-positive.\n"
			<< "EXECUTION ABORTED" << endl;
			exit(1);		
	}

    allocate();

    for (int a = 0; a < numbArcs; a++)
    {
        Arc arc(numbCommods);
		
		is >> arc.orig >> arc.dest >> fixedCost >> capac >> numbComm;
	
		if (arc.orig > numbNodes || arc.dest > numbNodes || arc.orig < 1 || arc.dest < 1)
		{	
			cerr << "\nERROR:\n"
				<< "At least one of the two node indices (" << arc.orig << ", " << arc.dest << ") stated\n"
				<< "in connection with an arc listed in instance file " << params->inputInstFileName << "\n"
				<< "is non-positive or exceeds the total number of nodes stated at beginning of file.\n"
				<< "EXECUTION ABORTED" << endl;
				exit(1);		
		}	
		if (fixedCost < 0)
		{
			cerr << "\nERROR:\n"
				<< "Fixed cost " << fixedCost << " stated in connection with the " << a+1 << "th arc row listed\n"
				<< "(starting from 1) in instance file " << params->inputInstFileName << " is negative.\n"
				<< "EXECUTION ABORTED" << endl;
				exit(1);		
		}
		if (numbCommods < numbComm)
		{
			cerr << "\nERROR:\n"
				<< "In instance file " << params->inputInstFileName << " the number of commodities " << numbComm << "\n"
				<< "stated in connection with arc (" << arc.orig << ", " << arc.dest << ") for which the\n"
				<< "capacities are non-zero exceeds the number of commodities stated at beginning of file.\n"
				<< "EXECUTION ABORTED" << endl;
				exit(1);		
		}	

        arc.orig--;
        arc.dest--;
		
		if (capac < 0)
			capac = ASSIGNED_CAPAC_IF_READ_NEG;
		
        stochElems[CAP_ARC][countArc] = capac;
        stochElems[FIXED_COST][countArc] = fixedCost;
        nodes[arc.orig].arcOut.push_back(countArc);
        nodes[arc.dest].arcIn.push_back(countArc);
        countArc++;

        for (int k = 0; k < numbComm; k++)
        {					
			is >> idxComm;
			
			if (idxComm > numbCommods)
			{
				cerr << "\nERROR:\n"
					<< "The commodity index " << idxComm << " stated while listing the costs and capacities\n"
					<< "in connection with arc ("<< (arc.orig + 1) << ", " << (arc.dest + 1) << ") in instance file\n"
					<< params->inputInstFileName << " exceeds the total number of commodities at beginning of file.\n"
					<< "EXECUTION ABORTED" << endl;
					exit(1);
			}
			
            idxComm--;			
						
			is >> varCost >> capac;			
				
			if (varCost < 0)
			{
				cerr << "\nERROR:\n"
					<< "The variable cost " << varCost << " stated in connection with the " << k+1 << "th commodity row\n"
					<< "listed (starting from 1) in instance file " << params->inputInstFileName << " is negative.\n"
					<< "EXECUTION ABORTED" << endl;
					exit(1);		
			}				
				
			if (capac < 0)
				capac = ASSIGNED_CAPAC_IF_READ_NEG;
            
			stochElems[CAP_COMMOD][a * numbCommods + (idxComm)] = capac;
            stochElems[VAR_COST][a * numbCommods + (idxComm)] = varCost;
			
            countComm++;
        }		

        arcs[a] = arc;
    }
	
	if (countArc != numbArcs)
	{
		cerr << "\nERROR:\n"
			<< "The number of arc rows detected in instance file " << params->inputInstFileName << "\n"
			<< "does not agree with the number of arcs stated at beginning of file.\n"
			<< "EXECUTION ABORTED" << endl;
			exit(1);		
	}	

    for (int k = 0; k < numbCommods; k++)
    {
        commods[k].outflows.resize(numbNodes, 0);
        commods[k].inflows.resize(numbNodes, 0);
		commods[k].adjOutflows.resize(numbNodes, 0);
        commods[k].adjInflows.resize(numbNodes, 0);
        commods[k].origVolume = 0;
    }

    while (is >> idxComm >> idxNode >> volume)
    {		
		if (idxComm > numbCommods || idxComm < 1 || idxNode > numbNodes || idxNode < 1)
		{	
			cerr << "\nERROR:\n"
				<< "At least one of the commodity index " << idxComm << " and the node index\n"
				<< idxNode << " stated in connection with a commodity volume in instance file\n"
				<< params->inputInstFileName << " is either non-positive or exceeds the total\n"
				<< "number of commodities or nodes stated at beginning of file.\n"
				<< "EXECUTION ABORTED" << endl;
				exit(1);
		}
		
        idxComm--;
        idxNode--;

        if (volume > 0)
        {
            commods[idxComm].outflows[idxNode] = volume;
            commods[idxComm].origVolume += volume;
        }
        else
            commods[idxComm].inflows[idxNode] = volume;
    }

    for (int k = 0; k < numbCommods; k++)
    {
        stochElems[DEMAND][k] = commods[k].origVolume;
    }
	
	if (params->verbose)
		cout << "\n****** End reading GSTD instance " << params->inputInstFileName << ".\n" << endl;
}

// Writes a generated MCFND instance to file in generic STD format.
void InstanceGStd::write(ostream &os, string feasibility, string validity) const
{
	//if (params->verbose)
	//	cout << "\n****** Begin writing GSTD instance in " << params->inputInstFileName << ".\n" << endl;	
	
    static int scen = 0;
	os << endl;
    os << "SCENARIO: " << scen << " Feasibility: " << feasibility << " Validity: " << validity << " Prob: " << (*vectProbsPtr)[scen-1] << endl;
    scen++;
	os << setw(15) << right << numbNodes << "   " << setw(15) << right << numbArcs << "   " << setw(15) << right
       << numbCommods << endl;

    for (int n = 0; n < numbNodes; n++)
    {
        for (int a = 0; a < numbArcs; a++)
        {
            if (arcs[a].orig == n)
            {
                int count = 0;

                for (int k = 0; k < numbCommods; k++)
                    if (*(arcs[a].comms[k].capac) != 0)
                        count++;

                os << fixed << setprecision(DECIMAL_PLACES) << setw(15) << right << arcs[a].orig + 1 << "   " << setw(15) << right
                   << arcs[a].dest + 1 << "   " << setw(15) << right << *(arcs[a].cost) << "   " << setw(15) << right
                   << *(arcs[a].capac) << "   " << setw(15) << right << count << endl;

                for (int k = 0; k < numbCommods; k++)
                    if (*(arcs[a].comms[k].capac) != 0)

                        os << fixed << setprecision(DECIMAL_PLACES) << setw(15) << right << k + 1 << "   " << setw(15) << right
                           << *(arcs[a].comms[k].cost) << "   " << setw(15) << right << *(arcs[a].comms[k].capac)
                           << endl;
            }
        }
    }

    for (int k = 0; k < numbCommods; k++)
    {
        for (int n = 0; n < numbNodes; n++)
        {

            if (commods[k].outflows[n] > 0)
                os << fixed << setprecision(DECIMAL_PLACES) << setw(15) << right << k + 1 << "   " << setw(15) << right << n + 1
                   << "   " << setw(15) << right
                   << commods[k].adjOutflows[n] << endl;
            if (commods[k].inflows[n] < 0)
                os << fixed << setprecision(DECIMAL_PLACES) << setw(15) << right << k + 1 << "   " << setw(15) << right << n + 1
                   << "   " << setw(15) << right
                   << commods[k].adjInflows[n] << endl;
        }
    }
	
	//if (params->verbose)
	//	cout << "\n****** End writing GSTD instance in " << params->inputInstFileName << ".\n" << endl;	
}
