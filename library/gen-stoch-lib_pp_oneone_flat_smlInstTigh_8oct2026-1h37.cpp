
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <regex>
#include <chrono>
#include <filesystem>

// Generates library of deterministic MCFNDPs
// Point-to-point primary topology, flat secondary topology, one source and one sinks for each commodity
// Handles generation of configuration files for a selection of small and tight instances (1.1)

using namespace std;
namespace fs = filesystem;

const string INROOT = "/local2/lareri/library/pointToPoint/oneSrcOneSnk/flat/";

const string OUTROOT = "./pointToPoint/oneSrcOneSnk/flat/smlInstTigh/";

const string CONFIGROOT = "/local2/lareri/library/pointToPoint/oneSrcOneSnk/flat/smlInstTigh/";

const int seedNumber = 50163844;

// pp:1, eqeq:1, flat:3, small:1, tight:1
const int initStreamNumber = 11311;

// A base deterministic instances

// 1 small instances

// 1.1 small and tight

vector<string> smInstTi = {"pp_vols-20-100_csts-20-100-20-100_oneSrcOneSnk_secTop-flat_secTopSpe-2-3-12-15-45-10-15-2_noPar-20_adjSpe-1~00-1~00_seed-50163844_strm-6313.std", "pp_vols-20-100_csts-20-100-20-100_oneSrcOneSnk_secTop-flat_secTopSpe-4-5-10-15-45-10-15-2_noPar-20_adjSpe-1~00-1~00_seed-50163844_strm-8813.std"};

vector<string> smInstTi_seedStream = {"base-det-seed-50163844_base-det-strm-6313", "base-det-seed-50163844_base-det-strm-8813"};

// 1.2 small and loose

vector<string> smInstLo = {"pp_vols-20-100_csts-20-100-20-100_oneSrcOneSnk_secTop-flat_secTopSpe-2-3-12-15-45-10-15-2_noPar-20_adjSpe-0~33-1~00_seed-50163844_strm-5313.std", "pp_vols-20-100_csts-20-100-20-100_oneSrcOneSnk_secTop-flat_secTopSpe-4-5-10-15-45-10-15-2_noPar-20_adjSpe-0~33-1~00_seed-50163844_strm-7813.stf"};

vector<string> smInstLo_seedStream = {"base-det-seed-50163844_base-det-strm-5313", "base-det-seed-50163844_base-det-strm-7813"};

// 2 large instances

// 2.1 large and tight

vector<string> laInstTi = {"pp_vols-20-100_csts-20-100-20-100_oneSrcOneSnk_secTop-flat_secTopSpe-6-7-21-15-45-10-15-2_noPar-20_adjSpe-1~00-1~00_seed-50163844_strm-16313.std", "pp_vols-20-100_csts-20-100-20-100_oneSrcOneSnk_secTop-flat_secTopSpe-8-9-36-15-45-10-15-2_noPar-20_adjSpe-1~00-1~00_seed-50163844_strm-23813.std"};

vector<string> laInstTi_seedStream = {"base-det-seed-50163844_base-det-strm-16313", "base-det-seed-50163844_base-det-strm-23813"};

// 2.2 large and loose

vector<string> laInstLo = {"pp_vols-20-100_csts-20-100-20-100_oneSrcOneSnk_secTop-flat_secTopSpe-6-7-21-15-45-10-15-2_noPar-20_adjSpe-0~33-1~00_seed-50163844_strm-15313.std", "pp_vols-20-100_csts-20-100-20-100_oneSrcOneSnk_secTop-flat_secTopSpe-8-9-36-15-45-10-15-2_noPar-20_adjSpe-0~33-1~00_seed-50163844_strm-22813.std"};

vector<string> laInstLo_seedStream = {"base-det-seed-50163844_base-det-strm-15313", "base-det-seed-50163844_base-det-strm-22813"};

// B option key-values 

// 1 small instances

vector<string> smInst_optVal_T = {"-T T", "-T U"};
vector<string> smInst_optName_T = {"T-T", "T-U"};

vector<string> smInst_optVal_n = {"-n 1500", "-n 3000"};
vector<string> smInst_optName_n = {"n-1500", "n-3000"};

vector<string> smInst_optVal_m = {"-m 0.001"};
vector<string> smInst_optName_m = {"m-0~001"};

vector<string> smInst_optVal_c = {"-c 0.001"};
vector<string> smInst_optName_c = {"c-0~001"};

vector<string> smInst_optVal_SKX = {"// randomizing volumes\n-S 1\n\n// correlations are hereafter specified for each block\n-K\n-XDD 0.5", "// randomizing overall capacities\n-S 2\n// correlations are hereafter specified for each block\n-K\n-XAA 0.5", "// randomizing volumes and overall capacities\n-S 3\n// correlations are hereafter specified for each block\n-K\n-XDD 0.5\n-XDA -0.2\n-XAA 0.5", "// randomizing volumes and commodity-specific capacities\n-S 5\n// correlations are hereafter specified for each block\n-K\n-XDD 0.5\n-XDC -0.2\n-XCC 0.5"};
vector<string> smInst_optName_SKX = {"S-1_K_XDD-0~5", "S-2_K_XAA-0~5", "S-3_K_XDD-0~5_XDA-0~2_XAA-0~5", "S-5_K_XDD-0~5_XDC--0~2_XCC-0~5"};

// 1.1 small and tight instances

vector<string> smInstTi_optVal_AB = {"// alpha\n-A 0.25\n\n// beta\n-B 0.25"};
vector<string> smInstTi_optName_AB = {"A-0~25_B-0~25"};

// 1.2 small and loose instances

vector<string> smInstLo_optVal_AB = {"// alpha\n-A 0.5\n\n// beta\n-B 0.5"};
vector<string> smInstLo_optName_AB = {"A-0~5_B-0~5"};

// 2 large instances

vector<string> laInst_optVal_T = {"-T T", "-T U"};
vector<string> laInst_optName_T = {"T-T", "T-U"};

vector<string> laInst_optVal_n = {"-n 2500", "-n 5000"};
vector<string> laInst_optName_n = {"n-2500", "n-5000"};

vector<string> laInst_optVal_m = {"-m 0.01"};
vector<string> laInst_optName_m = {"m-0~01"};

vector<string> laInst_optVal_c = {"-c 0.01"};
vector<string> laInst_optName_c = {"c-0~01"};

vector<string> laInst_optVal_SKX = {"// randomizing volumes\n-S 1\n\n// correlations are hereafter specified for each block\n-K\n-XDD 0.5", "// randomizing overall capacities\n-S 2\n// correlations are hereafter specified for each block\n-K\n-XAA 0.5", "// randomizing volumes and overall capacities\n-S 3\n// correlations are hereafter specified for each block\n-K\n-XDD 0.5\n-XDA -0.2\n-XAA 0.5"};
vector<string> laInst_optName_SKX = {"S-1_K_XDD-0~5", "S-2_K_XAA-0~5", "S-3_K_XDD-0~5_XDA-0~2_XAA-0~5"};

// 2.1 large and tight instances

vector<string> laInstTi_optVal_AB = {"// alpha\n-A 0.25\n// beta\n-B 0.25"};
vector<string> laInstTi_optName_AB = {"A-0~25_B-0~25"};

// 2.2 large and loose instances

vector<string> laInstLo_optVal_AB = {"// alpha\n-A 0.5\n// beta\n-B 0.5"};
vector<string> laInstLo_optName_AB = {"A-0~5_B-0~5"};


int main( int argc , char ** argv )
{		
	string folderPathStr, filePath;
	string instName, instName0, instName1, instName2, instName3, instName4, instName5;
	string optsStr, optsStr0, optsStr1, optsStr2, optsStr3, optsStr4, optsStr5;
	
	int streamNumber = initStreamNumber;	
	
	// 1.1 small and tight instances	
	
	for (int ind0 = 0; ind0 < smInstTi.size(); ind0++)
	{				
		optsStr0 = "// base deterministic instance file has GSTD format\n-F G\n\n";
		optsStr0 += ("// path of base deterministic instance\n-I " + INROOT + smInstTi[ind0] + "\n\n");
		
		instName0 = "pp_oneSrcOneSnk_flat_" + smInstTi_seedStream[ind0] + "_smlInstTigh_";
		instName0 += "F-G_";

		for (int ind1 = 0; ind1 < smInst_optVal_T.size(); ind1++)
		{			
			optsStr1 = ("// moments will be generated\nG\n\n");
			instName1 = ("G_");			
			
			optsStr1 += ("// distribution\n" + smInst_optVal_T[ind1] + "\n\n");
			instName1 += (smInst_optName_T[ind1] + "_");
			
			for (int ind2 = 0; ind2 < smInstTi_optVal_AB.size(); ind2++)
			{	
				optsStr2 = (smInstTi_optVal_AB[ind2] + "\n\n");
				instName2 = (smInstTi_optName_AB[ind2] + "_");	
		
				for (int ind3 = 0; ind3 < smInst_optVal_SKX.size(); ind3++)
				{				
					optsStr3 = (smInst_optVal_SKX[ind3] + "\n\n");
					instName3 = (smInst_optName_SKX[ind3] + "_");					
		
					for (int ind4 = 0; ind4 < smInst_optVal_n.size(); ind4++)
					{						
						optsStr4 = ("// number of scenarios\n" + smInst_optVal_n[ind4] + "\n\n");
						instName4 = (smInst_optName_n[ind4] + "_");	

						optsStr4 += ("// maximum number of trials in HKW algo.\n-t 20\n\n");
						instName4 += ("t-20_");
						
						optsStr4 += ("// maximum number of iterations per trial in HKW algo.\n-i 50\n\n");
						instName4 += ("i-50_");
						
						optsStr4 += ("// level of onscreen reporting\n-l 3\n\n");
						instName4 += ("l-3_");						
	
						for (int ind5 = 0; ind5 < smInst_optVal_m.size(); ind5++)
						{
							optsStr5 = ("// tolerance about errors in moments\n" + smInst_optVal_m[ind5] + "\n\n");
							instName5 = (smInst_optName_m[ind5] + "_");	
													
							optsStr5 += ("// tolerance about errors in correlations\n" + smInst_optVal_c[ind5] + "\n\n");
							instName5 += (smInst_optName_c[ind5] + "_");

							streamNumber += 100;
							instName5 += ("seed-" + to_string(seedNumber) + "_strm-" + to_string(streamNumber) + "_gstd");
														
							optsStr5 += ("// seed\n-r " + to_string(seedNumber) + "\n\n"); 							
							optsStr5 += ("// stream\n-R " + to_string(streamNumber) + "\n\n");

							instName = instName0 + instName1 + instName2 + instName3 + instName4 + instName5;

							optsStr5 += ("// path of file where correlations will be saved\n-C ");
							optsStr5 += (CONFIGROOT + "corrs_" + instName + ".dat\n\n");

							optsStr5 += ("// path of file where moments will be saved\n-M ");
							optsStr5 += (CONFIGROOT + "moms_" + instName + ".dat\n\n");
								
							optsStr5 += ("// path of file where raw scenarios will be saved\n-o ");
							optsStr5 += (CONFIGROOT + "raw-scens_" + instName + ".dat\n\n");
								
							optsStr5 += ("// path of file where generated stochastic instance will be saved\n-O ");
							optsStr5 += (CONFIGROOT + "stoch-inst_" + instName + ".txt\n\n");	
								
							optsStr = optsStr0 + optsStr1 + optsStr2 + optsStr3 + optsStr4 + optsStr5;
							
							optsStr = ("// this configuration file is called with this command line:\n// ./stoch_gen +v +F " + CONFIGROOT + instName + ".par\n\n" + optsStr);
								
							folderPathStr = OUTROOT;				
							fs::path folderPath = folderPathStr;
								
							try
							{        
								fs::create_directories(folderPathStr);            
							}
							catch (const fs::filesystem_error& e)
							{
								cerr << "Error: " << e.what() << '\n';
							}

							filePath = folderPathStr + instName + ".par";				
							ofstream outFile(filePath);				
							outFile << optsStr << endl;	

						}
					}
				}
			}			
		}
	}	
}


