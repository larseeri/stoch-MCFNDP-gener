/*********************************************************************/
/***  Front-end for the HKW scenario generation algorithm          ***/
/***                                                               ***/
/***  author: Michal Kaut                                          ***/
/*********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "../HOYLAND_KAUT_WALLACE/HKW_sg.h"
#include "../HOYLAND_KAUT_WALLACE/matrix.h"
#include "random.h"


// Parameters that can be changed from command line
static int MaxTrial = 10;
static int HKW_MaxIter = 20;
static int FormatOfTgMoms = 0; // see PrintUsageAndExit for details
static double MaxErrMom = 1e-3;
static double MaxErrCorr = 1e-3;
static int TestLevel = 2;
static char OutFileName[] = "out_scen.txt";
//
// The next six parameters are unused for the purposes of the
// stochastic generator; commented out to silence gcc.
// (EL: 11 Oct 2025) 
//// static int NmbScen = 4;
//// static char ProbFileName[] = "";
//// static char *StartFileName = NULL;
//// static unsigned char UseStartDistrib = 0;
//// static int RandomSeed = 122334;
//// static int StreamSeed = 455667;

double InMom[13]; /// Input moments: 12 moments + InMom[0]:=1
double TgMom[4];  /// Target moments

void PrintUsageAndExit([[maybe_unused]] char ExecName[])
{
  printf("\nScenario generation code based on paper by K. Høyland, M. Kaut & S.W. Wallace.\n");
  printf("Code by Michal Kaut (michal.kaut@iot.ntnu.no) & Diego Mathieu.\n");

  printf("\nUsage: %s nmb_scen [options]\n", ExecName);
  printf("\nList of options                                                 %15s\n", "default value");
  printf(" -t trials ... number of Trials (attempts to generate)          %15d\n", MaxTrial);
  printf(" -i iters  ... number of Iterations in HKW alg.                 %15d\n", HKW_MaxIter);
  printf(" -f format ... Format of input data - see below                 %15d\n", FormatOfTgMoms);
  printf(" -m mError ... maximal error for Moments (scaled to var=1)      %15g\n", MaxErrMom);
  printf(" -c cError ... maximal error for Correlations (scaled to var=1) %15g\n", MaxErrCorr);
  printf(" -l levOut ... output Level (amount of output)                  %15d\n", TestLevel);
  printf(" -o oFile  ... Output filename (for scenarios)                  %15s\n", OutFileName);
  printf(" -p pFile  ... read Probabilities from a given file             %15s\n", "p[s]=1/nmb_scen");
  printf(" -s sFile  ... read Starting distribution from a given file     %15s\n", "sampled values");
  printf(" -r rSeed  ... set a random seed                                %15s\n", "use comp. time");
  printf(" -R sSeed  ... set a (stream) random seed (pcg_random)          %15s\n", "use comp. time");
  printf(" -h        ... display an additional Help message\n");

  printf("\n%s:\n  %s\n  %s\n  %s\n  %s\n  %s\n",
	 "option -f: format is a sum of following bits (any number from 0 to 16)",
	 " 1 -> use population estimators (as in spreadsheets)",
	 " 2 -> 2nd moment is Var instead of StDev",
	 " 4 -> 4th moment is Kurtosis - 3",
	 " 8 -> Higher moments are not scaled by StDev",
	 "16 -> Use non-central moments, E{X^i} ... lower bits are ignored"
	 );
  exit(1);
}


void PrintHelpAndExit([[maybe_unused]] char ExecName[])
{
  printf("\n");
  printf("\nScenario generation code from paper 'A Heuristic for Moment-matching Scenario");
  printf("\nGeneration' by K. Høyland, M. Kaut & S.W. Wallace; Computational Optimization");
  printf("\nand Applications, 24 (2-3), pp. 169–185, 2003; doi:10.1023/A:1021853807313.\n");
  printf("Code by Michal Kaut (michal.kaut@iot.ntnu.no) & Diego Mathieu.\n");
  printf("\n");
  printf("The code generates scenarios for multivariate random variables.\n");
  printf("Distribution is described by the first four moments and correlations.\n");
  printf("Iterative algorithm -> reports distance from the target properties.\n");
  printf("No convergence guarantee -> runs in several trials, if needed.\n");

  printf("\nINPUT FILES\n");
  printf("target moments are by default in file named: 'tg_moms.txt'\n");
  printf("target correlations are by default in file: 'tg_corrs.txt'\n");
  printf("These files must include a matrix of numbers in the following format:\n");
  printf(" : number of rows\n");
  printf(" : number of columns\n");
  printf(" : data (by rows)\n");
  printf(" : rest of the file is ignored\n");
  printf("\n");
  printf("If the probabilities are given, the file has to be in a vector format:\n");
  printf(" : number of elements\n");
  printf(" : data\n");
  printf(" : rest of the file is ignored\n");
  printf("\n");
  exit(1);
}
