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
#include "generate.h"

// from random.h
Random rnd;

// Constructor
HoylandKautWallace::HoylandKautWallace(int numbMoments, int numbElems, const Configuration *_params,
                                       const VectElem &_stochElems)
    : params(_params), stochElems(_stochElems), moments(Mat_0), correlation(Mat_0), outMat(Mat_0)
{

    // moments are generated based on assumptions
    if (params->generMoments)
    {
        generMoments(numbMoments, numbElems);
        moments.Export(params->momentsFileName);
    }
    // moments are read from file
    else
    {
        moments.Import(params->momentsFileName);
        if (moments.nrow != 4 || moments.ncol != numbElems)
        {
            cerr << "\nERROR:\n"
                 << "While reading moments from file " << params->momentsFileName << "\n"
                 << "expected 4 rows and " << numbElems << " columns.\n"
                 << "Instead, read " << moments.nrow << " rows and " << moments.ncol << " columns.\n"
                 << "EXECUTION ABORTED\n"
                 << endl;
            exit(1);
        }
    }

    // correlations are generated based on assumptions
    if (params->specCorrBlocks)
    {
        generCorrelations(numbElems);
        correlation.Export(params->correlsFileName);
    }
    // correlations are read from file
    else
    {
        correlation.Import(params->correlsFileName);
        if (correlation.nrow != numbElems || correlation.ncol != numbElems)
        {
            cerr << "\nERROR:\n"
                 << "While reading correlations from file " << params->correlsFileName << "\n"
                 << "expected " << numbElems << " rows " << numbElems << " columns.\n"
                 << "Instead, read " << correlation.nrow << " rows and " << correlation.ncol << " columns.\n"
                 << "EXECUTION ABORTED\n"
                 << endl;
            exit(1);
        }
    }
}

// Destructor
HoylandKautWallace::~HoylandKautWallace()
{
    Mat_Kill(&moments);
    Mat_Kill(&correlation);
    // outMat is destroyed by ~Instance() since it must first be returned by GenerScenarios()
}

// Generates target moments based on specifications.
void HoylandKautWallace::generMoments(int numbMoments, int numbElems)
{
    Mat_Init(&moments, numbMoments, numbElems);

    switch (params->momentsCharacts)
    {
    case UNIFORM:
        generMomentsUnifDist(numbElems);
        break;
    case TRIANGULAR:
        generMomentsTriangDist(numbElems);
        break;
    default:
        break;
    }
}

// Generates the matrix of first four target moments of each stochastic element
// (in format expected by HKW code) when uniform(a, b) distributional
// characteristics are specified.
//
// order 1 moment (expectation): (1/2)(a + b),
// order 2 moment (variance): (1/12)(b - a)**2,
// order 3 moment (standardized skewness): 0,
// order 4 moment (excess standardized kurtosis): -(6/5),
// where,
// a = D - (params->alpha * D),
// b = D + (params->beta * D),
// D is the base value of the stochastic element taken from the base deterministic instance,
// params->alpha in [0, 1),
// params->beta in [0, \infty).
// The distribution is symmetric around D when params->alpha = params->beta.
void HoylandKautWallace::generMomentsUnifDist(int numbElems)
{
    double cst;

    // order 1 moment
    // m1 = (a + b) / 2  ==>  m1 = D * (2 - alpha + beta) / 2

    cst = (2 - params->alpha + params->beta) / 2.0;

    for (int i = 0; i < numbElems; i++)
        moments.val[0][i] = stochElems[ALL_TYPES][i] * cst;

    // order 2 moment
    // m2 = (b - a)^2 / 12  ==>  m2 = D^2 * (alpha + beta)^2 / 12

    cst = (params->alpha + params->beta) * (params->alpha + params->beta) / 12.0;

    for (int i = 0; i < numbElems; i++)
        moments.val[1][i] = stochElems[ALL_TYPES][i] * stochElems[ALL_TYPES][i] * cst;

    // order 3 moment

    for (int i = 0; i < numbElems; i++)
        moments.val[2][i] = 0.0;

    // order 4 moment

    for (int i = 0; i < numbElems; i++)
        moments.val[3][i] = -1.2;
}

// Generates the matrix of first four target moments of each stochastic element
// (in format expected by HKW code) when triangular(a, b, c) distributional
// characteristics are specified.
//
// order 1 moment (expectation): (1/3)(a + b + c),
// order 2 moment (variance): (1/18)(a**2 + b**2 + c**2 - a*b - a*c - b*c),
// order 3 moment (standardized skewness): 0, if \alpha = \beta
// otherwise, see https://mathworld.wolfram.com/TriangularDistribution.html equation (9),
// order 4 moment (excess standardized kurtosis): -(3/5),
// where,
// a = D - (params->alpha * D),
// b = D + (params->beta * D),
// c = D,
// D is the base value of the stochastic element taken from the base deterministic instance,
// params->alpha in [0, 1),
// params->beta in [0, \infty).
// The distribution is symmetric around D when params->alpha = params->beta.
void HoylandKautWallace::generMomentsTriangDist(int numbElems)
{
    double cst;

    // order 1 moment
    // m1 = (a + b + c) / 3  ==>  m = D * (3 - alpha + beta) / 3

    cst = (3 - params->alpha + params->beta) / 3.0;

    for (int i = 0; i < numbElems; i++)
        moments.val[0][i] = stochElems[ALL_TYPES][i] * cst;

    // order 2 moment
    // m2 = (a^2 + b^2 + c^2 - a*b - a*c - b*c) / 18
    // ==>  m2 = D^2 * (alpha^2 + beta^2 + (alpha * beta)) / 18
    // ==>  m2 = D^2 * [(alpha + beta)^2 - (alpha * beta)] / 18

    cst = ((params->alpha + params->beta) * (params->alpha + params->beta) - (params->alpha * params->beta)) / 18.0;

    for (int i = 0; i < numbElems; i++)
        moments.val[1][i] = stochElems[ALL_TYPES][i] * stochElems[ALL_TYPES][i] * cst;

    // order 3 moment
    // m3 = (2^0.5 * (a+b-2c) (2a-b-c) (a-2b+c)) / (5 * (a^2 + b^2 + c^2 - a*b - a*c - b*c)^1.5)
    // ==>  m3 = [(2^0.5) / 5] * [(beta - alpha) (-2 * alpha - beta) (-alpha - 2 * beta)] / [(alpha + beta)^2 - (alpha *
    // beta)]^1.5

    cst = pow(2, 0.5) / 5.0;

    double cst1 =
        (params->beta - params->alpha) * (-2 * params->alpha - params->beta) * (-2 * params->beta - params->alpha);
    double cst2 = ((params->alpha + params->beta) * (params->alpha + params->beta) - (params->alpha * params->beta));

    cst2 = pow(cst2, 1.5);

    cst *= cst1 / cst2;

    for (int i = 0; i < numbElems; i++)
        moments.val[2][i] = cst;

    // order 4 moment

    for (int i = 0; i < numbElems; i++)
        moments.val[3][i] = -0.6;
}

// Generates matrix of target correlation coefficients based on values
// specified for each block.
void HoylandKautWallace::generCorrelations(int numbElems)
{
    int count1 = 0, count2 = 0;

    Mat_Init(&correlation, numbElems, numbElems);

    for (uint i = 0; i < params->stochTypes.size(); i++)
    {
        for (uint ii = 0; ii < stochElems[params->stochTypes[i]].size; ii++)
        {
            for (uint j = 0; j < params->stochTypes.size(); j++)
            {
                for (uint jj = 0; jj < stochElems[params->stochTypes[j]].size; jj++)
                {
                    correlation.val[count1][count2] =
                        (count1 == count2) ? 1.0
                                           : params->correlCoeffsMatrix[params->stochTypes[i]][params->stochTypes[j]];
                    count2++;
                }
            }
            count1++;
            count2 = 0;
        }
    }
}

// Generates scenarios using Hoyland-Kaut-Wallace code.
TMatrix2 HoylandKautWallace::generScenarios()
{
    TVector probs = Vec_0;
    int useStartDistrib = 0;
    int errorlevel = 0;

    // define output matrix of scenarios
    Mat_Init(&outMat, moments.ncol, params->numbScenarios);

    // initialize output matrix of scenarios from file if specified
    if (params->readStartScenarios)
    {
        Mat_GetFromFile(&outMat, params->startScenariosFileName.c_str());
        useStartDistrib = 1;
    }

    // define vector of probabilities
    Vec_Init(&probs, params->numbScenarios);

    // initialize vector of probabilities from file if specified
    if (params->readProbs)
        Vec_GetFromFile(&probs, params->probsFileName.c_str());
    // otherwise, scenarios are equiprobable
    else
        for (int i = 0; i < params->numbScenarios; i++)
            probs.val[i] = 1.0 / params->numbScenarios;

    rnd.setSeed(params->randomSeed, params->randomStream);

    // launch generation of scenarios with Hoyland-Kaut-Wallace algorithm
    errorlevel = HKW_ScenGen(params->distCharactMoments, &moments, &correlation, &probs, &outMat,
                             params->maxErrorMoment, params->maxErrorCorrel, params->displayLevel, params->maxTrial,
                             params->maxIter, useStartDistrib, NULL, NULL, NULL, NULL);

    if (errorlevel > 0)
    {
		printf("\n");
		printf("*****************************************************************************\n");
		printf("\n");
        printf("ERROR: The Hoyland-Kaut-Wallace algorithm FAILED to satisfy the specified\n");
        printf("target moments and target correlations WITH THE REQUESTED NUMBER OF SCENARIOS\n");
        printf("WITHIN THE SPECIFIED MAXIMUM NUMBERS OF TRIALS AND ITERATION.\n");
		printf("\n");
        printf("The FOLLOWING COURSE OF ACTION IS RECOMMENDED:\n");
		printf("(reproduced from Section 6 of user guide)\n");
        printf("\n");
        printf("1- Verify that target moments and correlations are proper.\n");
		printf("\n");
		printf("If option G has not been selected, target moments are read from a file. It\n");
		printf("should be ascertained that they are compatible with a probability \n");
		printf("distribution. If this condition fails, convergence of the HKW algorithm is\n");
		printf("expected to fail. By construction, target moments that are generated when\n");
		printf("option G is selected are correctly specified.\n");
		printf("\n");
		printf("If option K has not been selected, a target correlation matrix is read from\n");
		printf("a file. it should be ascertained that this matrix is well-specified \n");
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
        exit(1);
    }

    if (params->displayLevel > 4)
        Mat_DisplayTransp(&outMat, "Scenarios");

    Vec_Kill(&probs);

    outMat.Export(params->outScenariosFileName);

    return outMat;
}
