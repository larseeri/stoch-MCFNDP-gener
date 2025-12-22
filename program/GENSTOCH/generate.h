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
#ifndef __GENERATE_H__
#define __GENERATE_H__

#include <cmath>
#include <fstream>
#include <iostream>

#include "../HOYLAND_KAUT_WALLACE/HKW_sg.h" // from Hoyland-Kaut-Wallace
#include "../HOYLAND_KAUT_WALLACE/matrix.h" // from Hoyland-Kaut-Wallace
#include "configuration.h"
#include "random.h"

using namespace std;

// Inheriting from TMatrix (scen-gen_HKW_src/matrix.h). Simplifies I/O.
struct TMatrix2 : public TMatrix
{

    void Export(string filename)
    {
        ofstream os(filename);
		
		if (!os.is_open())
        {
            cerr << "\nERROR:" << endl;
            cerr << "The specified output file " << filename << " could not be opened." << endl;
            cerr << "EXECUTION ABORTED\n" << endl;
            exit(1);
        }
		
        os << *this;
    }

    void Import(string filename)
    {

        // checking if input file exists
        if (!filesystem::is_regular_file(filename))
        {
            cerr << "\nERROR:" << endl;
            cerr << "The specified input file " << filename << " could not be found." << endl;
            cerr << "EXECUTION ABORTED\n" << endl;
            exit(1);
        }

        ifstream is(filename);

        if (!is.is_open())
        {
            cerr << "\nERROR:" << endl;
            cerr << "The specified input file " << filename << " could not be opened." << endl;
            cerr << "EXECUTION ABORTED\n" << endl;
            exit(1);
        }

        // features better diagnostics than is >> *this;
        Mat_GetFromFile(this, filename.c_str());

        // is >> *this;
    }

    TMatrix2(int nrow, int ncol, double **val) : TMatrix{nrow, ncol, val} {};

    friend istream &operator>>(istream &is, TMatrix2 &t);

    friend ostream &operator<<(ostream &os, const TMatrix2 &t);
};

inline istream &operator>>(istream &is, TMatrix2 &t)
{
    int nrow, ncol;

    is >> nrow >> ncol;
    Mat_Init(&t, nrow, ncol);

    for (int i = 0; i < nrow; i++)
        for (int j = 0; j < ncol; j++)
            is >> t.val[i][j];

    return is;
}

inline ostream &operator<<(ostream &os, const TMatrix2 &t)
{
    os << t.nrow << endl << t.ncol << endl;

    for (int i = 0; i < t.nrow; i++)
    {
        for (int j = 0; j < t.ncol; j++)
            os << t.val[i][j] << " ";
        os << endl;
    }
    os << endl;

    return os;
}

// Handling generation of scenarios. Launches Hoyland-Kaut-Wallace code.
struct HoylandKautWallace
{
    const Configuration *params;

    const VectElem stochElems;

    // matrix of first four target moments (in format expected by
    // Hoyland-Kaut-Wallace code) for each stochastic element: expectation,
    //  variance, standardized skewness, excess standardized kurtosis
    TMatrix2 moments;

    // matrix of target correlation coefficients (in format expected by
    // Hoyland-Kaut-Wallace code) for all stochastic elements
    TMatrix2 correlation;

    // output matrix of scenarios of size (number of stochastic elements x
    // number of scenarios) possibly initialized from file and updated by
    // Hoyland-Kaut-Wallace code
    TMatrix2 outMat;

    // Constructor
    HoylandKautWallace(int numbMoments, int numbElems, const Configuration *_params, const VectElem &_stochElems);

    // Destructor
    ~HoylandKautWallace();

    // Generates target moments based on program input parameters.
    void generMoments(int numbMoments, int numbElems);

    // Generates the matrix of first four target moments of each stochastic element
    // (in format expected by Hoyland-Kaut-Wallace code)
    // when uniform(a, b) distributional characteristics are specified.
    void generMomentsUnifDist(int numbElems);

    // Generates the matrix of first four target moments of each stochastic element
    // (in format expected by Hoyland-Kaut-Wallace code)
    // when triangular(a, b, c) distributional characteristics are specified.
    void generMomentsTriangDist(int numbElems);

    // Generates matrix of target correlation coefficients based on values
    // specified for each block.
    void generCorrelations(int numbElems);

    // Generates scenarios using Hoyland-Kaut-Wallace code.
    TMatrix2 generScenarios();
};

#endif
