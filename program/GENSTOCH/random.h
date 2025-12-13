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
#ifndef __RANDOM_H__
#define __RANDOM_H__

#include "pcg_random.hpp"
#include <random>

using namespace std;

typedef uint32_t SEED_TYPE;

// Handling random number generation by
// calling Permuted Congruential Generator:
// https://www.pcg-random.org/download.html
struct Random
{
    pcg32 rng;

    Random() : rng(pcg_extras::seed_seq_from<std::random_device>{})
    {
    }

    Random(SEED_TYPE _seed, SEED_TYPE _stream)
    {
        setSeed(_seed, _stream);
    }

    void setSeed(SEED_TYPE _seed, SEED_TYPE _stream)
    {
        rng.seed(_seed, _stream);
        rng(); // to initialize the generator properly
    }

    uint32_t getRandomInt()
    {
        return rng();
    }

    uint32_t getRandomInt(int upperBound)
    {
        return (upperBound > 0) ? rng(upperBound) : 0;
    }

    uint32_t getRandomInt(int lowerBound, int upperBound)
    {
        return getRandomInt(upperBound) + lowerBound;
    }

    double getRandomDoubl01()
    {
        return (double)rng() / rng.max();
    }
};

#endif
