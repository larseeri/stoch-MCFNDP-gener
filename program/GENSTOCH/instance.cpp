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

// Generates scenarios with HKW algorithm.
int Instance::getScenarios()
{

    printf("\nBegin generation of raw scenarios with HKW algorithm.\n");

    HoylandKautWallace hkw(numb_moments, stochElems[ALL_TYPES].size, params, stochElems);

    // matrix of HKW scenarios, one scenario per column, individual stochastic
    // elements on rows
    scenarios = hkw.generScenarios();

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

// Reads a base deterministic MCFND network stored in DOW format.
void InstanceDow::read(istream &is)
{

    string line, dummy;
    double capac, volume, varCost, fixedCost;
    int countArc = 0, countDemand = 0, countComm = 0;

    getline(is, line);

    is >> numbNodes >> numbArcs >> numbCommods;

    allocate();

    for (int a = 0; a < numbArcs; a++)
    {
        Arc arc(numbCommods);
        int idxArc;

        is >> arc.orig >> arc.dest >> varCost // arc.comms[0].cost
            >> capac >> fixedCost             // arc.cost
            >> dummy >> idxArc;

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

    for (int k = 0; k < numbCommods; k++)
    {
        is >> commods[k].orig >> commods[k].dest >> volume;

        commods[k].orig--;
        commods[k].dest--;
        stochElems[DEMAND][countDemand++] = volume;

        // TODO: EL: 29sept2025 doit probablement etre supprime; voir si applique dans SMS++
        // desactive le 7nov2025
        // for (int a = 0; a < numbArcs; a++)
        //    *(arcs[a].comms[k].capac) = min(*(arcs[a].comms[k].capac), volume);
    }
}

// Writes a generated MCFND network to file in DOW format.
void InstanceDow::write(ostream &os) const
{
    static int scen = 0;
    os << " SCENARIO        " << scen++ << endl;
    os << " MULTIGEN.DAT:" << endl;
    os << setw(12) << right << numbNodes << "   " << setw(12) << right << numbArcs << "   " << setw(12) << right
       << numbCommods << endl;

    for (int n = 0; n < numbNodes; n++)
        for (int a = 0; a < numbArcs; a++)
        {
            if (arcs[a].orig == n)
            {
                os << fixed << setprecision(3) << setw(12) << right << arcs[a].orig + 1 << "   " << setw(12) << right
                   << arcs[a].dest + 1 << "   " << setw(12) << right << *(arcs[a].comms[0].cost) << "   " << setw(12)
                   << right << *(arcs[a].capac) << "   " << setw(12) << right << *(arcs[a].cost) << "   " << setw(12)
                   << right << 1 << "   " << setw(12) << right << a + 1 << endl;
            }
        }

    for (int k = 0; k < numbCommods; k++)
    {
        os << fixed << setprecision(3) << setw(12) << right << commods[k].orig + 1 << "   " << setw(12) << right
           << commods[k].dest + 1 << "   " << setw(12) << right << *(commods[k].volume) << endl;
    }
}

// Reads a base deterministic MCFND network stored in restricted STD format.
//
// This is specialized to handle cases where there is a single source node
// and a single sink node for each commodity.
void InstanceRStd::read(istream &is)
{
    string line, dummy;
    int numbComm, idxComm, idxNode, capac, volume, varCost, fixedCost;
    int countArc = 0, countDemand = 0, countComm = 0;

    is >> numbNodes >> numbArcs >> numbCommods;

    allocate();

    for (int a = 0; a < numbArcs; a++)
    {
        Arc arc(numbCommods);

        is >> arc.orig >> arc.dest >> fixedCost // arc.cost
            >> capac >> numbComm;

        arc.orig--;
        arc.dest--;
        stochElems[CAP_ARC][countArc] = capac;
        stochElems[FIXED_COST][countArc] = fixedCost;
        nodes[arc.orig].arcOut.push_back(countArc);
        nodes[arc.dest].arcIn.push_back(countArc);
        countArc++;

        for (int k = 0; k < numbComm; k++)
        {
            is >> idxComm;
            idxComm--;
            is >> varCost // arc.comms[idxComm].cost
                >> capac;

            stochElems[CAP_COMMOD][countComm] = capac;
            stochElems[VAR_COST][countComm] = varCost;
            countComm++;
        }

        arcs[a] = arc;
    }

    for (int k = 0; k < numbCommods * 2; k++)
    {
        is >> idxComm >> idxNode >> volume;

        idxComm--;
        idxNode--;

        if (volume < 0)
            commods[idxComm].dest = idxNode;
        else
        {
            commods[idxComm].orig = idxNode;
            stochElems[DEMAND][countDemand++] = volume;
        }
    }
}

// Writes a generated MCFND networks to file in restricted STD format.
//
// This is specialized to handle cases where there is a single source
// node and a single sink node for each commodity.
void InstanceRStd::write(ostream &os) const
{
    static int scen = 0;
    os << " SCENARIO        " << scen++ << endl;
    os << setw(12) << right << numbNodes << "   " << setw(12) << right << numbArcs << "   " << setw(12) << right
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

                os << fixed << setprecision(3) << setw(12) << right << arcs[a].orig + 1 << "   " << setw(12) << right
                   << arcs[a].dest + 1 << "   " << setw(12) << right << *(arcs[a].cost) << "   " << setw(12) << right
                   << *(arcs[a].capac) << "   " << setw(12) << right << count << endl;

                for (int k = 0; k < numbCommods; k++)
                    if (arcs[a].comms[k].capac != 0)
                        os << fixed << setprecision(3) << setw(12) << right << k + 1 << "   " << setw(12) << right
                           << *(arcs[a].comms[k].cost) << "   " << setw(12) << right << *(arcs[a].comms[k].capac)
                           << endl;
            }
        }
    }

    for (int k = 0; k < numbCommods; k++)
    {
        os << fixed << setprecision(3) << setw(12) << right << k + 1 << "   " << setw(12) << right
           << commods[k].orig + 1 << "   " << setw(12) << right << *(commods[k].volume) << endl;
        os << fixed << setprecision(3) << setw(12) << right << k + 1 << "   " << setw(12) << right
           << commods[k].dest + 1 << "   " << setw(12) << right << -*(commods[k].volume) << endl;
    }
}

// Reads a base deterministic MCFND network stored in generic STD format.
void InstanceGStd::read(istream &is)
{
    string line, dummy;
    int numbComm, idxComm, idxNode, capac, volume, varCost, fixedCost;
    int countArc = 0, countComm = 0;

    is >> numbNodes >> numbArcs >> numbCommods;

    allocate();

    for (int a = 0; a < numbArcs; a++)
    {
        Arc arc(numbCommods);

        is >> arc.orig >> arc.dest >> fixedCost // arc.cost
            >> capac >> numbComm;

        arc.orig--;
        arc.dest--;
        stochElems[CAP_ARC][countArc] = capac;
        stochElems[FIXED_COST][countArc] = fixedCost;
        nodes[arc.orig].arcOut.push_back(countArc);
        nodes[arc.dest].arcIn.push_back(countArc);
        countArc++;

        for (int k = 0; k < numbComm; k++)
        {
            is >> idxComm;
            idxComm--;
            is >> varCost // arc.comms[idxComm].cost
                >> capac;

            stochElems[CAP_COMMOD][countComm] = capac;
            stochElems[VAR_COST][countComm] = varCost;
            countComm++;
        }

        arcs[a] = arc;
    }

    for (int k = 0; k < numbCommods; k++)
    {
        commods[k].outflows.resize(numbNodes, 0);
        commods[k].inflows.resize(numbNodes, 0);
        commods[k].origVolume = 0;
    }

    while (is >> idxComm >> idxNode >> volume)
    {
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
}

// Writes a generated MCFND networks to file in generic STD format.
void InstanceGStd::write(ostream &os) const
{
    static int scen = 0;
    os << " SCENARIO        " << scen++ << endl;
    os << setw(12) << right << numbNodes << "   " << setw(12) << right << numbArcs << "   " << setw(12) << right
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

                os << fixed << setprecision(3) << setw(12) << right << arcs[a].orig + 1 << "   " << setw(12) << right
                   << arcs[a].dest + 1 << "   " << setw(12) << right << *(arcs[a].cost) << "   " << setw(12) << right
                   << *(arcs[a].capac) << "   " << setw(12) << right << count << endl;

                for (int k = 0; k < numbCommods; k++)
                    if (arcs[a].comms[k].capac != 0)

                        os << fixed << setprecision(3) << setw(12) << right << k + 1 << "   " << setw(12) << right
                           << *(arcs[a].comms[k].cost) << "   " << setw(12) << right << *(arcs[a].comms[k].capac)
                           << endl;
            }
        }
    }

    for (int k = 0; k < numbCommods; k++)
    {
        for (int n = 0; n < numbNodes; n++)
        {

            // cout << "commod: " << k + 1 << " node: " << n + 1 << " orig outflow: " << commods[k].outflows[n]
            //      << " new volume: " << *(commods[k].volume) << " orig volume: " << commods[k].origVolume << endl;

            if (commods[k].outflows[n] > 0)
                os << fixed << setprecision(3) << setw(12) << right << k + 1 << "   " << setw(12) << right << n + 1
                   << "   " << setw(12) << right
                   << commods[k].outflows[n] * (*(commods[k].volume) / commods[k].origVolume) << endl;
            if (commods[k].inflows[n] < 0)
                os << fixed << setprecision(3) << setw(12) << right << k + 1 << "   " << setw(12) << right << n + 1
                   << "   " << setw(12) << right
                   << commods[k].inflows[n] * (*(commods[k].volume) / commods[k].origVolume) << endl;
        }
    }
}
