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
#ifndef __INSTANCE_H__
#define __INSTANCE_H__

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "configuration.h"
#include "generate.h"

using namespace std;

typedef unsigned int uint;
typedef vector<int> VectInt;

const int numb_moments = 4;

// Describes the characteristics of commodities as related to an arc.
struct CommodOnArc
{
    // pointer to static array holding variable costs for all arcs and
    // for all commodities over each arc
    static double *ptrCost;

    // pointer to static array holding capacities for all arcs and
    // for all commodities over each arc
    static double *ptrCap;

    // used to synchronize local values of cost and capac with static arrays
    static int count;

    // pointer to variable cost for this arc and this commodity
    double *cost;

    // pointer to capacity for this arc and this commodity
    double *capac;

    // unused
    friend ostream &operator<<(ostream &os, const CommodOnArc &c);

    CommodOnArc()
    {
        // connect capacity for this arc and this commodity with static array
        // holding commodity specific capacities for all arcs and commodities
        capac = &ptrCap[count];

        // connect variable cost for this arc and this commodity with static
        // array holding commodity specific variable costs for all arcs and
        // commodities
        cost = &ptrCost[count];
        count++;
    }

    // sets beginning of this static array if it is a subarray of static master
    // array (if CAP_COMMOD is among the stochastic elements)
    static void setPtrCapac(double *_ptr)
    {
        ptrCap = _ptr;
    }

    // sets beginning of this static array if it is a subarray of static master
    // array (if VAR_COST is among the stochastic elements)
    static void setPtrCost(double *_ptr)
    {
        ptrCost = _ptr;
    }
};

// will hold all commodities on an arc
typedef vector<CommodOnArc> VectCommArc;

// Defining the characteristics of a node.
struct Node
{
    // inbound arcs
    VectInt arcIn;

    // outbound arcs
    VectInt arcOut;
};

// will hold all nodes
typedef vector<Node> VectNode;

// Defining and handling the characteristics of the arcs.
struct Arc
{
    // pointer to static array holding fixed costs for all arcs
    static double *ptrCost;

    // pointer to static array holding capacities for all arcs
    static double *ptrCap;

    // used to synchronize local values of cost and capac with static arrays
    static int count;

    // origin of this arc
    int orig;

    // destination of this arc
    int dest;

    // pointer to fixed cost of this arc
    double *cost;

    // pointer to capacity of this arc
    double *capac;

    // vector of all CommodOnArcs associated with this arc
    VectCommArc comms;

    Arc() : orig(-1), dest(-1), cost(NULL), capac(NULL)
    {
    }

    Arc(int numbComm) : orig(-1), dest(-1)
    {
        // connect capacity for this arc with static array
        // holding capacities for all arcs
        capac = &ptrCap[count];
        // connect fixed cost for this arc with static array
        // holding fixed costs for all arcs
        cost = &ptrCost[count];
        comms.resize(numbComm);
        count++;
    }

    ~Arc()
    {
    }
    // sets beginning of this static array if it is a subarray of static
    // master array (if CAP_ARC is among the stochastic elements)
    static void setPtrCapac(double *_ptr)
    {
        ptrCap = _ptr;
    }

    // sets beginning of this static array if it is a subarray of static master array (if FIXED_COST is among the
    // stochastic elements)
    static void setPtrCost(double *_ptr)
    {
        ptrCost = _ptr;
    }
};

// will hold all arcs
typedef vector<Arc> VectArc;

// Defining and handling the characteristics of a commodity.
struct Commodity
{
    // pointer to static array holding volumes for all commodities
    static double *ptr;

    // used to synchronize local value of volume with static array
    static int count;

    // origin of this commodity
    int orig;

    // destination of this commodty
    int dest;

    // pointer to this volume
    double *volume;

    // (non-negative) outflows of this commodity from each node
    VectDbl outflows;

    // (non-positive) inflows of this commodity to each node
    VectDbl inflows;

    // original volume
    double origVolume;

    // unused
    // void connectToArray(int index)
    // {
    // // connect volume for this commodity with static array
    // // holding volumes for all commodities
    // volume = &ptr[index];
    // }

    Commodity() : orig(-1), dest(-1)
    {
        // connect volume for this commodity with static array
        // holding volumes for all commodities
        volume = &ptr[count++];
    }

    // sets beginning of this static array if it is a subarray of static
    // master array (if DEMAND is among the stochastic elements)
    static void setPtr(double *_ptr)
    {
        ptr = _ptr;
    }
};

// will hold all commodities
typedef vector<Commodity> VectCommod;

// Defining and handling the characteristics of a computing instance.
struct Instance
{
    // 1 Variables specifying characteristics of base deterministic instance.
    //
    // configuration parameters extracted from configuration file and
    // command line
    Configuration *params;

    // number of nodes
    int numbNodes;

    // number of arcs
    int numbArcs;

    // number of commodities
    int numbCommods;

    // array of Arc objects
    VectArc arcs;

    // array of Node objects
    VectNode nodes;

    // array of Commodity objects
    VectCommod commods;

    // 2 Variables characterizing scenarios deviating from base deterministic
    // instance.

    // Vector of StochElem objects. Each StochElem object is associated to a
    // particular type in STOCH_ELEMS_TYPE and each of these objects points
    // to an array of the stochastic elements that are meant to vary
    // stochastically between generated scenarios. Elements in this
    // array can be accessed directly using the [] operator defined in class
    // StochElem.
    VectElem stochElems;

    // stores scenarios generated by Hoyland-Kaut-Wallace algorithm
    TMatrix2 scenarios;

    // constructor
    Instance(Configuration *_params) : params(_params), numbNodes(-1), numbCommods(-1), scenarios(Mat_0)
    {
    }

    // destructor
    virtual ~Instance()
    {
        for (uint i = 0; i < params->nonStochTypes.size(); i++)
            free((void *)stochElems[params->nonStochTypes[i]].ptr);

        free((void *)stochElems[ALL_TYPES].ptr);

        Mat_Kill(&scenarios); // Mat_Kill() from scen-gen_HKW_src/matrix.h
    }

    friend istream &operator>>(istream &is, Instance &inst);

    friend ostream &operator<<(ostream &os, const Instance &inst);

    virtual void write([[maybe_unused]] ostream &os) const
    {
    }

    virtual void read([[maybe_unused]] istream &is)
    {
    }

    // Allocates arrays.
    void allocate();

    // Assigns sizes to subarrays comprising master array of stochastic
    // elements.
    void assignSizes();

    int getSize(STOCH_ELEMS_TYPE elem)
    {
        return stochElems[elem].size;
    }

    // Copies scenario to element array for purpose of verifying
    // feasibility with solver.
    void copyScenario(int idx);

    // Generates scenarios with Hoyland-Kaut-Wallace algorithm.
    int getScenarios();
};

// Reads instance from stream.
inline istream &operator>>(istream &is, Instance &inst)
{
    inst.read(is);
    return is;
}

// Writes instance to stream.
inline ostream &operator<<(ostream &os, const Instance &inst)
{
    inst.write(os);
    return os;
}

// Subclass handling DOW format
struct InstanceDow : public Instance
{
    InstanceDow(Configuration *params) : Instance(params)
    {
    }
    virtual void write(ostream &os) const override;
    virtual void read(istream &is) override;
};

// Subclass handling generic STD format
struct InstanceGStd : public Instance
{
    InstanceGStd(Configuration *params) : Instance(params)
    {
    }
    virtual void write(ostream &os) const override;
    virtual void read(istream &is) override;
};

// Subclass handling restricted STD format
struct InstanceRStd : public Instance
{
    InstanceRStd(Configuration *params) : Instance(params)
    {
    }
    virtual void write(ostream &os) const override;
    virtual void read(istream &is) override;
};

#endif
