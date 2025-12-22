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
#include "./model.h"

Model::Model(const Instance *_inst) : inst(_inst)
{
}

#if CPLEX

// Constructor
CPXModel::CPXModel(const Instance *_inst) : Model(_inst)
{
    iloModel = IloModel(iloEnv);

    cplex = IloCplex(iloModel);

    // set up 2-dimensional array of flow variables
    createXVars();

    // set up constraints maintaining flow conservation at each node
    createConstrFlowCons();

    // set up constraints restricting overall flow on each arc
    createConstrCapArcs();

    // silence CPLEX
    cplex.setOut(iloEnv.getNullStream());
    cplex.setWarning(iloEnv.getNullStream());
    cplex.setError(iloEnv.getNullStream());

    // cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0);

    // select the network algorithm
    cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Network);

    objExpr = IloExpr(iloEnv);

    for (int a = 0; a < inst->numbArcs; a++)
    {
        for (int k = 0; k < inst->numbCommods; k++)
        {
            objExpr += X_ak[a][k];
        }
    }

    IloObjective objective = IloMinimize(iloEnv, objExpr);
    iloModel.add(objective);
}

// Verifies if scenario is feasible.
bool CPXModel::testScenario()
{
    modify();

    cplex.solve();

    return (cplex.getStatus() == IloAlgorithm::Optimal);
}

// Modifies model to reflect scenario currently under test.
// Cases DOW or RSTD: one source node and one sink node per commodity.
void CPXModel_DowRStd::modify()
{
    int countComm = 0;
    // update source and sink volumes of each commodity
    for (int n = 0; n < inst->numbNodes; n++)
        for (int k = 0; k < inst->numbCommods; k++)
        {
            int bound = (int)inst->stochElems[DEMAND][k];

            if (inst->commods[k].dest == n)
                bound = -bound;
            else if (inst->commods[k].orig != n)
                bound = 0;

            constr_flowCons[n][k].setBounds(bound, bound);
        }

    // update overall capacity of each arc
    for (int a = 0; a < inst->numbArcs; a++)
        constr_capArcs[a].setUB(inst->stochElems[CAP_ARC][a]);

    // update commodity specific capacities of each arc
    for (int a = 0; a < inst->numbArcs; a++)
        for (int k = 0; k < inst->numbCommods; k++)
        {
            X_ak[a][k].setUB(inst->stochElems[CAP_COMMOD][countComm++]);
        }
}

// Modifies model to reflect scenario currently under test.
// Case GSTD: not one source node and one sink node per commodity.
void CPXModel_GStd::modify()
{
    int countComm = 0;
    // update source and sink volumes of each commodity
    if (std::find(inst->params->nonStochTypes.begin(), inst->params->nonStochTypes.end(), DEMAND) !=
        inst->params->nonStochTypes.end())
    {
        for (int n = 0; n < inst->numbNodes; n++)
            for (int k = 0; k < inst->numbCommods; k++)
            {
                // for given k and n, only one term on RHS is non-zero
                double bound = inst->commods[k].outflows[n] + inst->commods[k].inflows[n];

                constr_flowCons[n][k].setBounds(bound, bound);
            }
    }
    else
    {
        for (int n = 0; n < inst->numbNodes; n++)
            for (int k = 0; k < inst->numbCommods; k++)
            {
                // for given k and n, only one term on RHS is non-zero
                double bound =
                    inst->commods[k].outflows[n] * (inst->stochElems[DEMAND][k] / inst->commods[k].origVolume) +
                    inst->commods[k].inflows[n] * (inst->stochElems[DEMAND][k] / inst->commods[k].origVolume);

                constr_flowCons[n][k].setBounds(bound, bound);
            }
    }

    // update overall capacity of each arc
    for (int a = 0; a < inst->numbArcs; a++)
        constr_capArcs[a].setUB(inst->stochElems[CAP_ARC][a]);

    // update commodity specific capacities of each arc
    for (int a = 0; a < inst->numbArcs; a++)
        for (int k = 0; k < inst->numbCommods; k++)
        {
            X_ak[a][k].setUB(inst->stochElems[CAP_COMMOD][countComm++]);
        }
}

// Sets up 2-dimensional array of flow variables (first dim is arcs,
// second dim is commodities).
void CPXModel::createXVars()
{
    X_ak = IloNumVarArray2(iloEnv, inst->numbArcs);

    for (int a = 0; a < inst->numbArcs; a++)
    {
        X_ak[a] = IloNumVarArray(iloEnv, inst->numbCommods, 0, IloInfinity);
    }
}

// Sets up constraints restricting net flow on each node for each
// commodity.
void CPXModel::createConstrFlowCons()
{
    constr_flowCons = IloRangeArray2(iloEnv, inst->numbNodes);

    for (int n = 0; n < inst->numbNodes; n++)
    {
        constr_flowCons[n] = IloRangeArray(iloEnv);

        for (int k = 0; k < inst->numbCommods; k++)
        {
            IloExpr exp(iloEnv);

            for (uint a = 0; a < inst->nodes[n].arcOut.size(); a++)
                exp += X_ak[inst->nodes[n].arcOut[a]][k];

            for (uint a = 0; a < inst->nodes[n].arcIn.size(); a++)
                exp -= X_ak[inst->nodes[n].arcIn[a]][k];

            constr_flowCons[n].add(exp == 0.0);
            exp.end();
        }

        iloModel.add(constr_flowCons[n]);
    }
}

// Sets up constraints restricting total flow on each arc.
void CPXModel::createConstrCapArcs()
{
    constr_capArcs = IloRangeArray(iloEnv);

    for (int a = 0; a < inst->numbArcs; a++)
        constr_capArcs.add(IloSum(X_ak[a]) <= 0.0);

    iloModel.add(constr_capArcs);
}

#elif HIGHS

// Constructor
HIGHSModel::HIGHSModel(const Instance *_inst) : Model(_inst)
{
    highs.setOptionValue("output_flag", false);
    highs.setOptionValue("log_to_console", false);

    highs.setOptionValue("dual_feasibility_tolerance", 0.000001);

    numbRows1 = inst->numbNodes * inst->numbCommods;
    numbRows2 = inst->numbArcs;
    numbRows = numbRows1 + numbRows2;
    numbCols = inst->numbCommods * inst->numbArcs;

    prelAMatrix.resize(numbRows);

    for (int row = 0; row < numbRows; row++)

        prelAMatrix[row].resize(numbCols, 0);

    setUpperRows();
    setLowerRows();

    highsModel.lp_.num_col_ = numbCols;
    highsModel.lp_.num_row_ = numbRows;
    highsModel.lp_.sense_ = ObjSense::kMinimize;
    highsModel.lp_.offset_ = 0;
    highsModel.lp_.col_cost_.assign(numbCols, 1);

    highsModel.lp_.a_matrix_.format_ = MatrixFormat::kRowwise;

    highsModel.lp_.a_matrix_.index_ = {};
    highsModel.lp_.a_matrix_.value_ = {};
    highsModel.lp_.a_matrix_.start_ = {};

    setHIGHSAMatrix();
}

// Verifies if scenario is feasible.
bool HIGHSModel::testScenario()
{

    modify();

    highs.passModel(highsModel);

    highs.run();

    const HighsModelStatus &model_status = highs.getModelStatus();

    return (model_status == HighsModelStatus::kOptimal);
}

void HIGHSModel::setUpperRows()
{
    for (int n = 0; n < inst->numbNodes; n++)
    {
        for (int k = 0; k < inst->numbCommods; k++)
        {
            for (uint a0 = 0; a0 < inst->nodes[n].arcOut.size(); a0++)

                prelAMatrix[n * inst->numbCommods + k][inst->nodes[n].arcOut[a0] * inst->numbCommods + k] = 1;

            for (uint a1 = 0; a1 < inst->nodes[n].arcIn.size(); a1++)

                prelAMatrix[n * inst->numbCommods + k][inst->nodes[n].arcIn[a1] * inst->numbCommods + k] = -1;
        }
    }
}

void HIGHSModel::setLowerRows()
{
    for (int a = 0; a < inst->numbArcs; a++)
    {
        for (int k = 0; k < inst->numbCommods; k++)
        {
            prelAMatrix[a + numbRows1][a * inst->numbCommods + k] = 1;
        }
    }
}

void HIGHSModel::setHIGHSAMatrix()
{

    int countAll = 0;
    for (int row = 0; row < numbRows; row++)
    {
        bool first = true;
        for (int col = 0; col < numbCols; col++)
        {
            if (prelAMatrix[row][col] != 0)
            {
                highsModel.lp_.a_matrix_.index_.push_back(col);
                highsModel.lp_.a_matrix_.value_.push_back(prelAMatrix[row][col]);

                if (first)
                {
                    highsModel.lp_.a_matrix_.start_.push_back(countAll);
                    first = false;
                }
                countAll++;
            }
        }
    }
    highsModel.lp_.a_matrix_.start_.push_back(countAll);
}

void HIGHSModel_DowRStd::modify()
{

    double max_double_value = numeric_limits<double>::max();

    highsModel.lp_.row_lower_.assign(numbRows, 0);
    highsModel.lp_.row_upper_.assign(numbRows, 0);

    // update source and sink volumes of each commodity
    for (int n = 0; n < inst->numbNodes; n++)
    {
        for (int k = 0; k < inst->numbCommods; k++)
        {
            int bound = (int)inst->stochElems[DEMAND][k];

            if (inst->commods[k].dest == n)
                bound = -bound;
            else if (inst->commods[k].orig != n)
                bound = 0;

            highsModel.lp_.row_upper_[n * inst->numbCommods + k] = bound;
            highsModel.lp_.row_lower_[n * inst->numbCommods + k] = bound;
        }
    }

    // update overall capacity of each arc
    for (int a = 0; a < inst->numbArcs; a++)
        highsModel.lp_.row_upper_[numbRows1 + a] = inst->stochElems[CAP_ARC][a];

    highsModel.lp_.col_lower_.assign(numbCols, 0);
    highsModel.lp_.col_upper_.assign(numbCols, max_double_value);

    // update commodity specific capacities of each arc
    for (int a = 0; a < inst->numbArcs; a++)
    {
        for (int k = 0; k < inst->numbCommods; k++)
        {
            highsModel.lp_.col_upper_[a * inst->numbCommods + k] =
                inst->stochElems[CAP_COMMOD][a * inst->numbCommods + k];
        }
    }
}

void HIGHSModel_GStd::modify()
{
    double max_double_value = numeric_limits<double>::max();

    highsModel.lp_.row_lower_.assign(numbRows, 0);
    highsModel.lp_.row_upper_.assign(numbRows, 0);

    // update source and sink volumes of each commodity
    if (std::find(inst->params->nonStochTypes.begin(), inst->params->nonStochTypes.end(), DEMAND) !=
        inst->params->nonStochTypes.end())
    {
        for (int n = 0; n < inst->numbNodes; n++)
        {
            for (int k = 0; k < inst->numbCommods; k++)
            {
                // for given k and n, only one term on RHS is non-zero
                double bound = inst->commods[k].outflows[n] + inst->commods[k].inflows[n];

                highsModel.lp_.row_upper_[n * inst->numbCommods + k] = bound;
                highsModel.lp_.row_lower_[n * inst->numbCommods + k] = bound;
            }
        }
    }
    else
    {
        for (int n = 0; n < inst->numbNodes; n++)
        {
            for (int k = 0; k < inst->numbCommods; k++)
            {
                // for given k and n, only one term on RHS is non-zero
                double bound =
                    inst->commods[k].outflows[n] * (inst->stochElems[DEMAND][k] / inst->commods[k].origVolume) +
                    inst->commods[k].inflows[n] * (inst->stochElems[DEMAND][k] / inst->commods[k].origVolume);

                highsModel.lp_.row_upper_[n * inst->numbCommods + k] = bound;
                highsModel.lp_.row_lower_[n * inst->numbCommods + k] = bound;
            }
        }
    }

    // update overall capacity of each arc
    for (int a = 0; a < inst->numbArcs; a++)
        highsModel.lp_.row_upper_[numbRows1 + a] = inst->stochElems[CAP_ARC][a];

    highsModel.lp_.col_lower_.assign(numbCols, 0);
    highsModel.lp_.col_upper_.assign(numbCols, max_double_value);

    // update commodity specific capacities of each arc
    for (int a = 0; a < inst->numbArcs; a++)
    {
        for (int k = 0; k < inst->numbCommods; k++)
        {
            highsModel.lp_.col_upper_[a * inst->numbCommods + k] =
                inst->stochElems[CAP_COMMOD][a * inst->numbCommods + k];
        }
    }
}

#endif