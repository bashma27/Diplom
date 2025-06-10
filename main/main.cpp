#include "globals.h"
#include "solver_slae.h"
#include "gener_calc_area.h"
#include "flow_balancing.h"
#include "cell_state_recalc.h"
#include "test.h"

int main()
{
    setlocale(LC_ALL, "Russian");
    GenEndElGrid();
    GenArrayParallelepipeds();
    GenPortraitMatr();
    GenFirstBoundCondit();
    GenSecBoundCondit();
    GenArrayFictEndEl();
    BuildMatrVec();
    ConsiderBoundCondit();
    q.resize(NUM_NODES, 0);
    vector<double> r(NUM_NODES);
    vector<double> z(NUM_NODES);
    vector<double> Mult(NUM_NODES);
    vector<double> Az(NUM_NODES);
    int max_iter = 1000;
    double eps = 1e-15;
    MSG::LU_sq_MSG(q, r, z, Az, Mult, NUM_NODES, eps, max_iter);

    GenFictNodes();
    //RelativeErrorNorm();
    //EnergyNorm();
    //VecAnalitP();

    //BalancingFlows();
    //RecalcCellState();

    //OutputResult();
}


