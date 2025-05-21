#include "globals.h"
#include "solver_slae.h"
#include "gener_calc_area.h"
#include "flow_balancing.h"
#include "cell_state_recalc.h"

int main()
{
    setlocale(LC_ALL, "Russian");
    cout << "GenGrid" << endl;
    GenEndElGrid();
    cout << "GenArrayParallelepipeds" << endl;
    GenArrayParallelepipeds();
    cout << "GenPortraitMatr" << endl;
    GenPortraitMatr();
    cout << "GenFirstBoundCondit" << endl;
    GenFirstBoundCondit();
    cout << "GenSecBoundCondit" << endl;
    GenSecBoundCondit();
    cout << "BuildMatrVec" << endl;
    BuildMatrVec();
    cout << "ConsiderBoundCondit" << endl;
    ConsiderBoundCondit();
    q.resize(NUM_NODES, 0);
    vector<double> r(NUM_NODES);
    vector<double> z(NUM_NODES);
    vector<double> Mult(NUM_NODES);
    vector<double> Az(NUM_NODES);
    int max_iter = 1000;
    double eps = 1e-15;
    MSG::LU_sq_MSG(q, r, z, Az, Mult, NUM_NODES, eps, max_iter);
    GenArrayFictEndEl();
    BalancingFlows();
    RecalcCellState();

    /*GenFictNodes();
    Test();*/

    /*vector<double> vec_analit_P;
    VecAnalitP();*/

    //OutputResult();
}


