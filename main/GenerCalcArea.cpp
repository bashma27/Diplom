#include <fstream>
#include <vector>
#include <algorithm>
#define EPS 1e-16 //Для сравнения double

using namespace std;

extern struct Coord2 {
    double x, y;
    Coord2(double _x, double _y) {
        x = _x; y = _y;
    }
    Coord2() {
        x = 0; y = 0;
    }
};

extern struct Coord3 {
    double x, y, z;
    Coord3(double _x, double _y, double _z) {
        x = _x; y = _y; z = _z;
    }
    Coord3() {
        x = 0; y = 0; z = 0;
    }
};

extern vector<vector<Coord2>> zone_perf; // массив зон перфорации
extern vector<Coord3> nodes;
extern vector<double> set_flow_zp; // заданный поток зон перфорации
extern int num_ph;
extern int NUM_SPLIT_X, NUM_SPLIT_Y, NUM_SPLIT_Z; // суммарное количество разбиений по x, y, z
extern int NUM_ZONE_PERF; // количество зон перфорации
extern int NUM_NODES_IN_EDGE_X, NUM_NODES_IN_EDGE_Y, NUM_NODES_IN_EDGE_Z, NUM_NODES; // количество зон перфорации
extern vector<pair<int, vector<double>>> W; // массив подобластей (первое значение пары - номер подобласти, второе - массив значений границ подобласти по x, y, z)
extern vector<vector<double>> k_ph; //массив коэффициентов множителей структурной проницаемости,
extern vector<double> K, eta_ph; // массивы коэффициентов структурной проницаемости,    
//         коэффициентов динамической вязкости соответственно

int Nx, Ny, Nz, L; // Nx - число границ по x, Ny - число границ по y, Nz - число границ по z, L - кол-во подобластей
vector<double> Xw, Yw, Zw; // границы области
vector<vector<pair<int, double>>> coef; // Коэффициеты n, q  для каждой области в виде пары

int InputBorders() { // ввод координат границ подобластей
    ifstream File("grid_coord.txt");
    if (!File.is_open()) return 1;
    int n;
    double bound;
    File >> n;
    L = n - 1;
    for (int i = 0; i < n; i++) {
        File >> bound;
        Xw.push_back(bound);
    }
    File >> n;
    L *= n - 1;
    for (int i = 0; i < n; i++) {
        File >> bound;
        Yw.push_back(bound);
    }
    File >> n;
    L *= n - 1;
    for (int i = 0; i < n; i++) {
        File >> bound;
        Zw.push_back(bound);
    }
    W.resize(L);
    for (int i = 0; i < L; i++) {
        File >> W[i].first;
        W[i].first -= 1;
        W[i].second.resize(6);
        for (int j = 0; j < 6; j++) {
            File >> W[i].second[j];
        }
    }
    File.close();
    return 0;
}

bool InVectorInt(int value, vector<int>& temp) { // для записи в более коротком формате
    for (auto var : temp) {
        if (var == value)
            return false;
    }
    return true;
}

bool InVectorDouble(int value, vector<double>& temp) { // для записи в более коротком формате
    for (auto var : temp) {
        if (abs(var - value) < EPS)
            return true;
    }
    return false;
}

int InputZonePerf() { // ввод координат зон перфорации
    ifstream File("coord_zone_perf.txt");
    if (!File.is_open()) return 1;
    File >> NUM_ZONE_PERF;
    zone_perf.resize(NUM_ZONE_PERF);
    for (int i = 0; i < NUM_ZONE_PERF; i++) {
        zone_perf[i].resize(2);
        double value;
        File >> value; if (!InVectorDouble(value, Xw)) Xw.push_back(value);
        zone_perf[i][0].x = value;
        File >> value; if (!InVectorDouble(value, Xw)) Xw.push_back(value);
        zone_perf[i][1].x = value;
        File >> value; if (!InVectorDouble(value, Yw)) Yw.push_back(value);
        zone_perf[i][0].y = value;
        File >> value; if (!InVectorDouble(value, Yw)) Yw.push_back(value);
        zone_perf[i][1].y = value;
    }
    File.close();
    sort(Xw.begin(), Xw.end(), [](const double& a, const double& b) {return a < b; });
    sort(Yw.begin(), Yw.end(), [](const double& a, const double& b) {return a < b; });
    Nx = Xw.size(); Ny = Yw.size(); Nz = Zw.size();
    return 0;
}

int InputSetFlowZonePerf() { // ввод заданных потоков на зонах перфорации
    ifstream File("set_flow_zone_perf.txt");
    if (!File.is_open()) return 1;
    for (int i = 0; i < NUM_ZONE_PERF; i++) {
        set_flow_zp.resize(NUM_ZONE_PERF);
        File >> set_flow_zp[i];
    }
    return 0;
}

int CalcSum_ni(int i, int start, int stop) { // Посчитать сумму ni по порядку от start до end
    int value = 0;
    for (int k = start; k < stop; k++) {
        value += coef[i][k].first;
    }
    return value;
}

int Input_splits() {
    ifstream File("num_splits.txt");
    if (!File.is_open()) return 1;
    coef.resize(3);
    coef[0].resize(Nx - 1); coef[1].resize(Ny - 1); coef[2].resize(1);
    int ni;
    double qi;
    for (int i = 0; i < Nx - 1; i++) {
        File >> ni >> qi;
        coef[0][i] = { ni, qi };
    }
    for (int i = 0; i < Ny - 1; i++) {
        File >> ni >> qi;
        coef[1][i] = { ni, qi };
    }
    for (int i = 0; i < Nz - 1; i++) {
        File >> ni >> qi;
        coef[2][i] = { ni, qi };
    }
    int n1, n2, n3;
    File >> n1 >> n2 >> n3;
    for (int i = 0; i < Nx - 1; i++) {
        coef[0][i].first *= pow(2, n1);
    }
    for (int i = 0; i < Ny - 1; i++) {
        coef[1][i].first *= pow(2, n2);
    }
    for (int i = 0; i < Nz - 1; i++) {
        coef[2][i].first *= pow(2, n3);
    }

    for (int i = 0; i < Nx - 1; i++) {
        coef[0][i].second = pow(coef[0][i].second, 1. / pow(2, n1));
    }
    for (int i = 0; i < Ny - 1; i++) {
        coef[1][i].second = pow(coef[1][i].second, 1. / pow(2, n2));
    }
    for (int i = 0; i < Nz - 1; i++) {
        coef[2][i].second = pow(coef[2][i].second, 1. / pow(2, n3));
    }
    File.close();
    NUM_SPLIT_X = CalcSum_ni(0, 0, Nx - 1);
    NUM_SPLIT_Y = CalcSum_ni(1, 0, Ny - 1);
    NUM_SPLIT_Z = CalcSum_ni(2, 0, Nz - 1);
    NUM_NODES_IN_EDGE_X = 2 * CalcSum_ni(0, 0, Nx - 1) + 1;
    NUM_NODES_IN_EDGE_Y = 2 * CalcSum_ni(1, 0, Ny - 1) + 1;
    NUM_NODES_IN_EDGE_Z = 2 * CalcSum_ni(2, 0, Nz - 1) + 1;
    NUM_NODES = NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y * NUM_NODES_IN_EDGE_Z;
    nodes.resize(NUM_NODES);
    return 0;
}

int Input_coef() {
    ifstream File1("K.txt");
    if (!File1.is_open()) return 1;
    K.resize(L);
    for (int i = 0; i < L; i++) {
        File1 >> K[i];
        //K[i] *= 1e-15; // !!! перевод милиДарси в м2 (м2 - проницаемость в СИ)
    }
    File1.close();

    ifstream File2("ph.txt");
    if (!File2.is_open()) return 1;
    File2 >> num_ph;
    k_ph.assign(L, vector<double>(num_ph));
    eta_ph.resize(num_ph);
    for (int i = 0; i < num_ph; i++) {
        File2 >> eta_ph[i];
    }
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < num_ph; j++) {
            File2 >> k_ph[i][j];
        }
    }
    File2.close();
    return 0;
}
#pragma endregion

#pragma region Генерация сетки
void GenEndElGrid() {
    InputBorders();
    InputZonePerf();
    Input_splits();
    Input_coef();
    InputSetFlowZonePerf();
    int nx = 0, ny = 0, nz = 0, j = 0, i = 0; //nx и ny для подсчёта кол-ва разбиений до нужной точки
    int it = 0; // текущий номер последнего элемента в nodes
    double qx = 1.0, qy = 1.0, qz = 1.0, hy = 0.0, hx = 0.0, hz = 0.0, y = 0.0, z = 0.0;
    for (int _i = 0; _i < Nz - 1; _i++) {
        nz = coef[2][_i].first;
        qz = coef[2][_i].second;
        if (abs(qz - 1.0) < EPS)
            hz = (Zw[_i + 1] - Zw[_i]) / nz;
        else
            hz = (Zw[_i + 1] - Zw[_i]) * (qz - 1) / (pow(qz, nz) - 1);
        for (int m = 0; m < nz; m++) {
            if (abs(qz - 1.0) < EPS)
                z = Zw[_i] + hz * m;
            else
                z = Zw[_i] + hz * (pow(qz, m) - 1) / (qz - 1);
            for (i = 0; i < Ny - 1; i++) { // перебираем все коэфициенты по y
                ny = coef[1][i].first;
                qy = coef[1][i].second;
                if (abs(qy - 1.0) < EPS)
                    hy = (Yw[i + 1] - Yw[i]) / ny; // вычисления шага для равномерной сетки
                else
                    hy = (Yw[i + 1] - Yw[i]) * (qy - 1) / (pow(qy, ny) - 1); //вычисление h0 для неравномерной сетки
                for (int l = 0; l < ny; l++) { // цикл по каждому разбиению по y
                    if (abs(qy - 1.0) < EPS) // подсчёт текущего y
                        y = Yw[i] + hy * l;
                    else
                        y = Yw[i] + hy * (pow(qy, l) - 1) / (qy - 1);
                    for (j = 0; j < Nx - 1; j++) { // перебор всех коэфициентов по x
                        nx = coef[0][j].first;
                        qx = coef[0][j].second;
                        if (abs(qx - 1.0) < EPS)
                            hx = (Xw[j + 1] - Xw[j]) / nx; // вычисление шага для равномерной сетки
                        else
                            hx = (Xw[j + 1] - Xw[j]) * (qx - 1) / (pow(qx, nx) - 1); // вычисление h0 для неравномерной сетки
                        for (int k = 0; k < nx; k++) {
                            if (abs(qx - 1.0) < EPS) { // считаем текущий x и кладём всё в nodes
                                nodes[it].x = Xw[j] + k * hx;
                                nodes[it].y = y;
                                nodes[it].z = z;
                            }

                            else {
                                nodes[it].x = Xw[j] + hx * (pow(qx, k) - 1) / (qx - 1);
                                nodes[it].y = y;
                                nodes[it].z = z;
                            }

                            it += 2;
                        }
                    }
                    nodes[it].x = Xw[j];
                    nodes[it].y = y;  // кладётся конец для x, иначе бы элемент бы повторялся
                    nodes[it].z = z;

                    // добавление фиктивных узлов по x
                    for (int _i = it - NUM_NODES_IN_EDGE_X + 2; _i < it;) {
                        nodes[_i].x = (nodes[_i - 1].x + nodes[_i + 1].x) / 2.0;
                        nodes[_i].y = y;
                        nodes[_i].z = z;
                        _i += 2;
                    }
                    it += NUM_NODES_IN_EDGE_X + 1;
                }
            }
            //---------------------------------------------------------
            // дальше идёт тоже самое для последней строчки
            //---------------------------------------------------------
            y = Yw[i];
            for (j = 0; j < Nx - 1; j++) {
                nx = coef[0][j].first;
                qx = coef[0][j].second;
                if (abs(qx - 1.0) < EPS)
                    hx = (Xw[j + 1] - Xw[j]) / nx;
                else
                    hx = (Xw[j + 1] - Xw[j]) * (qx - 1) / (pow(qx, nx) - 1);
                for (int k = 0; k < nx; k++) {
                    if (abs(qx - 1.0) < EPS) {
                        nodes[it].x = Xw[j] + k * hx;
                        nodes[it].y = y;
                        nodes[it].z = z;
                    }
                    else {
                        nodes[it].x = Xw[j] + hx * (pow(qx, k) - 1) / (qx - 1);
                        nodes[it].y = y;
                        nodes[it].z = z;
                    }
                    it += 2;
                }
            }
            nodes[it].x = Xw[j];
            nodes[it].y = y;
            nodes[it].z = z;
            // добавление фиктивных узлов по x
            for (int _i = it - NUM_NODES_IN_EDGE_X + 2; _i < it;) {
                nodes[_i].x = (nodes[_i - 1].x + nodes[_i + 1].x) / 2.0;
                nodes[_i].y = y;
                nodes[_i].z = z;
                _i += 2;
            }
            it += NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 1;
        }
    }

    // добавляем в конец z
    z = Zw[Zw.size() - 1];
    for (i = 0; i < Ny - 1; i++) { // перебираем все коэфициенты по y
        ny = coef[1][i].first;
        qy = coef[1][i].second;
        if (abs(qy - 1.0) < EPS)
            hy = (Yw[i + 1] - Yw[i]) / ny; // вычисления шага для равномерной сетки
        else
            hy = (Yw[i + 1] - Yw[i]) * (qy - 1) / (pow(qy, ny) - 1); //вычисление h0 для неравномерной сетки
        for (int l = 0; l < ny; l++) { // цикл по каждому разбиению по y
            if (abs(qy - 1.0) < EPS) // подсчёт текущего y
                y = Yw[i] + hy * l;
            else
                y = Yw[i] + hy * (pow(qy, l) - 1) / (qy - 1);
            for (j = 0; j < Nx - 1; j++) { // перебор всех коэфициентов по x
                nx = coef[0][j].first;
                qx = coef[0][j].second;
                if (abs(qx - 1.0) < EPS)
                    hx = (Xw[j + 1] - Xw[j]) / nx; // вычисление шага для равномерной сетки
                else
                    hx = (Xw[j + 1] - Xw[j]) * (qx - 1) / (pow(qx, nx) - 1); // вычисление h0 для неравномерной сетки
                for (int k = 0; k < nx; k++) {
                    if (abs(qx - 1.0) < EPS) { // считаем текущий x и кладём всё в nodes
                        nodes[it].x = Xw[j] + k * hx;
                        nodes[it].y = y;
                        nodes[it].z = z;
                    }

                    else {
                        nodes[it].x = Xw[j] + hx * (pow(qx, k) - 1) / (qx - 1);
                        nodes[it].y = y;
                        nodes[it].z = z;
                    }

                    it += 2;
                }
            }
            nodes[it].x = Xw[j];
            nodes[it].y = y;  // кладётся конец для x, иначе бы элемент бы повторялся
            nodes[it].z = z;

            // добавление фиктивных узлов по x
            for (int _i = it - NUM_NODES_IN_EDGE_X + 2; _i < it;) {
                nodes[_i].x = (nodes[_i - 1].x + nodes[_i + 1].x) / 2.0;
                nodes[_i].y = y;
                nodes[_i].z = z;
                _i += 2;
            }
            it += NUM_NODES_IN_EDGE_X + 1;
        }
    }
    //---------------------------------------------------------
    // дальше идёт тоже самое для последней строчки
    //---------------------------------------------------------
    y = Yw[i];
    for (j = 0; j < Nx - 1; j++) {
        nx = coef[0][j].first;
        qx = coef[0][j].second;
        if (abs(qx - 1.0) < EPS)
            hx = (Xw[j + 1] - Xw[j]) / nx;
        else
            hx = (Xw[j + 1] - Xw[j]) * (qx - 1) / (pow(qx, nx) - 1);
        for (int k = 0; k < nx; k++) {
            if (abs(qx - 1.0) < EPS) {
                nodes[it].x = Xw[j] + k * hx;
                nodes[it].y = y;
                nodes[it].z = z;
            }
            else {
                nodes[it].x = Xw[j] + hx * (pow(qx, k) - 1) / (qx - 1);
                nodes[it].y = y;
                nodes[it].z = z;
            }
            it += 2;
        }
    }
    nodes[it].x = Xw[j];
    nodes[it].y = y;
    nodes[it].z = z;

    // добавление фиктивных узлов по x
    for (int _i = it - NUM_NODES_IN_EDGE_X + 2; _i < it;) {
        nodes[_i].x = (nodes[_i - 1].x + nodes[_i + 1].x) / 2.0;
        nodes[_i].y = y;
        nodes[_i].z = z;
        _i += 2;
    }

    // добавление фиктивных узлов по y
    for (int m = 0; m < nz + 1; m++) {
        for (int _i = NUM_NODES_IN_EDGE_X; _i < NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y - NUM_NODES_IN_EDGE_X;) {
            for (int __i = _i; __i < _i + NUM_NODES_IN_EDGE_X; __i++) {
                nodes[__i + 2 * m * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y].x = nodes[__i - NUM_NODES_IN_EDGE_X].x;
                nodes[__i + 2 * m * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y].y = (nodes[_i - NUM_NODES_IN_EDGE_X].y + nodes[_i + NUM_NODES_IN_EDGE_X].y) / 2.0;
                nodes[__i + 2 * m * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y].z = nodes[2 * m * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y].z;
            }
            _i += 2 * NUM_NODES_IN_EDGE_X;
        }
    }

    // добавление фиктивных узлов по z
    it = NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y;
    for (int m = 0; m < NUM_NODES_IN_EDGE_Z - nz - 1; m++) {
        for (int _i = 0; _i < NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y; _i++) {
            for (int __i = 0; __i < NUM_NODES_IN_EDGE_X; __i++) {
                nodes[it + _i + __i].x = nodes[it + _i + __i - NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y].x;
                nodes[it + _i + __i].y = nodes[it + _i + __i - NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y].y;
                nodes[it + _i + __i].z = (nodes[it - NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y].z + nodes[it + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y].z) / 2.0;
            }
            _i += NUM_NODES_IN_EDGE_X - 1;
        }
        it += 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y;
    }
}
#pragma endregion