#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <set>
#include <iomanip>
#include <math.h>
#include "SolversSLAE.h"
#include "GenerCalcArea.h"
#define EPS 1e-308

using namespace std;

#pragma region Объявление глобальных переменных
struct Coord2 {
    double x, y;
    Coord2(double _x, double _y) {
        x = _x; y = _y;
    }
    Coord2() {
        x = 0; y = 0;
    }
};

struct Coord3 {
    double x, y, z;
    Coord3(double _x, double _y, double _z) {
        x = _x; y = _y; z = _z;
    }
    Coord3() {
        x = 0; y = 0; z = 0;
    }
};

vector<vector<Coord2>> zone_perf; // массив зон перфорации
vector<pair<int, vector<double>>> W;
int NUM_ZONE_PERF; // количество зон перфорации
int NUM_SPLIT_X, NUM_SPLIT_Y, NUM_SPLIT_Z; // суммарное количество разбиений по x, y, z
int NUM_NODES_IN_EDGE_X, NUM_NODES_IN_EDGE_Y, NUM_NODES_IN_EDGE_Z, NUM_NODES;
int num_ph;
vector<double> set_flow_zp; // заданный поток зон перфорации
vector<int> ia, ja, choice;
vector<double> aal, di, b, q, L_sq, di_sq, normal;
vector<vector<double>> k_ph; //массив коэффициентов множителей структурной проницаемости,
vector<double> K, eta_ph; // массивы коэффициентов структурной проницаемости,    
//         коэффициентов динамической вязкости соответственно
vector<Coord3> nodes; // узлы сетки, заданные координатами 
unordered_set<int> face_1; // массив с узлами первых краевых
vector<pair<int, vector<int>>> face_2_zp; // массив граней с краевыми 2 рода зон перфорации
vector<pair<int, vector<int>>> array_p; // массив параллелепипедов
vector<vector<double>> G, M; // локальные матрицы масс и жесткости
vector<int> ident_fict; // идентификаторы фиктивности (области или конечного элемента)
vector<int> i_ident_fict;
vector<int> fict_nodes; // массив фиктивных узлов
vector<double> face_V; // массив граней со скорректированным значением потока
vector<pair<int, int>> num_faces_zp; // массив с номерами граней с краевыми 2 рода зон перфорации
vector<int> fict_faces; // номера фиктивных граней
vector<int> fict_el; // номера фиктивных конечных элементов
vector<double> q_V;
vector<vector<int>> list_face;

#pragma endregion

#pragma region Константы для аналитического решения
//double _theta = 180;
//double r_w = 0.5;
//double R = 105;
double P_g = 130; // !!! [атмосферы] * (коэф для перевода атм в Па) 
//double _theta = 180;
//double _theta = 0.79 / 86400.0 / (2. * 3.1415926 * r_w * 1.); // !!! [м3/сут] / (секунды в 24 часах) / (2*PI*R_well*h_well)
#pragma endregion

#pragma region Функции краевых условий, их функций и вектора правой части f
double u_g(double x, double y, double z) { // краевое условие первого рода
    //return x * x * y * y;
    //return x * x + y * y;
    return P_g;
    //return y * y * y;
    //return x + y;
    //return y * y * y * y;
}

//vector<double> grad_u(double x, double y) { // градиент функции u
//    //return { 2 * x * y * y, 2 * y * x * x };
//    //return { 2 * x, 2 * y };
//    //return { 0, 3 * y * y };
//    //return { 1, 1 };
//    //return { 0, 4 * y * y * y};
//}

double lambda(int num_sub) {
    double res = 0;
    for (int i = 0; i < num_ph; i++) {
        res += k_ph[num_sub][i] / eta_ph[i];
    }
    res *= K[num_sub];
    return res;
}

//double theta(int num_sub, int num, double x, double y) { // краевое условие второго рода
//    switch (num)
//    {
//    /*case 0:
//        return grad_u(x, y)[1] * 1 * lambda(num_sub);
//        break;
//    case 1:
//        return grad_u(x, y)[1] * (-1) * lambda(num_sub);
//        break;
//    case 2:
//        return grad_u(x, y)[0] * 1 * lambda(num_sub);
//        break;
//    case 3:
//        return grad_u(x, y)[0] * (-1) * lambda(num_sub);
//        break;*/
//
//    case 0:
//        return _theta;
//        break;
//    case 1:
//        return _theta;
//        break;
//    case 2:
//        return _theta;
//        break;
//    case 3:
//        return _theta;
//        break;
//    }
//}

vector<double> theta_V(int num_face_2_zp) { // 0 - тетта, 1 - поток (для зоны перфорации)
    int num_sub = face_2_zp[num_face_2_zp].second[9];
    int num_zone_perf = face_2_zp[num_face_2_zp].second[10];
    int num = face_2_zp[num_face_2_zp].first;

    double h_x, h_y, h_z; // длинна, ширина и высота зоны перфорации
    h_x = zone_perf[num_zone_perf][1].x - zone_perf[num_zone_perf][0].x;
    h_y = zone_perf[num_zone_perf][1].y - zone_perf[num_zone_perf][0].y;
    h_z = W[num_sub].second[5] - W[num_sub].second[4]; // высота по z (так как все зоны перфорации проходят через всю область )
    double S_area_zp = 2. * h_x * h_z + 2. * h_y * h_z; // площадь поверхности зоны перфорации
    double S_curr_face;

    double h_1, h_2;
    if (num == 0 || num == 1) {
        h_1 = nodes[face_2_zp[num_face_2_zp].second[2]].x - nodes[face_2_zp[num_face_2_zp].second[0]].x;
        h_2 = nodes[face_2_zp[num_face_2_zp].second[6]].z - nodes[face_2_zp[num_face_2_zp].second[0]].z;
    }
    else {
        h_1 = nodes[face_2_zp[num_face_2_zp].second[2]].y - nodes[face_2_zp[num_face_2_zp].second[0]].y;
        h_2 = nodes[face_2_zp[num_face_2_zp].second[6]].z - nodes[face_2_zp[num_face_2_zp].second[0]].z;
    }

    S_curr_face = h_1 * h_2;
    return { set_flow_zp[num_zone_perf] * S_curr_face / (S_area_zp * S_area_zp), set_flow_zp[num_zone_perf] * S_curr_face / S_area_zp };
}

double f(int num_sub, double x, double y) {
    //return - 2 * y * y * lambda(num_sub) - 2 * x * x * lambda(num_sub);
    //return -4 * lambda(num_sub);
    //return -6 * y * lambda(num_sub);
    return 0;
    //return - 12 * y * y * lambda(num_sub);
}
#pragma endregion

#pragma region Создание массива конечных элементов 
int GetNumSubarea(vector<int> end_el) { // получить номер подобласти
    for (int i = 0; i < W.size(); i++) {
        double min_x = nodes[end_el[0]].x; double min_y = nodes[end_el[0]].y; double min_z = nodes[end_el[0]].z;  //берём индексы координат прямоугольника
        double max_x = nodes[end_el[2]].x; double max_y = nodes[end_el[6]].y; double max_z = nodes[end_el[18]].z;
        if ((min_x > W[i].second[0] || fabs(W[i].second[0] - min_x) < EPS) &&
            (W[i].second[1] > max_x || fabs(max_x - W[i].second[1]) < EPS) &&
            (min_y > W[i].second[2] || fabs(W[i].second[2] - min_y) < EPS) &&
            (W[i].second[3] > max_y || fabs(max_y - W[i].second[3]) < EPS) &&
            (min_z > W[i].second[4] || fabs(W[i].second[4] - min_z) < EPS) &&
            (W[i].second[5] > max_z || fabs(max_z - W[i].second[5]) < EPS))
            return W[i].first;
    }
}

void ArrayParallelepipeds() { // массив параллелепипедов
    array_p.resize(NUM_SPLIT_X * NUM_SPLIT_Y * NUM_SPLIT_Z);
    int j, k;
    k = 0;
    for (int i = 0; i < array_p.size();) {
        for (int m = 0, j = 0; j < NUM_SPLIT_X * NUM_SPLIT_Y; m++) {
            for (int l = 0; l < NUM_SPLIT_X; l++) {
                int l_d_node = k + m * 2 * NUM_NODES_IN_EDGE_X + 2 * l; // нижний левый номер узла конечного элемента    
                array_p[i + j + l].second.push_back(l_d_node);
                array_p[i + j + l].second.push_back(l_d_node + 1);
                array_p[i + j + l].second.push_back(l_d_node + 2);
                array_p[i + j + l].second.push_back(l_d_node + NUM_NODES_IN_EDGE_X);
                array_p[i + j + l].second.push_back(l_d_node + 1 + NUM_NODES_IN_EDGE_X);
                array_p[i + j + l].second.push_back(l_d_node + 2 + NUM_NODES_IN_EDGE_X);
                array_p[i + j + l].second.push_back(l_d_node + 2 * NUM_NODES_IN_EDGE_X);
                array_p[i + j + l].second.push_back(l_d_node + 1 + 2 * NUM_NODES_IN_EDGE_X);
                array_p[i + j + l].second.push_back(l_d_node + 2 + 2 * NUM_NODES_IN_EDGE_X);

                array_p[i + j + l].second.push_back(l_d_node + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
                array_p[i + j + l].second.push_back(l_d_node + 1 + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
                array_p[i + j + l].second.push_back(l_d_node + 2 + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
                array_p[i + j + l].second.push_back(l_d_node + NUM_NODES_IN_EDGE_X + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
                array_p[i + j + l].second.push_back(l_d_node + 1 + NUM_NODES_IN_EDGE_X + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
                array_p[i + j + l].second.push_back(l_d_node + 2 + NUM_NODES_IN_EDGE_X + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
                array_p[i + j + l].second.push_back(l_d_node + 2 * NUM_NODES_IN_EDGE_X + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
                array_p[i + j + l].second.push_back(l_d_node + 1 + 2 * NUM_NODES_IN_EDGE_X + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
                array_p[i + j + l].second.push_back(l_d_node + 2 + 2 * NUM_NODES_IN_EDGE_X + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);

                array_p[i + j + l].second.push_back(l_d_node + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
                array_p[i + j + l].second.push_back(l_d_node + 1 + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
                array_p[i + j + l].second.push_back(l_d_node + 2 + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
                array_p[i + j + l].second.push_back(l_d_node + NUM_NODES_IN_EDGE_X + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
                array_p[i + j + l].second.push_back(l_d_node + 1 + NUM_NODES_IN_EDGE_X + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
                array_p[i + j + l].second.push_back(l_d_node + 2 + NUM_NODES_IN_EDGE_X + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
                array_p[i + j + l].second.push_back(l_d_node + 2 * NUM_NODES_IN_EDGE_X + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
                array_p[i + j + l].second.push_back(l_d_node + 1 + 2 * NUM_NODES_IN_EDGE_X + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
                array_p[i + j + l].second.push_back(l_d_node + 2 + 2 * NUM_NODES_IN_EDGE_X + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
                array_p[i + j + l].first = GetNumSubarea(array_p[i + j + l].second);
            }
            j += NUM_SPLIT_X;
        }
        j = 0;
        k += 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y;
        i += NUM_SPLIT_X * NUM_SPLIT_Y;
    }
    /*for (int i = 0; i < array_p.size(); i++) {
        for (int j = 0; j < 27; j++) {
            cout << nodes[array_p[i][j]].x << " " << nodes[array_p[i][j]].y << " " << nodes[array_p[i][j]].z << endl;
        }
        cout << endl;
    }*/
}
#pragma endregion

#pragma region Генерация портрета и подпр-мы добавления в глоб. матрицу
void GeneratePortrait() { // генерация портрета матрицы
    vector<set<int>> list(NUM_NODES);
    for (int i = 0; i < NUM_SPLIT_X * NUM_SPLIT_Y * NUM_SPLIT_Z; i++) {
        for (int j = 26; j >= 0; j--) {
            for (int m = j - 1; m >= 0; m--) {
                list[array_p[i].second[j]].insert(array_p[i].second[m]);
            }
        }
    }
    ia.resize(NUM_NODES + 1, 0);
    for (int i = 2; i <= NUM_NODES; i++) {
        ia[i] = ia[i - 1] + list[i - 1].size();
    }
    ja.resize(ia[NUM_NODES]);
    auto iter = ja.begin();
    for (int i = 0; i < NUM_NODES; i++) {
        copy(list[i].begin(), list[i].end(), iter);
        iter += list[i].size();
    }
    aal.resize(ia[NUM_NODES]);
    di.resize(NUM_NODES);
    b.resize(NUM_NODES);
    L_sq.resize(ia[NUM_NODES]);
    di_sq.resize(NUM_NODES);
}

// !!! все сложные структуры стоит передавать по указател/ссылке, чтобы не создавалась копия; если надо быть уверенным, что ничего в этом объекте не изменится, можно использовать const
void AddLocalMatr(vector<int>& node_num, vector<vector<double>>& loc_matr) { // добавление матрицы в глобальную
    for (int i = 0; i < 27; i++) {
        di[node_num[i]] += loc_matr[i][i];
    }
    for (int i = 0; i < 27; i++) {
        int i_beg = ia[node_num[i]];
        for (int j = 0; j < i; j++) {
            int i_end = ia[node_num[i] + 1];
            while (ja[i_beg] != node_num[j]) {
                int ind = (i_beg + i_end) / 2;
                if (ja[ind] <= node_num[j]) i_beg = ind;
                else i_end = ind;
            }
            aal[i_beg] += loc_matr[i][j];
            i_beg++;
        }
    }
}

void AddLocalVec(vector<int> node_num, vector<double>& loc_vec) { // добавление вектора в глобальный
    for (int i = 0; i < 27; i++) {
        b[node_num[i]] += loc_vec[i];
    }
}

void AddLocalVecBound(int num_face_2_zp, vector<double>& vec) { // добавление вектора из второго краевого в глобальный
    for (int i = 0; i < 9; i++) {
        b[face_2_zp[num_face_2_zp].second[i]] += vec[i];
    }
}
#pragma endregion

#pragma region Работа с краевыми условиями
void GenFirstBoundCondit() { // создание массива с первыми краевыми условиями
    int start_nodes = 0;
    int offset = 0;
    for (int i = 0; i < NUM_SPLIT_X * NUM_SPLIT_Z; i++) {
        face_1.insert(offset + start_nodes);
        face_1.insert(offset + start_nodes + 1);
        face_1.insert(offset + start_nodes + 2);

        face_1.insert(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
        face_1.insert(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 1);
        face_1.insert(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 2);

        face_1.insert(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
        face_1.insert(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 1);
        face_1.insert(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 2);
        start_nodes += 2;
        if ((i + 1) % NUM_SPLIT_X == 0) {
            start_nodes = 0;
            offset += 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y;
        }
    }
    start_nodes = NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y - NUM_NODES_IN_EDGE_X;
    offset = 0;
    for (int i = NUM_SPLIT_X * NUM_SPLIT_Z; i < 2 * NUM_SPLIT_X * NUM_SPLIT_Z; i++) {
        face_1.insert(offset + start_nodes);
        face_1.insert(offset + start_nodes + 1);
        face_1.insert(offset + start_nodes + 2);

        face_1.insert(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
        face_1.insert(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 1);
        face_1.insert(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 2);

        face_1.insert(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
        face_1.insert(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 1);
        face_1.insert(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 2);
        start_nodes += 2;
        if ((i + 1) % NUM_SPLIT_X == 0) {
            start_nodes = NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y - NUM_NODES_IN_EDGE_X;
            offset += 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y;
        }
    }
    start_nodes = 0;
    offset = 0;
    for (int i = 2 * NUM_SPLIT_X * NUM_SPLIT_Z; i < 2 * NUM_SPLIT_X * NUM_SPLIT_Z + NUM_SPLIT_Y * NUM_SPLIT_Z; i++) {
        face_1.insert(offset + start_nodes);
        face_1.insert(offset + start_nodes + NUM_NODES_IN_EDGE_X);
        face_1.insert(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X);

        face_1.insert(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
        face_1.insert(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + NUM_NODES_IN_EDGE_X);
        face_1.insert(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 2 * NUM_NODES_IN_EDGE_X);

        face_1.insert(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
        face_1.insert(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + NUM_NODES_IN_EDGE_X);
        face_1.insert(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 2 * NUM_NODES_IN_EDGE_X);
        start_nodes += 2 * NUM_NODES_IN_EDGE_X;
        if ((i + 1) % NUM_SPLIT_Y == 0) {
            start_nodes = 0;
            offset += 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y;
        }
    }
    start_nodes = NUM_NODES_IN_EDGE_X - 1;
    offset = 0;
    for (int i = 2 * NUM_SPLIT_X * NUM_SPLIT_Z + NUM_SPLIT_Y * NUM_SPLIT_Z; i < 2 * NUM_SPLIT_X * NUM_SPLIT_Z + 2 * NUM_SPLIT_Y * NUM_SPLIT_Z; i++) {
        face_1.insert(offset + start_nodes);
        face_1.insert(offset + start_nodes + NUM_NODES_IN_EDGE_X);
        face_1.insert(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X);

        face_1.insert(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
        face_1.insert(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + NUM_NODES_IN_EDGE_X);
        face_1.insert(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 2 * NUM_NODES_IN_EDGE_X);

        face_1.insert(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
        face_1.insert(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + NUM_NODES_IN_EDGE_X);
        face_1.insert(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 2 * NUM_NODES_IN_EDGE_X);
        start_nodes += 2 * NUM_NODES_IN_EDGE_X;
        if ((i + 1) % NUM_SPLIT_Y == 0) {
            start_nodes = NUM_NODES_IN_EDGE_X - 1;
            offset += 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y;
        }
    }
}

int GetNumNodes(vector<Coord2> zone_perf_i, int i_x, int i_y) { // получить номер узла для конкретной зоны перфорации
    for (int i = 0; i < NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y - 1; i++) {
        if ((zone_perf_i[i_x].x == nodes[i].x || abs(zone_perf_i[i_x].x - nodes[i].x) < EPS) &&
            (zone_perf_i[i_y].y == nodes[i].y || abs(zone_perf_i[i_y].y - nodes[i].y) < EPS)) return i;
    }
}

int GetNumSubareaBound(int ind, vector<int> face_2_zp) { // получить номер подобласти для второго краевого (первый индекс 0/1 - грань параллельная оси y/x)
    for (int i = 0; i < W.size(); i++) {
        if (ind == 0) {
            double x1 = nodes[face_2_zp[0]].x; double z1 = nodes[face_2_zp[0]].z;  //берём индексы координат прямоугольника
            double x2 = nodes[face_2_zp[2]].x; double z2 = nodes[face_2_zp[6]].z;
            if ((EPS < x1 - W[i].second[0]) &&
                (W[i].second[1] - x2 > EPS) &&
                (W[i].second[4] == z1 || EPS < z1 - W[i].second[4]) &&
                (W[i].second[5] == z2 || W[i].second[5] - z2) > EPS) return W[i].first;
        }
        else {
            double y1 = nodes[face_2_zp[0]].y; double z1 = nodes[face_2_zp[0]].z;  //берём индексы координат прямоугольника
            double y2 = nodes[face_2_zp[2]].y; double z2 = nodes[face_2_zp[6]].z;
            if ((EPS < y1 - W[i].second[2]) &&
                (W[i].second[3] - y2 > EPS) &&
                (W[i].second[4] == z1 || EPS < z1 - W[i].second[4]) &&
                (W[i].second[5] == z2 || W[i].second[5] - z2) > EPS) return W[i].first;
        }
    }
}

void GenSecBoundCondit() { // создание массива с вторыми краевыми условиями
    int size_face_2 = 0;
    int num_mis_finite_el = 0; // количество отсутствующих конечных элементов
    for (int ii = 0; ii < NUM_ZONE_PERF; ii++) {
        int start_nodes_xy, end_nodes_x, end_nodes_y, num_split_x, num_split_y;
        start_nodes_xy = GetNumNodes(zone_perf[ii], 0, 0);
        end_nodes_x = GetNumNodes(zone_perf[ii], 1, 0); end_nodes_y = GetNumNodes(zone_perf[ii], 0, 1);
        num_split_x = (end_nodes_x - start_nodes_xy) / 2; num_split_y = (end_nodes_y - start_nodes_xy) / (2 * NUM_NODES_IN_EDGE_X);
        num_mis_finite_el += num_split_x * num_split_y;
        size_face_2 += 2 * num_split_x * NUM_SPLIT_Z + 2 * num_split_y * NUM_SPLIT_Z;
        face_2_zp.resize(size_face_2);

        int prev_size_fase_2 = size_face_2 - 2 * num_split_x * NUM_SPLIT_Z - 2 * num_split_y * NUM_SPLIT_Z;
        int start_nodes = start_nodes_xy;
        int offset = 0;
        for (int i = 0; i < num_split_x * NUM_SPLIT_Z; i++) {
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + 1);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + 2);

            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 1);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 2);

            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 1);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 2);
            int num_sub_area = GetNumSubareaBound(0, face_2_zp[i + prev_size_fase_2].second);
            face_2_zp[i + prev_size_fase_2].second.push_back(num_sub_area);
            face_2_zp[i + prev_size_fase_2].second.push_back(ii);
            face_2_zp[i + prev_size_fase_2].first = 0;
            start_nodes += 2;
            if ((i + 1) % num_split_x == 0) {
                start_nodes = start_nodes_xy;
                offset += 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y;
            }
        }
        start_nodes = end_nodes_y;
        offset = 0;
        for (int i = num_split_x * NUM_SPLIT_Z; i < 2 * num_split_x * NUM_SPLIT_Z; i++) {
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + 1);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + 2);

            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 1);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 2);

            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 1);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 2);
            int num_sub_area = GetNumSubareaBound(0, face_2_zp[i + prev_size_fase_2].second);
            face_2_zp[i + prev_size_fase_2].second.push_back(num_sub_area);
            face_2_zp[i + prev_size_fase_2].second.push_back(ii);
            face_2_zp[i + prev_size_fase_2].first = 1;
            start_nodes += 2;
            if ((i + 1) % num_split_x == 0) {
                start_nodes = end_nodes_y;
                offset += 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y;
            }
        }
        start_nodes = start_nodes_xy;
        offset = 0;
        for (int i = 2 * num_split_x * NUM_SPLIT_Z; i < 2 * num_split_x * NUM_SPLIT_Z + num_split_y * NUM_SPLIT_Z; i++) {
            ident_fict.push_back(offset + start_nodes);
            i_ident_fict.push_back(num_split_x);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + NUM_NODES_IN_EDGE_X);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X);

            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + NUM_NODES_IN_EDGE_X);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 2 * NUM_NODES_IN_EDGE_X);

            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + NUM_NODES_IN_EDGE_X);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 2 * NUM_NODES_IN_EDGE_X);
            int num_sub_area = GetNumSubareaBound(1, face_2_zp[i + prev_size_fase_2].second);
            face_2_zp[i + prev_size_fase_2].second.push_back(num_sub_area);
            face_2_zp[i + prev_size_fase_2].second.push_back(ii);
            face_2_zp[i + prev_size_fase_2].first = 2;
            start_nodes += 2 * NUM_NODES_IN_EDGE_X;
            if ((i + 1) % num_split_y == 0) {
                start_nodes = start_nodes_xy;
                offset += 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y;
            }
        }
        start_nodes = end_nodes_x;
        offset = 0;
        for (int i = 2 * num_split_x * NUM_SPLIT_Z + num_split_y * NUM_SPLIT_Z; i < 2 * num_split_x * NUM_SPLIT_Z + 2 * num_split_y * NUM_SPLIT_Z; i++) {
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + NUM_NODES_IN_EDGE_X);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X);

            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + NUM_NODES_IN_EDGE_X);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 2 * NUM_NODES_IN_EDGE_X);

            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + NUM_NODES_IN_EDGE_X);
            face_2_zp[i + prev_size_fase_2].second.push_back(offset + start_nodes + 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y + 2 * NUM_NODES_IN_EDGE_X);
            int num_sub_area = GetNumSubareaBound(1, face_2_zp[i + prev_size_fase_2].second);
            face_2_zp[i + prev_size_fase_2].second.push_back(num_sub_area);
            face_2_zp[i + prev_size_fase_2].second.push_back(ii);
            face_2_zp[i + prev_size_fase_2].first = 3;
            start_nodes += 2 * NUM_NODES_IN_EDGE_X;
            if ((i + 1) % num_split_y == 0) {
                start_nodes = end_nodes_x;
                offset += 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y;
            }
        }
    }
}

void BuildLocalMatrices() { // построение локальных матриц массы и жесткости
    G.assign(3, vector<double>(3));
    M.assign(3, vector<double>(3));
    G[0][0] = 7, M[0][0] = 4;
    G[0][1] = G[1][0] = -8, M[0][1] = M[1][0] = 2;
    G[0][2] = G[2][0] = 1, M[0][2] = M[2][0] = -1;
    G[1][1] = 16, M[1][1] = 16;
    G[1][2] = G[2][1] = -8, M[1][2] = M[2][1] = 2;
    G[2][2] = 7, M[2][2] = 4;
}

void ConsiderBoundConditFirstType(int n) { // учет краевых условий первого типа
    double x, y, z;
    x = nodes[n].x, y = nodes[n].y, z = nodes[n].z;
    b[n] = u_g(x, y, z);
    di[n] = 1;
    for (int i = ia[n]; i < ia[n + 1]; i++) {
        int _i = ja[i];
        if (face_1.find(_i) != face_1.end()) {
            aal[i] = 0;
            continue;
        }
        b[_i] -= b[n] * aal[i];
        aal[i] = 0;
    }
    for (int i = n; i < NUM_NODES; i++) {
        int k = 0;
        for (int j = ia[i]; j < ia[i + 1]; j++) {
            if (ja[j] == n) {
                if (face_1.find(i) != face_1.end()) {
                    aal[j] = 0;
                    continue;
                }
                b[i] -= b[n] * aal[j];
                aal[j] = 0;
            }
        }
    }
}

bool IsFictitious(int num) { // определение фиктивности узла
    double x = nodes[num].x, y = nodes[num].y;
    int x1, x2, y1, y2;
    for (int k = 0; k < NUM_ZONE_PERF; k++) {
        x1 = zone_perf[k][0].x; y1 = zone_perf[k][0].y; //берём индексы координат прямоугольника
        x2 = zone_perf[k][1].x; y2 = zone_perf[k][1].y;
        if ((x1 - x) < EPS && x1 != x && EPS < (x2 - x) && x2 != x && (y1 - y) < EPS && y1 != y && EPS < (y2 - y) && y2 != y) // проверка на то что узел находится в прямоугольнике
            return true;
    }
    return false;
}

void CreateVecFict() {
    for (int i = 0; i < NUM_NODES; i++) {
        if (IsFictitious(i)) {
            fict_nodes.push_back(i);
        }
    }
}

void ConsiderFictitiousNodes() { // учет фиктивных узлов
    // !!! вот эта проверка может быть долгой, особенно с дроблением сетки; 
    // как вариант - можно просто проверять, что диагональный элемент 0 (если была инициальзация нулями и не было вкладов - то элемент диагонали останется нулем)
    // или собрать список фиктивных узлов, как и элементов - их будет мало относительно полной сетки
    for (int i = 0; i < NUM_NODES; i++) {
        //if (IsFictitious(i)) {
        if (di[i] == 0) {
            b[i] = 0;
            di[i] = double(1);
        }
    }
}

void ConsiderBoundConditSecType(int num_face_2_zp) { // учет краевых условий второго типа
    vector<double> b_s2(9, 0);
    double _theta = theta_V(num_face_2_zp)[0];
    //for (int i = 0; i < 9; i++) {
    //    /*_theta[i] = theta(face_2_zp[num_face_2_zp].second[9], face_2_zp[num_face_2_zp].first,
    //        nodes[face_2_zp[num_face_2_zp].second[i]].x,
    //        nodes[face_2_zp[num_face_2_zp].second[i]].y);*/
    //}
    double h_1, h_2;
    if (face_2_zp[num_face_2_zp].first == 0 || face_2_zp[num_face_2_zp].first == 1) {
        h_1 = nodes[face_2_zp[num_face_2_zp].second[2]].x - nodes[face_2_zp[num_face_2_zp].second[0]].x;
        h_2 = nodes[face_2_zp[num_face_2_zp].second[6]].z - nodes[face_2_zp[num_face_2_zp].second[0]].z;
    }
    else {
        h_1 = nodes[face_2_zp[num_face_2_zp].second[2]].y - nodes[face_2_zp[num_face_2_zp].second[0]].y;
        h_2 = nodes[face_2_zp[num_face_2_zp].second[6]].z - nodes[face_2_zp[num_face_2_zp].second[0]].z;
    }

    vector<vector<double>> C;
    C.assign(9, vector<double>(9));

    for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 9; j++) {
            C[i][j] = h_1 * h_2 * M[i % 3][j % 3] * M[i / 3][j / 3] / 900.0;
        }
    }
    for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 9; j++) {
            //b_s2[i] += C[i][j] * _theta[j];
            b_s2[i] += C[i][j] * _theta;
        }
    }
    AddLocalVecBound(num_face_2_zp, b_s2);
}

void ConsiderBoundCondit() { // учет всех краевых
    for (int i = 0; i < face_2_zp.size(); i++) { // учет вторых краевых
        ConsiderBoundConditSecType(i);
    }
    cout << "second" << endl;
    for (int n : face_1) { // учет первых краевых
        ConsiderBoundConditFirstType(n);
    }
    // !!! первые краевые подозрительно долго применяются; учитывая поиски в цикле возможно стоит или делать обычный set (чтобы он сам сортировался), или что-то с хэш-поиском использовать (вроде map) 
    cout << "first" << endl;
    ConsiderFictitiousNodes();
    cout << "fiction" << endl;

}
#pragma endregion

#pragma region Работа с локальными матрицами и векторами
void LocalMatrStiff(int num_sub, vector<int> node_num) { // локаланая матрица жесткости
    double h_x, h_y, h_z;
    h_x = nodes[node_num[2]].x - nodes[node_num[0]].x;
    h_y = nodes[node_num[6]].y - nodes[node_num[0]].y;
    h_z = nodes[node_num[18]].z - nodes[node_num[0]].z;
    vector<vector<double>> loc_matr_stiff(27);
    for (int i = 0; i < 27; i++) {
        loc_matr_stiff[i].resize(27);
    }
    for (int i = 0; i < 27; i++) {
        for (int j = 0; j < 27; j++) {
            loc_matr_stiff[i][j] = lambda(num_sub) *
                ((h_y * h_z * G[i % 3][j % 3] * M[(i / 3) % 3][(j / 3) % 3] * M[i / 9][j / 9]) / (2700.0 * h_x) +
                    (h_x * h_z * M[i % 3][j % 3] * G[(i / 3) % 3][(j / 3) % 3] * M[i / 9][j / 9]) / (2700.0 * h_y) +
                    (h_x * h_y * M[i % 3][j % 3] * M[(i / 3) % 3][(j / 3) % 3] * G[i / 9][j / 9]) / (2700.0 * h_z));
        }
    }
    AddLocalMatr(node_num, loc_matr_stiff);
}

void localVecB(int num_sub, vector<int> node_num) {
    double h_x, h_y, h_z;
    h_x = nodes[node_num[2]].x - nodes[node_num[0]].x;
    h_y = nodes[node_num[6]].y - nodes[node_num[0]].y;
    h_z = nodes[node_num[18]].z - nodes[node_num[0]].z;
    vector<double> loc_vec_b(27, 0);
    for (int i = 0; i < 27; i++) {
        for (int j = 0; j < 27; j++) {
            loc_vec_b[i] += M[i % 3][j % 3] * M[(i / 3) % 3][(j / 3) % 3] * M[i / 9][j / 9];
        }
        loc_vec_b[i] *= h_x * h_y * h_z * f(num_sub, nodes[node_num[i]].x, nodes[node_num[i]].y) / (27000.0);
    }
    AddLocalVec(node_num, loc_vec_b);
}
#pragma endregion

#pragma region Построение глобальной матрицы и вектора
void BuildMatrA_VecB() {
    for (int i = 0; i < NUM_SPLIT_X * NUM_SPLIT_Y * NUM_SPLIT_Z; i++) {

        auto result = find(begin(ident_fict), end(ident_fict), array_p[i].second[0]);
        if (result != end(ident_fict)) i += i_ident_fict[result - begin(ident_fict)];

        vector<int> node_num(27);
        for (int j = 0; j < 27; j++) {
            node_num[j] = array_p[i].second[j];
        }
        int num_sub = array_p[i].first;
        LocalMatrStiff(num_sub, node_num);
        localVecB(num_sub, node_num);
    }
}
#pragma endregion

#pragma region Тестирование
bool InVector(int num) {
    for (int i = 0; i < fict_nodes.size(); i++) {
        if (num == fict_nodes[i]) return true;
    }
    return false;
}

void Test() {
    vector<double> q_u(NUM_NODES, 0);
    for (int i = 0; i < NUM_NODES; i++) {
        q_u[i] = u_g(nodes[i].x, nodes[i].y, nodes[i].z);
    }
    double norm_vec_err = 0, norm_vec_q_u = 0; // норма вектора погрешности и q_u
    for (int i = 0; i < NUM_NODES; i++) {
        if (InVector(i)) continue;
        norm_vec_err += (q[i] - q_u[i]) * (q[i] - q_u[i]);
        norm_vec_q_u += (q_u[i]) * (q_u[i]);
    }
    cout << endl;
    cout << "Относительная норма вектора погрешности полученного решения:" << endl;
    cout << sqrt(norm_vec_err) / sqrt(norm_vec_q_u) << endl << endl;
}
#pragma endregion

#pragma region Возврат значения функции в любой точке
double BasicFunc(int i, double input_point, double point_1, double point_2, double point_3) {
    double h = point_3 - point_1;
    switch (i)
    {
    case 0:
        return 2 * (input_point - point_2) * (input_point - point_3) / (h * h);
        break;
    case 1:
        return -4 * (input_point - point_1) * (input_point - point_3) / (h * h);
        break;
    case 2:
        return 2 * (input_point - point_1) * (input_point - point_2) / (h * h);
        break;
    }
}

int GetNumEndEl(Coord3 point) { // получить номер конечного элемента
    for (int i = 0; i < array_p.size(); i++) {
        double min_x = nodes[array_p[i].second[0]].x; double min_y = nodes[array_p[i].second[0]].y; double min_z = nodes[array_p[i].second[0]].z;  //берём значения координат параллелепипеда
        double max_x = nodes[array_p[i].second[2]].x; double max_y = nodes[array_p[i].second[6]].y; double max_z = nodes[array_p[i].second[18]].z;
        if ((point.x > min_x || fabs(point.x - min_x) < EPS) &&
            (max_x > point.x || fabs(max_x - point.x) < EPS) &&
            (point.y > min_y || fabs(point.y - min_y) < EPS) &&
            (max_y > point.y || fabs(max_y - point.y) < EPS) &&
            (point.z > min_z || fabs(point.z - min_z) < EPS) &&
            (max_z > point.z || fabs(max_z - point.z) < EPS))
            return i;
    }
}
double ResUInPoint(Coord3 point) {
    int num_end_el = GetNumEndEl(point);
    double res = 0;

    for (int i = 0; i < array_p[num_end_el].second.size(); i++) {

        auto result = find(begin(ident_fict), end(ident_fict), array_p[i].second[0]);
        if (result != end(ident_fict)) i += i_ident_fict[result - begin(ident_fict)];

        int ind = array_p[num_end_el].second[i];

        int ind_x_left = array_p[num_end_el].second[0];
        int ind_x_middle = array_p[num_end_el].second[1];
        int ind_x_right = array_p[num_end_el].second[2];

        int ind_y_left = array_p[num_end_el].second[0];
        int ind_y_middle = array_p[num_end_el].second[3];
        int ind_y_right = array_p[num_end_el].second[6];

        int ind_z_left = array_p[num_end_el].second[0];
        int ind_z_middle = array_p[num_end_el].second[9];
        int ind_z_right = array_p[num_end_el].second[18];

        res += q[ind] * BasicFunc(i % 3, point.x, nodes[ind_x_left].x, nodes[ind_x_middle].x, nodes[ind_x_right].x) *
            BasicFunc((i / 3) % 3, point.y, nodes[ind_y_left].y, nodes[ind_y_middle].y, nodes[ind_y_right].y) *
            BasicFunc(i / 9, point.z, nodes[ind_z_left].z, nodes[ind_z_middle].z, nodes[ind_z_right].z);
    }
    return res;
}
#pragma endregion

#pragma region Выводы нужных значений
void Output_result() {
    ofstream f_result;
    f_result.imbue(locale("Russian"));
    f_result.open("f_result.txt");
    for (int i = 0; i < NUM_NODES; i++) {
        f_result << setprecision(15) << q[i] << endl;
    }
    f_result.close();
}

void Output_u() {
    ofstream f_result;
    f_result.imbue(locale("Russian"));
    f_result.open("u.txt");
    for (int i = 0; i < NUM_NODES; i++) {
        f_result << setprecision(15) << u_g(nodes[i].x, nodes[i].y, nodes[i].z) << endl;
    }
    f_result.close();
}

void Output_x() {
    ofstream f_result;
    f_result.imbue(locale("Russian"));
    f_result.open("x.txt");
    for (int i = 0; i < NUM_NODES; i++) {
        f_result << setprecision(15) << nodes[i].x << endl;
    }
    f_result.close();
}

void Output_y() {
    ofstream f_result;
    f_result.imbue(locale("Russian"));
    f_result.open("y.txt");
    for (int i = 0; i < NUM_NODES; i++) {
        f_result << setprecision(15) << nodes[i].y << endl;
    }
    f_result.close();
}

void Output_z() {
    ofstream f_result;
    f_result.imbue(locale("Russian"));
    f_result.open("z.txt");
    for (int i = 0; i < NUM_NODES; i++) {
        f_result << setprecision(15) << nodes[i].z << endl;
    }
    f_result.close();
}
#pragma endregion

#pragma region Расчет объема смеси, небаланса
double DerivativeLocalBasicFuncs(int i, double psi) {
    switch (i)
    {
    case 0:
        return 4 * psi - 3;
        break;
    case 1:
        return -8 * psi + 4;
        break;
    case 2:
        return 4 * psi - 1;
        break;
    }
}

double F(vector<int> node_num, int num_sub) {
    double res = 0;
    double h_x, h_y, h_z;
    h_x = nodes[node_num[2]].x - nodes[node_num[0]].x;
    h_y = nodes[node_num[6]].y - nodes[node_num[0]].y;
    h_z = nodes[node_num[18]].z - nodes[node_num[0]].z;
    vector<double> line_matr_mass(27);
    line_matr_mass = { 13.0 / 140.0, 1.0 / 21.0, -1.0 / 140.0, 1.0 / 21.0,
                      4.0 / 105.0, -2.0 / 105.0, -1.0 / 140.0, -2.0 / 105.0,
                      -1.0 / 140.0, 1.0 / 21.0, 4.0 / 105.0, -2.0 / 105.0,
                      4.0 / 105.0, 16.0 / 35.0, 4.0 / 105.0, -2.0 / 105.0,
                      4.0 / 105.0, 1.0 / 21.0, -1.0 / 140.0, -2.0 / 105.0,
                      -1.0 / 140.0, -2.0 / 105.0, 4.0 / 105.0, 1.0 / 21.0,
                      -1.0 / 140.0, 1.0 / 21.0, 13.0 / 140.0 };
    for (int i = 0; i < 27; i++) {
        res += line_matr_mass[i] * f(num_sub, nodes[node_num[i]].x, nodes[node_num[i]].y);
    }
    return res * h_x * h_y * h_z;
}

double CalcMixtureVFaceX(vector<int> node_num, int num_sub, int n) {
    double h_x, h_y, h_z;
    h_x = nodes[node_num[2]].x - nodes[node_num[0]].x;
    h_y = nodes[node_num[6]].y - nodes[node_num[0]].y;
    h_z = nodes[node_num[18]].z - nodes[node_num[0]].z;
    vector<double> line_matr_mass(27);
    line_matr_mass = { 4, 4, 4, 2, 2, 2, -1, -1, -1, 2, 2, 2, 16, 16, 16,
                       2, 2, 2, -1, -1, -1, 2, 2, 2, 4, 4, 4 };
    double res = 0;
    for (int i = 0; i < 27; i++) {
        if (n == -1) {
            double psi = 0;
            res += q[node_num[i]] * line_matr_mass[i] * DerivativeLocalBasicFuncs(i % 3, psi);
        }
        else {
            double psi = 1;
            res += -q[node_num[i]] * line_matr_mass[i] * DerivativeLocalBasicFuncs(i % 3, psi);
        }
    }
    res *= h_y * h_z * lambda(num_sub) / (30.0 * h_x);
    return res;
}

double CalcMixtureVFaceY(vector<int> node_num, int num_sub, int n) {
    double h_x, h_y, h_z;
    h_x = nodes[node_num[2]].x - nodes[node_num[0]].x;
    h_y = nodes[node_num[6]].y - nodes[node_num[0]].y;
    h_z = nodes[node_num[18]].z - nodes[node_num[0]].z;
    vector<double> line_matr_mass(27);
    line_matr_mass = { 4, 2, -1, 4, 2, -1, 4, 2, -1, 2, 16, 2, 2, 16, 2,
                       2, 16, 2, -1, 2, 4, -1, 2, 4, -1, 2, 4 };
    double res = 0;
    for (int i = 0; i < 27; i++) {
        if (n == -1) {
            double psi = 0;
            res += q[node_num[i]] * line_matr_mass[i] * DerivativeLocalBasicFuncs((i / 3) % 3, psi);
        }
        else {
            double psi = 1;
            res += -q[node_num[i]] * line_matr_mass[i] * DerivativeLocalBasicFuncs((i / 3) % 3, psi);
        }
    }
    res *= h_x * h_z * lambda(num_sub) / (30.0 * h_y);
    return res;
}

double CalcMixtureVFaceZ(vector<int> node_num, int num_sub, int n) {
    double h_x, h_y, h_z;
    h_x = nodes[node_num[2]].x - nodes[node_num[0]].x;
    h_y = nodes[node_num[6]].y - nodes[node_num[0]].y;
    h_z = nodes[node_num[18]].z - nodes[node_num[0]].z;
    vector<double> line_matr_mass(27);
    line_matr_mass = { 4, 2, -1, 2, 16, 2, -1, 2, 4, 4, 2, -1, 2, 16, 2,
                       -1, 2, 4, 4, 2, -1, 2, 16, 2, -1, 2, 4 };
    double res = 0;
    for (int i = 0; i < 27; i++) {
        if (n == -1) {
            double psi = 0;
            res += q[node_num[i]] * line_matr_mass[i] * DerivativeLocalBasicFuncs(i / 9, psi);
        }
        else {
            double psi = 1;
            res += -q[node_num[i]] * line_matr_mass[i] * DerivativeLocalBasicFuncs(i / 9, psi);
        }
    }
    res *= h_x * h_y * lambda(num_sub) / (30.0 * h_z);
    return res;
}

vector<int> GetNodeNum(int num_end_el) { // создать массив номеров узлов конечного элемента
    double res = 0;
    vector<int> node_num(27);
    for (int j = 0; j < 27; j++) {
        node_num[j] = array_p[num_end_el].second[j];
    }
    return node_num;
}

void GetArrayFictEndEl() { // массив фиктивных конечных элементов
    for (int i = 0; i < NUM_SPLIT_X * NUM_SPLIT_Y * NUM_SPLIT_Z; i++) {
        auto result = find(begin(ident_fict), end(ident_fict), array_p[i].second[0]);
        if (result != end(ident_fict)) {
            int num_fict_el = i_ident_fict[result - begin(ident_fict)];
            for (int l = 0; l < num_fict_el; l++) {
                fict_el.push_back(i + l);
            }
            i += num_fict_el - 1;
            continue;
        }
    }
}

bool IsFictEl(int num_end_el) { // является ли конечный элемент фиктивным
    for (int i = 0; i < fict_el.size(); i++) {
        if (fict_el[i] == num_end_el)
            return true;
    }
    return false;
}

int GetNumFace2ZP(vector<int> _face_2_zp) {
    for (int i = 0; i < face_2_zp.size(); i++) {
        int a = 0;
        for (int j = 0; j < 9; j++) {
            if (face_2_zp[i].second[j] == _face_2_zp[j]) {
                a++;
            }
        }
        if (a == 9) {
            return i;
        }
    }
}

void CreateFaceV() { // генерация массива граней со значением потока
    // порядок заполнения: 
    // сначала грани по y, потом по x (слева - направо; от ближайшей до дальней (если смотреть в сечении xy то снизу вверх)), потом вверх по z

    GetArrayFictEndEl();


    vector<int> node_num(27); int num_sub;
    vector<int> el_for_face; // конечный элемент, соответствующий грани
    bool left = true;
    bool face_y = true;

    int num_faces = ((2 * NUM_SPLIT_X + 1) * NUM_SPLIT_Y + NUM_SPLIT_X) * NUM_SPLIT_Z; // пока что берем только x и y грани
    // (верхняя и нижняя z известна и не нужна фактически, а внутренние по z практически нулевые - можно не учитывать (скважина полностью проходит по глубине z))
    int start_el = 0;



    for (int i = 0; i < NUM_SPLIT_Z; i++) {
        //самая передняя грань по y
        for (int j = 0; j < NUM_SPLIT_X; j++) {
            node_num = GetNodeNum(start_el + j);
            num_sub = array_p[start_el + j].first;
            face_V.push_back(-CalcMixtureVFaceY(node_num, num_sub, -1));
            el_for_face.push_back(start_el + j);
        }
        //самая передняя грань по y

        //самая левая грань по x
        node_num = GetNodeNum(start_el);
        num_sub = array_p[start_el].first;
        face_V.push_back(-CalcMixtureVFaceX(node_num, num_sub, -1));
        el_for_face.push_back(start_el);
        //самая левая грань по x

        //внутренние грани по x
        for (int j = 0; j < NUM_SPLIT_X - 1; j++) {
            double lambda_left_el, lambda_right_el;
            node_num = GetNodeNum(start_el + j);
            num_sub = array_p[start_el + j].first;
            double right_V_left_el = CalcMixtureVFaceX(node_num, num_sub, 1); // правый поток левого конечного элемента
            lambda_left_el = lambda(num_sub);

            node_num = GetNodeNum(start_el + j + 1);
            num_sub = array_p[start_el + j + 1].first;
            double left_V_right_el = CalcMixtureVFaceX(node_num, num_sub, -1); // левый поток правого конечного элемента
            lambda_right_el = lambda(num_sub);

            double res = (right_V_left_el * lambda_right_el - left_V_right_el * lambda_left_el) / (lambda_right_el + lambda_left_el);
            face_V.push_back(res);
            el_for_face.push_back(start_el + j);
        }
        //внутренние грани по x

        //самая правая грань по x
        node_num = GetNodeNum(start_el + NUM_SPLIT_X - 1);
        num_sub = array_p[start_el + NUM_SPLIT_X - 1].first;
        face_V.push_back(CalcMixtureVFaceX(node_num, num_sub, 1));
        el_for_face.push_back(start_el + NUM_SPLIT_X - 1);
        //самая правая грань по x

        start_el += NUM_SPLIT_X;
        int num_face_2_zp = 0;

        //основной цикл
        for (int j = 1; j < NUM_SPLIT_Y; j++) {

            //внутренние грани по y
            for (int k = 0; k < NUM_SPLIT_X; k++) {

                if (IsFictEl(start_el + k) && !IsFictEl(start_el + k - NUM_SPLIT_X)) {
                    node_num = GetNodeNum(start_el + k);
                    vector<int> _face_2_zp(9);
                    _face_2_zp[0] = node_num[0]; _face_2_zp[1] = node_num[1]; _face_2_zp[2] = node_num[2];
                    _face_2_zp[3] = node_num[9]; _face_2_zp[4] = node_num[10]; _face_2_zp[5] = node_num[11];
                    _face_2_zp[6] = node_num[18]; _face_2_zp[7] = node_num[19]; _face_2_zp[8] = node_num[20];
                    int num_face_2_zp = GetNumFace2ZP(_face_2_zp);
                    face_V.push_back(-theta_V(num_face_2_zp)[1]);
                    num_faces_zp.push_back({ face_V.size() - 1 , 0 });
                    num_face_2_zp++;
                }
                else if (!IsFictEl(start_el + k) && IsFictEl(start_el + k - NUM_SPLIT_X)) {
                    node_num = GetNodeNum(start_el + k);
                    vector<int> _face_2_zp(9);
                    _face_2_zp[0] = node_num[0]; _face_2_zp[1] = node_num[1]; _face_2_zp[2] = node_num[2];
                    _face_2_zp[3] = node_num[9]; _face_2_zp[4] = node_num[10]; _face_2_zp[5] = node_num[11];
                    _face_2_zp[6] = node_num[18]; _face_2_zp[7] = node_num[19]; _face_2_zp[8] = node_num[20];
                    int num_face_2_zp = GetNumFace2ZP(_face_2_zp);
                    face_V.push_back(theta_V(num_face_2_zp)[1]);
                    num_faces_zp.push_back({ face_V.size() - 1 , 1 });
                    num_face_2_zp++;
                }
                else if (IsFictEl(start_el + k) && IsFictEl(start_el + k - NUM_SPLIT_X)) {
                    face_V.push_back(0);
                    fict_faces.push_back(face_V.size() - 1);
                }
                else {
                    double lambda_front_el, lambda_back_el;
                    node_num = GetNodeNum(start_el + k - NUM_SPLIT_X);
                    num_sub = array_p[start_el + k - NUM_SPLIT_X].first;
                    double back_V_front_el = CalcMixtureVFaceY(node_num, num_sub, 1); // задний поток переднего конечного элемента
                    lambda_front_el = lambda(num_sub);

                    node_num = GetNodeNum(start_el + k);
                    num_sub = array_p[start_el + k].first;
                    double front_V_back_el = CalcMixtureVFaceY(node_num, num_sub, -1); // передний поток заднего конечного элемента
                    lambda_back_el = lambda(num_sub);

                    double res = (back_V_front_el * lambda_back_el - front_V_back_el * lambda_front_el) / (lambda_back_el + lambda_front_el);
                    face_V.push_back(res);
                    el_for_face.push_back(start_el + k);
                }
            }
            //внутренние грани по y

            //самая левая грань по x
            node_num = GetNodeNum(start_el);
            num_sub = array_p[start_el].first;
            face_V.push_back(-CalcMixtureVFaceX(node_num, num_sub, -1));
            el_for_face.push_back(start_el);
            //самая левая грань по x

            //внутренние грани по x
            for (int k = 0; k < NUM_SPLIT_X - 1; k++) {
                if (!IsFictEl(start_el + k) && IsFictEl(start_el + k + 1)) {
                    node_num = GetNodeNum(start_el + k);
                    vector<int> _face_2_zp(9);
                    _face_2_zp[0] = node_num[2]; _face_2_zp[1] = node_num[5]; _face_2_zp[2] = node_num[8];
                    _face_2_zp[3] = node_num[11]; _face_2_zp[4] = node_num[14]; _face_2_zp[5] = node_num[17];
                    _face_2_zp[6] = node_num[20]; _face_2_zp[7] = node_num[23]; _face_2_zp[8] = node_num[26];
                    int num_face_2_zp = GetNumFace2ZP(_face_2_zp);
                    face_V.push_back(-theta_V(num_face_2_zp)[1]);
                    num_faces_zp.push_back({ face_V.size() - 1 , 2 });
                    num_face_2_zp++;
                }
                else if (IsFictEl(start_el + k) && !IsFictEl(start_el + k + 1)) {
                    node_num = GetNodeNum(start_el + k);
                    vector<int> _face_2_zp(9);
                    _face_2_zp[0] = node_num[2]; _face_2_zp[1] = node_num[5]; _face_2_zp[2] = node_num[8];
                    _face_2_zp[3] = node_num[11]; _face_2_zp[4] = node_num[14]; _face_2_zp[5] = node_num[17];
                    _face_2_zp[6] = node_num[20]; _face_2_zp[7] = node_num[23]; _face_2_zp[8] = node_num[26];
                    int num_face_2_zp = GetNumFace2ZP(_face_2_zp);
                    face_V.push_back(theta_V(num_face_2_zp)[1]);
                    num_faces_zp.push_back({ face_V.size() - 1 , 3 });
                    num_face_2_zp++;
                }
                else if (IsFictEl(start_el + k) && IsFictEl(start_el + k + 1)) {
                    face_V.push_back(0);
                    fict_faces.push_back(face_V.size() - 1);
                }
                else {
                    double lambda_left_el, lambda_right_el;
                    node_num = GetNodeNum(start_el + k);
                    num_sub = array_p[start_el + k].first;
                    double right_V_left_el = CalcMixtureVFaceX(node_num, num_sub, 1); // правый поток левого конечного элемента
                    lambda_left_el = lambda(num_sub);

                    node_num = GetNodeNum(start_el + k + 1);
                    num_sub = array_p[start_el + k + 1].first;
                    double left_V_right_el = CalcMixtureVFaceX(node_num, num_sub, -1); // левый поток правого конечного элемента
                    lambda_right_el = lambda(num_sub);

                    double res = (right_V_left_el * lambda_right_el - left_V_right_el * lambda_left_el) / (lambda_right_el + lambda_left_el);
                    face_V.push_back(res);
                    el_for_face.push_back(start_el + k);
                }
            }
            //внутренние грани по x

            //самая правая грань по x
            node_num = GetNodeNum(start_el + NUM_SPLIT_X - 1);
            num_sub = array_p[start_el + NUM_SPLIT_X - 1].first;
            face_V.push_back(CalcMixtureVFaceX(node_num, num_sub, 1));
            el_for_face.push_back(start_el + NUM_SPLIT_X - 1);
            //самая правая грань по x

            if (j == NUM_SPLIT_Y - 1) continue;
            start_el += NUM_SPLIT_X;
        }
        //основной цикл

        //самая дальняя грань по y
        for (int j = 0; j < NUM_SPLIT_X; j++) {
            node_num = GetNodeNum(start_el + j);
            num_sub = array_p[start_el + j].first;
            face_V.push_back(CalcMixtureVFaceY(node_num, num_sub, 1));
            el_for_face.push_back(start_el + j);
        }
        //самая дальняя грань по y

        start_el += NUM_SPLIT_X;
    }
}

void ClearAndGenPortMatrB() { // отчистить старое и создать новый портрет для матрицы B (балансировка)
    ia.clear();
    ja.clear();
    aal.clear();
    di.clear();
    b.clear();
    L_sq.clear();
    di_sq.clear();
    list_face.resize(face_V.size());
    bool face_y = true;
    int num_face_layer = face_V.size() / NUM_SPLIT_Z; // кол-во ребер в слое (в плоскости xy в одном разбиении по z)
    for (int k = 0; k < NUM_SPLIT_Z; k++) {
        for (int i = 0; i < num_face_layer;) {
            if (face_y) {
                for (int j = 0; j < NUM_SPLIT_X; j++) {
                    if (i != 0)
                        list_face[k * num_face_layer + i + j].push_back(k * num_face_layer + i + j - 2 * NUM_SPLIT_X - 1);

                }
                for (int j = 0; j < NUM_SPLIT_X; j++) {
                    if (i != 0) {
                        list_face[k * num_face_layer + i + j].push_back(k * num_face_layer + i + j - NUM_SPLIT_X - 1);
                        list_face[k * num_face_layer + i + j].push_back(k * num_face_layer + i + j - NUM_SPLIT_X);
                    }
                }
                i += NUM_SPLIT_X;
                face_y = !face_y;
            }
            else {
                for (int j = 0; j < NUM_SPLIT_X + 1; j++) {
                    if (j == 0) {
                        list_face[k * num_face_layer + i + j].push_back(k * num_face_layer + i + j - NUM_SPLIT_X);
                    }
                    else if (j == NUM_SPLIT_X) {
                        list_face[k * num_face_layer + i + j].push_back(k * num_face_layer + i + j - NUM_SPLIT_X - 1);
                    }
                    else {
                        list_face[k * num_face_layer + i + j].push_back(k * num_face_layer + i + j - NUM_SPLIT_X - 1);
                        list_face[k * num_face_layer + i + j].push_back(k * num_face_layer + i + j - NUM_SPLIT_X);
                    }
                }
                for (int j = 0; j < NUM_SPLIT_X + 1; j++) {
                    if (j == 0) {
                        continue;
                    }
                    list_face[k * num_face_layer + i + j].push_back(k * num_face_layer + i + j - 1);
                }
                i += NUM_SPLIT_X + 1;
                face_y = !face_y;
            }
        }
        face_y = !face_y;
    }


    ia.resize(face_V.size() + 1, 0);
    for (int i = 2; i <= face_V.size(); i++) {
        ia[i] = ia[i - 1] + list_face[i - 1].size();
    }
    ja.resize(ia[face_V.size()]);
    auto iter = ja.begin();
    for (int i = 0; i < face_V.size(); i++) {
        copy(list_face[i].begin(), list_face[i].end(), iter);
        iter += list_face[i].size();
    }
    //aal.resize(ia[face_V.size()]);
    di.resize(face_V.size());
    b.resize(face_V.size());
    L_sq.resize(ia[face_V.size()]);
    di_sq.resize(face_V.size());
}

int FindFaceZP(int num_face) { // найти грань в массиве с номерами граней с краевыми 2 рода зон перфорации
    for (int i = 0; i < num_faces_zp.size(); i++) {
        if (num_faces_zp[i].first == num_face)
            return num_faces_zp[i].second;
    }
    return -1;
}

bool FindFaceBool(int num_face) {
    for (int i = 0; i < num_faces_zp.size(); i++) {
        if (num_faces_zp[i].first == num_face)
            return true;
    }
    return false;
}

void GenMatrBAndVecB() {

    auto max_it = max_element(face_V.begin(), face_V.end(),
        [](double a, double b) { return fabs(a) < fabs(b); });
    double max_V = fabs(*max_it);

    double beta = 1e+8;
    double alpha;
    int num_face_layer = face_V.size() / NUM_SPLIT_Z; // кол-во ребер в слое (в плоскости xy в одном разбиении по z)

    for (int k = 0; k < NUM_SPLIT_Z; k++) {

        // начальные компоненты вектора
        for (int i = 0; i < NUM_SPLIT_X; i++) { // не работает для разбиений по z
            int S_g_1, S_g_2, S_g_3, S_g_4;

            if (fabs(face_V[k * num_face_layer + i]) > (max_V * 1e-4))
                alpha = 1. / max_V;
            else
                alpha = 1. / (max_V * 1e-4);

            S_g_1 = (face_V[k * num_face_layer + i] < 0) ? 1 : -1;
            S_g_2 = (face_V[k * num_face_layer + i + NUM_SPLIT_X] < 0) ? 1 : -1;
            S_g_3 = (face_V[k * num_face_layer + i + NUM_SPLIT_X + 1] < 0) ? -1 : 1;
            S_g_4 = (face_V[k * num_face_layer + i + 2 * NUM_SPLIT_X + 1] < 0) ? -1 : 1;
            di[k * num_face_layer + i] = beta + alpha;
            b[k * num_face_layer + i] -= beta * S_g_1 * (S_g_1 * fabs(face_V[k * num_face_layer + i]) + S_g_2 * fabs(face_V[k * num_face_layer + i + NUM_SPLIT_X]) +
                S_g_3 * fabs(face_V[k * num_face_layer + i + NUM_SPLIT_X + 1]) + S_g_4 * fabs(face_V[k * num_face_layer + i + 2 * NUM_SPLIT_X + 1]));
        }

        for (int i = 0; i < num_face_layer; i++) {
            // расчет диагонали
            if (fabs(face_V[k * num_face_layer + i]) > (max_V * 1e-4))
                alpha = 1. / max_V;
            else
                alpha = 1. / (max_V * 1e-4);

            int S_g_1, S_g_2;

            // расчет компонент матрицы
            for (int j = 0; j < list_face[k * num_face_layer + i].size(); j++) {

                if (list_face[k * num_face_layer + i].size() == 1 || (k * num_face_layer + i) == list_face[k * num_face_layer + i][list_face[k * num_face_layer + i].size() - 1] + 1) { // грани по x
                    if (list_face[k * num_face_layer + i].size() == 1) { // самый левый x
                        di[k * num_face_layer + i] = beta + alpha;
                        S_g_1 = (face_V[k * num_face_layer + i] < 0) ? 1 : -1;
                        S_g_2 = (face_V[list_face[k * num_face_layer + i][j]] < 0) ? 1 : -1;

                        // расчет компонент вектора правой части
                        int dif_S_g_1, dif_S_g_2;
                        dif_S_g_1 = (face_V[k * num_face_layer + i + 1] < 0) ? -1 : 1;
                        dif_S_g_2 = (face_V[list_face[k * num_face_layer + i][j] + 2 * NUM_SPLIT_X + 1] < 0) ? -1 : 1;

                        b[k * num_face_layer + i] -= beta * S_g_1 * (S_g_1 * fabs(face_V[k * num_face_layer + i]) + S_g_2 * fabs(face_V[list_face[k * num_face_layer + i][j]]) +
                            dif_S_g_1 * fabs(face_V[k * num_face_layer + i + 1]) + dif_S_g_2 * fabs(face_V[list_face[k * num_face_layer + i][j] + 2 * NUM_SPLIT_X + 1]));
                    }
                    else if (list_face[k * num_face_layer + i].size() == 2) { // самый правый x
                        di[k * num_face_layer + i] = beta + alpha;
                        S_g_1 = (face_V[k * num_face_layer + i] < 0) ? -1 : 1;
                        S_g_2 = (face_V[list_face[k * num_face_layer + i][j]] < 0) ? 1 : -1;

                        // расчет компонент вектора правой части
                        if (j == 0) {
                            int dif_S_g_1, dif_S_g_2;
                            dif_S_g_1 = (face_V[k * num_face_layer + i - 1] < 0) ? 1 : -1;
                            dif_S_g_2 = (face_V[list_face[k * num_face_layer + i][j] + 2 * NUM_SPLIT_X + 1] < 0) ? -1 : 1;

                            b[k * num_face_layer + i] -= beta * S_g_1 * (S_g_1 * fabs(face_V[k * num_face_layer + i]) + S_g_2 * fabs(face_V[list_face[k * num_face_layer + i][j]]) +
                                dif_S_g_1 * fabs(face_V[k * num_face_layer + i - 1]) + dif_S_g_2 * fabs(face_V[list_face[k * num_face_layer + i][j] + 2 * NUM_SPLIT_X + 1]));
                        }
                    }
                    else {
                        di[k * num_face_layer + i] = 2 * beta + alpha;
                        if (j == 0 || j == 2) {
                            S_g_1 = (face_V[k * num_face_layer + i] < 0) ? -1 : 1;
                            S_g_2 = (face_V[list_face[k * num_face_layer + i][j]] < 0) ? 1 : -1;

                            // расчет компонент вектора правой части
                            if (j == 0) {
                                int dif_S_g_1, dif_S_g_2;
                                dif_S_g_1 = (face_V[k * num_face_layer + i - 1] < 0) ? 1 : -1;
                                dif_S_g_2 = (face_V[list_face[k * num_face_layer + i][j] + 2 * NUM_SPLIT_X + 1] < 0) ? -1 : 1;

                                b[k * num_face_layer + i] -= beta * S_g_1 * (S_g_1 * fabs(face_V[k * num_face_layer + i]) + S_g_2 * fabs(face_V[list_face[k * num_face_layer + i][j]]) +
                                    dif_S_g_1 * fabs(face_V[k * num_face_layer + i - 1]) + dif_S_g_2 * fabs(face_V[list_face[k * num_face_layer + i][j] + 2 * NUM_SPLIT_X + 1]));
                            }
                        }
                        else {
                            S_g_1 = (face_V[k * num_face_layer + i] < 0) ? 1 : -1;
                            S_g_2 = (face_V[list_face[k * num_face_layer + i][j]] < 0) ? 1 : -1;

                            // расчет компонент вектора правой части
                            int dif_S_g_1, dif_S_g_2;
                            dif_S_g_1 = (face_V[k * num_face_layer + i + 1] < 0) ? -1 : 1;
                            dif_S_g_2 = (face_V[list_face[k * num_face_layer + i][j] + 2 * NUM_SPLIT_X + 1] < 0) ? -1 : 1;

                            b[k * num_face_layer + i] -= beta * S_g_1 * (S_g_1 * fabs(face_V[k * num_face_layer + i]) + S_g_2 * fabs(face_V[list_face[k * num_face_layer + i][j]]) +
                                dif_S_g_1 * fabs(face_V[k * num_face_layer + i + 1]) + dif_S_g_2 * fabs(face_V[list_face[k * num_face_layer + i][j] + 2 * NUM_SPLIT_X + 1]));
                        }
                    }
                }
                else { // грани по y  
                    S_g_1 = (face_V[k * num_face_layer + i] < 0) ? -1 : 1;
                    if (j == 0 || j == 1) {
                        S_g_2 = (face_V[list_face[k * num_face_layer + i][j]] < 0) ? 1 : -1;
                    }
                    else {
                        S_g_2 = (face_V[list_face[k * num_face_layer + i][j]] < 0) ? -1 : 1;
                    }

                    // расчет компонент вектора правой части для граней в конце
                    if (i >= (2 * NUM_SPLIT_X + 1) * NUM_SPLIT_Y && j == 0) {
                        di[k * num_face_layer + i] = beta + alpha;
                        int dif_S_g_1, dif_S_g_2;
                        dif_S_g_1 = (face_V[list_face[k * num_face_layer + i][j] + NUM_SPLIT_X] < 0) ? 1 : -1;
                        dif_S_g_2 = (face_V[list_face[k * num_face_layer + i][j] + NUM_SPLIT_X + 1] < 0) ? -1 : 1;

                        b[k * num_face_layer + i] -= beta * S_g_1 * (S_g_1 * fabs(face_V[k * num_face_layer + i]) + S_g_2 * fabs(face_V[list_face[k * num_face_layer + i][j]]) +
                            dif_S_g_1 * fabs(face_V[list_face[k * num_face_layer + i][j] + NUM_SPLIT_X]) + dif_S_g_2 * fabs(face_V[list_face[k * num_face_layer + i][j] + NUM_SPLIT_X + 1]));
                    }
                    else if (j == 0) { // не в конце
                        di[k * num_face_layer + i] = 2 * beta + alpha;
                        int dif_S_g_1, dif_S_g_2;
                        dif_S_g_1 = (face_V[list_face[k * num_face_layer + i][j] + NUM_SPLIT_X] < 0) ? 1 : -1;
                        dif_S_g_2 = (face_V[list_face[k * num_face_layer + i][j] + NUM_SPLIT_X + 1] < 0) ? -1 : 1;

                        b[k * num_face_layer + i] -= beta * S_g_1 * (S_g_1 * fabs(face_V[k * num_face_layer + i]) + S_g_2 * fabs(face_V[list_face[k * num_face_layer + i][j]]) +
                            dif_S_g_1 * fabs(face_V[list_face[k * num_face_layer + i][j] + NUM_SPLIT_X]) + dif_S_g_2 * fabs(face_V[list_face[k * num_face_layer + i][j] + NUM_SPLIT_X + 1]));

                        int d_S_g_1 = (face_V[k * num_face_layer + i] < 0) ? 1 : -1;
                        int d_S_g_2 = (face_V[k * num_face_layer + i + 2 * NUM_SPLIT_X + 1] < 0) ? -1 : 1;
                        dif_S_g_1 = (face_V[k * num_face_layer + i + NUM_SPLIT_X] < 0) ? 1 : -1;
                        dif_S_g_2 = (face_V[k * num_face_layer + i + NUM_SPLIT_X + 1] < 0) ? -1 : 1;

                        b[k * num_face_layer + i] -= beta * d_S_g_1 * (d_S_g_1 * fabs(face_V[k * num_face_layer + i]) + d_S_g_2 * fabs(face_V[k * num_face_layer + i + 2 * NUM_SPLIT_X + 1]) +
                            dif_S_g_1 * fabs(face_V[k * num_face_layer + i + NUM_SPLIT_X]) + dif_S_g_2 * fabs(face_V[k * num_face_layer + i + NUM_SPLIT_X + 1]));
                    }
                }
                double a = beta * S_g_1 * S_g_2;
                aal.push_back(beta * S_g_1 * S_g_2);
            }
            if (FindFaceBool(i)) {
                di[k * num_face_layer + i] = beta + alpha;
            }
        }
    }
}

bool FindFictFaceBool(int num_face) {
    for (int i = 0; i < fict_faces.size(); i++) {
        if (fict_faces[i] == num_face)
            return true;
    }
    return false;
}

void ConsiderKnownFlows(int n) { // учет известных потоков
    b[n] = 0;
    di[n] = 1;
    for (int i = ia[n]; i < ia[n + 1]; i++) {
        int _i = ja[i];
        if (FindFaceBool(_i)) {
            aal[i] = 0;
            continue;
        }
        b[_i] -= b[n] * aal[i];
        aal[i] = 0;
    }
    for (int i = n; i < face_V.size(); i++) {
        int k = 0;
        for (int j = ia[i]; j < ia[i + 1]; j++) {
            if (ja[j] == n) {
                if (FindFaceBool(i)) {
                    aal[j] = 0;
                    continue;
                }
                b[i] -= b[n] * aal[j];
                aal[j] = 0;
            }
        }
    }
}

void ConsiderFictFace(int n) { // учет фиктивных потоков
    b[n] = 0;
    di[n] = 1;
    for (int i = ia[n]; i < ia[n + 1]; i++) {
        int _i = ja[i];
        if (FindFictFaceBool(_i)) {
            aal[i] = 0;
            continue;
        }
        b[_i] -= b[n] * aal[i];
        aal[i] = 0;
    }
    for (int i = n; i < face_V.size(); i++) {
        int k = 0;
        for (int j = ia[i]; j < ia[i + 1]; j++) {
            if (ja[j] == n) {
                if (FindFictFaceBool(i)) {
                    aal[j] = 0;
                    continue;
                }
                b[i] -= b[n] * aal[j];
                aal[j] = 0;
            }
        }
    }
}

double CalcSumNonBalance() {
    double res = 0;
    int start_face = 0;
    int num_face_layer = face_V.size() / NUM_SPLIT_Z; // кол-во ребер в слое (в плоскости xy в одном разбиении по z)
    for (int i = 0; i < NUM_SPLIT_Z; i++) {
        for (int j = 0; j < NUM_SPLIT_Y; j++) {
            for (int k = 0; k < NUM_SPLIT_X; k++) {
                int curr_el = i * NUM_SPLIT_X * NUM_SPLIT_Y + j * NUM_SPLIT_Y + k;

                if (IsFictEl(curr_el))
                    continue;

                int S_g_1, S_g_2, S_g_3, S_g_4;
                S_g_1 = (face_V[start_face + k] < 0) ? 1 : -1;
                S_g_2 = (face_V[start_face + k + NUM_SPLIT_X] < 0) ? 1 : -1;
                S_g_3 = (face_V[start_face + k + NUM_SPLIT_X + 1] < 0) ? -1 : 1;
                S_g_4 = (face_V[start_face + k + 2 * NUM_SPLIT_X + 1] < 0) ? -1 : 1;
                double a = S_g_1 * fabs(face_V[start_face + k]) + S_g_2 * fabs(face_V[start_face + k + NUM_SPLIT_X]) +
                    S_g_3 * fabs(face_V[start_face + k + NUM_SPLIT_X + 1]) + S_g_4 * fabs(face_V[start_face + k + 2 * NUM_SPLIT_X + 1]);
                res += a;
            }
            start_face += 2 * NUM_SPLIT_X + 1;
        }
        start_face += NUM_SPLIT_X;
    }
    return res;
}

void FindBalancedFlows() { // найти сбалансированные потоки
    for (int i = 0; i < face_V.size(); i++) {
        int S_g = (face_V[i] < 0) ? -1 : 1;
        face_V[i] = S_g * fabs(face_V[i]) + S_g * q_V[i];
    }
}

void BalancingFlows() { // балансировка потоков
    CreateFaceV();
    cout << endl << "Стартовый небаланс на всех элементах: " << CalcSumNonBalance() << endl;
    ClearAndGenPortMatrB();
    GenMatrBAndVecB();
    for (int i = 0; i < num_faces_zp.size(); i++) {
        ConsiderKnownFlows(num_faces_zp[i].first);
    }
    for (int i = 0; i < fict_faces.size(); i++) {
        ConsiderKnownFlows(fict_faces[i]);
    }
    q_V.resize(face_V.size(), 0);
    vector<double> r(face_V.size());
    vector<double> z(face_V.size());
    vector<double> Mult(face_V.size());
    vector<double> Az(face_V.size());
    int max_iter = 1000;
    double eps = 1e-15;
    cout << endl;
    MSG::LU_sq_MSG(q_V, r, z, Az, Mult, face_V.size(), eps, max_iter);
    FindBalancedFlows();
    cout << endl << "Небаланс на всех элементах после балансировки: " << CalcSumNonBalance() << endl;
}
#pragma endregion

//#pragma region Сравнение с аналитическим решением
//double AnalitP(double r) {
//    double res = _theta * r_w / lambda(0) * log(R / r) + P_g;
//
//    return res;
//}
//
//void VecAnalitP() {
//    double h_r;
//    double curr_r = 1;
//    double norm_true = 0, norm_err = 0;
//    int i = 23;
//    //int i = 45;
//    while (curr_r < R) {
//        h_r = nodes[i].x - nodes[i - 1].x;
//        double analyt_value_P = AnalitP(curr_r);
//        Coord3 p(curr_r, 0, 0.5);
//        double num_value_P = ResUInPoint(p); // !!! т.к. аналитическое тоже в паскалях, лучше перевести давление в атмосферы исключительно на выдачу
//        cout << curr_r << '\t' << num_value_P / 101325. << '\t' << analyt_value_P / 101325. << endl;
//        norm_err += (analyt_value_P - num_value_P) * (analyt_value_P - num_value_P);
//        norm_true += num_value_P * num_value_P;
//        curr_r += h_r;
//        i++;
//    }
//    curr_r = R;
//    double analyt_value_P = AnalitP(curr_r);
//    Coord3 p(curr_r, 0, 0.5);
//    double num_value_P = ResUInPoint(p); // т.к. аналитическое тоже в паскалях, лучше перевести давление в атмосферы исключительно на выдачу
//    cout << curr_r << '\t' << num_value_P / 101325. << '\t' << analyt_value_P / 101325. << endl;
//    norm_err += (analyt_value_P - num_value_P) * (analyt_value_P - num_value_P);
//    norm_true += num_value_P * num_value_P;
//
//    cout << endl << "Относительная норма вектора погрешности для аналитического решения:" << endl;
//    cout << sqrt(norm_err) / sqrt(norm_true) << endl << endl; // !!! тут проcто опечатка, нормировка не на то была
//}
//#pragma endregion

int main()
{
    setlocale(LC_ALL, "Russian");
    cout << "GenGrid" << endl;
    GenEndElGrid();
    cout << "ArrayParallelepipeds" << endl;
    ArrayParallelepipeds();
    cout << "GeneratePortrait" << endl;
    GeneratePortrait();
    cout << "GenFirstBoundCondit" << endl;
    GenFirstBoundCondit();
    cout << "GenSecBoundCondit" << endl;
    GenSecBoundCondit();
    cout << "BuildMatrA_VecB" << endl;
    BuildLocalMatrices();
    BuildMatrA_VecB();
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
    //cout << endl << "Суммарный небаланс на всех элементах: " << CalcSumNonBalance() << endl;
    BalancingFlows();

    /*CreateVecFict();
    Test();*/

    /*vector<double> vec_analit_P;
    VecAnalitP();*/

    //Output_result();
    //cout << endl << "Численное значение функции: " << ResUInPoint() << endl;
    //Output_u();
    //Output_x();
    //Output_y();
    //Output_z();
}


