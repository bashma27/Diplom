#include <iostream>
#include <fstream>
#include <vector>
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
vector<int> ia, ja, choice;
vector<double> aal, di, b, q, L_sq, di_sq, normal;
vector<double> k, k_ph, eta_ph; // массивы коэффициентов структурной проницаемости,
                                //         коэффициентов относительной фазовой проницаемости,
                                //         коэффициентов динамической вязкости соответственно
vector<Coord3> nodes; // узлы сетки, заданные координатами 
set<int> face_1; // множество с узлами первых краевых
vector<pair<int, vector<int>>> face_2_zp; // массив граней с краевыми 2 рода зон перфорации
vector<pair<int, vector<int>>> array_p; // массив параллелепипедов
vector<vector<double>> G, M; // локальные матрицы масс и жесткости
vector<int> ident_fict; // идентификаторы фиктивности (области или конечного элемента)
vector<int> i_ident_fict; 
vector<int> fict_nodes; // массив фиктивных узлов
#pragma endregion

#pragma region Функции краевых условий, их функций и вектора правой части f
double u_g(double x, double y, double z) { // краевое условие первого рода
    return x * y;
}

vector<double> grad_u(double x, double y) { // градиент функции u
    return { y, x };
}

double lambda(int num_sub) {
    return k[num_sub] * k_ph[num_sub] / eta_ph[num_sub];
}

double theta(int num_sub, int num, double x, double y) { // краевое условие второго рода
    switch (num)
    {
    case 0:
        return grad_u(x, y)[1] * 1 * lambda(num_sub);
        break;
    case 1:
        return grad_u(x, y)[1] * (-1) * lambda(num_sub);
        break;
    case 2:
        return grad_u(x, y)[0] * (1) * lambda(num_sub);
        break;
    case 3:
        return grad_u(x, y)[0] * (-1) * lambda(num_sub);
        break;
    }
}

double f(int num_sub, double x, double y) {
    return 0;
}
#pragma endregion

#pragma region Создание массива конечных элементов 
int GetNumSubarea(vector<int> end_el) { // получить номер подобласти
    for (int i = 0; i < W.size(); i++) {
        double min_x = nodes[end_el[0]].x; double min_y = nodes[end_el[0]].y; double min_z = nodes[end_el[0]].z;  //берём индексы координат прямоугольника
        double max_x = nodes[end_el[2]].x; double max_y = nodes[end_el[6]].y; double max_z = nodes[end_el[18]].z;
        if ( (min_x > W[i].second[0] || fabs(W[i].second[0] - min_x) < EPS) &&
             (W[i].second[1] > max_x || fabs(max_x - W[i].second[1]) < EPS) &&
             (min_y > W[i].second[2] || fabs(W[i].second[2] - min_y) < EPS) &&
             (W[i].second[3] > max_y || fabs(max_y - W[i].second[3]) < EPS) &&
             (min_z > W[i].second[4] || fabs(W[i].second[4] - min_z) < EPS) &&
             (W[i].second[5] > max_z || fabs(max_z - W[i].second[5]) < EPS) )
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

void AddLocalMatr(vector<int> node_num, vector<vector<double>> loc_matr) { // добавление матрицы в глобальную
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

void AddLocalVec(vector<int> node_num, vector<double> loc_vec) { // добавление вектора в глобальный
    for (int i = 0; i < 27; i++) {
        b[node_num[i]] += loc_vec[i];
    }
}

void AddLocalVecBound(int num_face_2_zp, vector<double> vec) { // добавление вектора из второго краевого в глобальный
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
        if ( (zone_perf_i[i_x].x == nodes[i].x || abs(zone_perf_i[i_x].x - nodes[i].x) < EPS) &&
             (zone_perf_i[i_y].y == nodes[i].y || abs(zone_perf_i[i_y].y - nodes[i].y) < EPS) ) return i;
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
        if (face_1.count(_i)) {
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
                if (face_1.count(i)) {
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
    for (int i = 0; i < NUM_NODES; i++) {
        if (IsFictitious(i)) {
            b[i] = 0;
            di[i] = double(1);
        }
    }
}

void ConsiderBoundConditSecType(int num_face_2_zp) { // учет краевых условий второго типа
    vector<double> b_s2(9, 0);
    vector<double> _theta(9);
    for (int i = 0; i < 9; i++) {
        _theta[i] = theta(face_2_zp[num_face_2_zp].second[9], face_2_zp[num_face_2_zp].first,
                           nodes[face_2_zp[num_face_2_zp].second[i]].x,
                           nodes[face_2_zp[num_face_2_zp].second[i]].y );
    }
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
            b_s2[i] += C[i][j] * _theta[j];
        }
    }
    AddLocalVecBound(num_face_2_zp, b_s2);
}

void ConsiderBoundCondit() { // учет всех краевых
    for (int i = 0; i < face_2_zp.size(); i++) { // учет вторых краевых
        ConsiderBoundConditSecType(i);
    }
    for (int n : face_1) { // учет первых краевых
        ConsiderBoundConditFirstType(n);
    }
    ConsiderFictitiousNodes();
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
                                   ( (h_y * h_z * G[i % 3][j % 3] * M[(i / 3) % 3][(j / 3) % 3] * M[i / 9][j / 9]) / (2700.0 * h_x) +
                                     (h_x * h_z * M[i % 3][j % 3] * G[(i / 3) % 3][(j / 3) % 3] * M[i / 9][j / 9]) / (2700.0 * h_y) +
                                     (h_x * h_y * M[i % 3][j % 3] * M[(i / 3) % 3][(j / 3) % 3] * G[i / 9][j / 9]) / (2700.0 * h_z) );
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
        if ( (point.x > min_x || fabs(point.x - min_x) < EPS) &&
             (max_x > point.x || fabs(max_x - point.x) < EPS) &&
             (point.y > min_y || fabs(point.y - min_y) < EPS) &&
             (max_y > point.y || fabs(max_y - point.y) < EPS) &&
             (point.z > min_z || fabs(point.z - min_z) < EPS) &&
             (max_z > point.z || fabs(max_z - point.z) < EPS) )
             return i;
    }
}
double ResUInPoint() {
    Coord3 point;
    cout << "Введите точку в которой необходимо узнать численное значение функции:" << endl;
    cout << "-----------------------------------" << endl;
    double min_x, max_x, min_y, max_y, min_z, max_z;
    min_x = W[0].second[0]; max_x = W[W.size() - 1].second[1];
    min_y = W[0].second[2]; max_y = W[W.size() - 1].second[3];
    min_z = W[0].second[4]; max_z = W[W.size() - 1].second[5];

    while (true) {
        cout << "(Возможный диапазон значений по x:)" << endl;        
        cout << "min:" << min_x << "    " << "max:" << max_x << endl;
        cin >> point.x;
        if ( (point.x > min_x || fabs(point.x - min_x) < EPS) &&
             (max_x > point.x || fabs(max_x - point.x) < EPS) ) {
             break;
        }
        else {
            cout << endl << "ВВЕДЕНОЕ ЧИСЛО НЕ ПОПАДАЕТ В ДИАПАЗОН ДОПУСТИМЫХ ЗНАЧЕНИЙ!!!" << endl;
        }            
    }

    while (true) {
        cout << endl << "(Возможный диапазон значений по y:)" << endl;
        cout << "min:" << min_y << "    " << "max:" << max_y << endl;
        cin >> point.y;
        if ( (point.y > min_y || fabs(point.y - min_y) < EPS) &&
             (max_y > point.y || fabs(max_y - point.y) < EPS) ) {
             break;
        }
        else {
            cout << endl << "ВВЕДЕНОЕ ЧИСЛО НЕ ПОПАДАЕТ В ДИАПАЗОН ДОПУСТИМЫХ ЗНАЧЕНИЙ!!!" << endl;
        }
    }

    while (true) {
        cout << endl << "(Возможный диапазон значений по z:)" << endl;
        cout << "min:" << min_z << "    " << "max:" << max_z << endl;
        cin >> point.z;
        if ( (point.z > min_z || fabs(point.z - min_z) < EPS) &&
             (max_z > point.z || fabs(max_z - point.z) < EPS) ) {           
             break;
        }
        else {
            cout << endl << "ВВЕДЕНОЕ ЧИСЛО НЕ ПОПАДАЕТ В ДИАПАЗОН ДОПУСТИМЫХ ЗНАЧЕНИЙ!!!" << endl;
        }
    }
    
    int num_end_el = GetNumEndEl(point);
    double res = 0;

    for (int i = 0; i < array_p[num_end_el].second.size(); i++) {
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

int main()
{
    setlocale(LC_ALL, "Russian");
    GenEndElGrid();
    ArrayParallelepipeds();
    GeneratePortrait();
    GenFirstBoundCondit();
    GenSecBoundCondit();
    BuildLocalMatrices();
    BuildMatrA_VecB();
    ConsiderBoundCondit();
    q.resize(NUM_NODES, 0);
    vector<double> r(NUM_NODES);
    vector<double> z(NUM_NODES);
    vector<double> Mult(NUM_NODES);
    vector<double> Az(NUM_NODES);
    int max_iter = 1000;
    double eps = 1e-15;
    MSG::LU_sq_MSG(q, r, z, Az, Mult, NUM_NODES, eps, max_iter);
    CreateVecFict();
    Test();
    Output_result();
    cout << endl << "Численное значение функции: " << ResUInPoint() << endl;
    //Output_u();
    //Output_x();
    //Output_y();
    //Output_z();
}


