#include "fem.h"
#include "get_funcs.h"

#pragma region Вспомогательные функции
double BasicFunc(int i, double input_point, double point_1, double point_2, double point_3) { // используемые базисные функции
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

double DerivativeLocalBasicFuncs(int i, double psi) { // производная по базисным функциям
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

void GenFictNodes() { // генерация фиктивных узлов
    for (int i = 0; i < NUM_NODES; i++) {
        if (IsFictitious(i)) {
            fict_nodes.push_back(i);
        }
    }
}

void GenArrayFictEndEl() { // генерация массива фиктивных конечных элементов
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

bool IsFictEl(int num_end_el) { // определение фиктивности конечного элемента
    for (int i = 0; i < fict_el.size(); i++) {
        if (fict_el[i] == num_end_el)
            return true;
    }
    return false;
}

bool IsFindFaceZP(int num_face) { // нахождение грани в массиве граней зон перфорации
    for (int i = 0; i < num_faces_zp.size(); i++) {
        if (num_faces_zp[i].first == num_face)
            return true;
    }
    return false;
}

bool IsFindFictFace(int num_face) { // нахождение грани в массиве фиктивных граней
    for (int i = 0; i < fict_faces.size(); i++) {
        if (fict_faces[i] == num_face)
            return true;
    }
    return false;
}
#pragma endregion

#pragma region Задание/вычисление краевых условий, потока и правой части
double u_g(double x, double y, double z) { // краевое условие первого рода
    return 130.;
}

double lambda(int num_sub) { // расчет лямбды как итогового коэфиициента
    double res = 0;
    for (int i = 0; i < num_ph; i++) {
        res += k_ph[num_sub][i] / eta_ph[i];
    }
    res *= K[num_sub];
    return res;
}

double theta(int num_face_2_zp) { // вычисление тетты по заданному потоку
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
    return set_flow_zp[num_zone_perf] * S_curr_face / (S_area_zp * S_area_zp);
}

double SpecFlow(int num_face_2_zp) { // заданный поток на нужной грани для зоны перфорации
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
    return set_flow_zp[num_zone_perf] * S_curr_face / S_area_zp;
}

double f(int num_sub, double x, double y) { // функция правой части диф. уравнения
    return 0;
}
#pragma endregion

#pragma region Генерация массива конечных элементов и портрета матрицы
void GenArrayParallelepipeds() { // генерация массива параллелепипедов
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
}

void GenPortraitMatr() { // генерация портрета матрицы
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
#pragma endregion

#pragma region Функции добавления в глоб. матрицу и вектор правой части
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

void AddLocalVec(vector<int>& node_num, vector<double>& loc_vec) { // добавление вектора в глобальный
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

#pragma region Работа с краевыми условиями и фиктивными узлами (генерация и учет)
void GenFirstBoundCondit() { // создание массива с первыми краевыми условиями
    int start_nodes = 0;
    int offset = 0;
    for (int i = 0; i < NUM_SPLIT_Z; i++) {
        for (int j = 0; j < NUM_SPLIT_X; j++) {
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
        }
        start_nodes = 0;
        offset += 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y;
    }
    start_nodes = NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y - NUM_NODES_IN_EDGE_X;
    offset = 0;
    for (int i = 0; i < NUM_SPLIT_Z; i++) {
        for (int j = 0; j < NUM_SPLIT_X; j++) {
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
        }
        start_nodes = NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y - NUM_NODES_IN_EDGE_X;
        offset += 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y;
    }
    start_nodes = 0;
    offset = 0;
    for (int i = 0; i < NUM_SPLIT_Z; i++) {
        for (int j = 0; j < NUM_SPLIT_Y; j++) {
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
        }
        start_nodes = 0;
        offset += 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y;
    }
    start_nodes = NUM_NODES_IN_EDGE_X - 1;
    offset = 0;
    for (int i = 0; i < NUM_SPLIT_Z; i++) {
        for (int j = 0; j < NUM_SPLIT_Y; j++) {
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
        }
        start_nodes = NUM_NODES_IN_EDGE_X - 1;
        offset += 2 * NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y;
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

vector<vector<pair<int, int>>> BuildColumnToRows() { // построить связь колонок к строкам
    // (в каких строках и по каким индексам матрицы встречается каждый столбец)
    vector<vector<pair<int, int>>> columnToRows(NUM_NODES);
    for (int i = 0; i < NUM_NODES; i++) {
        for (int j = ia[i]; j < ia[i + 1]; j++) {
            int col = ja[j];
            columnToRows[col].emplace_back(i, j);
        }
    }
    return columnToRows;
}

void ConsiderBoundConditFirstType(int n, vector<vector<pair<int, int>>>& columnToRows) { // учет краевых условий первого типа
    double x = nodes[n].x, y = nodes[n].y, z = nodes[n].z;
    b[n] = u_g(x, y, z);
    di[n] = 1.;

    for (int i = ia[n]; i < ia[n + 1]; i++) {
        int col = ja[i];
        if (face_1.count(col)) {
            aal[i] = 0.;
            continue;
        }
        b[col] -= b[n] * aal[i];
        aal[i] = 0.;
    }

    for (auto& [row, idx] : columnToRows[n]) {

        if (face_1.count(row)) continue;

        b[row] -= b[n] * aal[idx];
        aal[idx] = 0.;
    }
}

void ConsiderBoundConditSecType(int num_face_2_zp) { // учет краевых условий второго типа
    vector<double> b_s2(9, 0);
    double _theta = theta(num_face_2_zp);
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
            b_s2[i] += C[i][j] * _theta;
        }
    }
    AddLocalVecBound(num_face_2_zp, b_s2);
}

void ConsiderFictitiousNodes() { // учет фиктивных узлов
    for (int i = 0; i < NUM_NODES; i++) {
        if (di[i] == 0) {
            b[i] = 0;
            di[i] = double(1);
        }
    }
}

void ConsiderBoundCondit() { // учет всех краевых
    for (int i = 0; i < face_2_zp.size(); i++) {
        ConsiderBoundConditSecType(i);
    }
    cout << "second" << endl;
    auto columnToRows = BuildColumnToRows();
    for (int n : face_1) {
        ConsiderBoundConditFirstType(n, columnToRows);
    }
    cout << "first" << endl;
    ConsiderFictitiousNodes();
    cout << "fiction" << endl;
}
#pragma endregion

#pragma region Работа с локальными матрицами и векторами
void CreateLocalMatrStiff(int num_sub, vector<int>& node_num) { // локальная матрица жесткости
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

void CreateLocalVec(int num_sub, vector<int>& node_num) { // локальный вектор правой части
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
void BuildMatrVec() { // построить глобальную матрицу и вектор
    for (int i = 0; i < NUM_SPLIT_X * NUM_SPLIT_Y * NUM_SPLIT_Z; i++) {

        auto result = find(begin(ident_fict), end(ident_fict), array_p[i].second[0]);
        if (result != end(ident_fict)) i += i_ident_fict[result - begin(ident_fict)];

        vector<int> node_num(27);
        for (int j = 0; j < 27; j++) {
            node_num[j] = array_p[i].second[j];
        }
        int num_sub = array_p[i].first;
        CreateLocalMatrStiff(num_sub, node_num);
        CreateLocalVec(num_sub, node_num);
    }
}
#pragma endregion

#pragma region Тестирование для аналитически заданных функций
void Test() { // тест
    vector<double> q_u(NUM_NODES, 0);
    for (int i = 0; i < NUM_NODES; i++) {
        q_u[i] = u_g(nodes[i].x, nodes[i].y, nodes[i].z);
    }
    double norm_vec_err = 0, norm_vec_q_u = 0; // норма вектора погрешности и q_u
    for (int i = 0; i < NUM_NODES; i++) {
        if (IsFictitious(i)) continue;
        norm_vec_err += (q[i] - q_u[i]) * (q[i] - q_u[i]);
        norm_vec_q_u += (q_u[i]) * (q_u[i]);
    }
    cout << endl;
    cout << "Относительная норма вектора погрешности полученного решения:" << endl;
    cout << sqrt(norm_vec_err) / sqrt(norm_vec_q_u) << endl << endl;
} 
#pragma endregion

#pragma region Возврат значения функции в любой точке
double GetResUInPoint(Coord3 point) { // получить значение полученной функции в точке
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

#pragma region Вывод вектора результата решения СЛАУ в файл
void OutputResult() { // вывод результата в файл
    ofstream f_result;
    f_result.imbue(locale("Russian"));
    f_result.open("result.txt");
    for (int i = 0; i < NUM_NODES; i++) {
        f_result << setprecision(15) << q[i] << endl;
    }
    f_result.close();
}
#pragma endregion

#pragma region Сравнение с аналитическим решением
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
#pragma endregion