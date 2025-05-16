#include "flow_balancing.h"
#include "get_funcs.h"
#include "solver_slae.h"

#pragma region Расчет потоков
double CalcFlowF(vector<int>& node_num, int num_sub) { // расчет потока, создаваемого источником
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

double CalcFlowFaceYZ(vector<int>& node_num, int num_sub, int n) { // расчет потока через грань yz
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

double CalcFlowFaceXZ(vector<int>& node_num, int num_sub, int n) { // расчет потока через грань xz
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

double CalcFlowFaceXY(vector<int>& node_num, int num_sub, int n) { // расчет потока через грань xy
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
#pragma endregion

#pragma region Генерация граней со значением потока
void GenFacesFlowValue() { // генерация граней со значением потока

    // порядок заполнения: 
    // сначала грани по y, потом по x (слева - направо; от ближайшей до дальней (если смотреть в сечении xy то снизу вверх)), потом вверх по z
    // (потоки по верхней и нижней z известны и не нужна фактически, а внутренние потоки по z практически нулевые - можно не учитывать (скважина полностью проходит по глубине z))

    vector<int> node_num(27); int num_sub;
    vector<int> el_for_face; // конечный элемент, соответствующий грани
    bool left = true;
    bool face_y = true;
    int start_el = 0;

    for (int i = 0; i < NUM_SPLIT_Z; i++) {

        //самая передняя грань по y
        for (int j = 0; j < NUM_SPLIT_X; j++) {
            node_num = array_p[start_el + j].second;
            num_sub = array_p[start_el + j].first;
            faces_flow_value.push_back(-CalcFlowFaceXZ(node_num, num_sub, -1));
            el_for_face.push_back(start_el + j);
        }
        //самая передняя грань по y

        //самая левая грань по x
        node_num = array_p[start_el].second;
        num_sub = array_p[start_el].first;
        faces_flow_value.push_back(-CalcFlowFaceYZ(node_num, num_sub, -1));
        el_for_face.push_back(start_el);
        //самая левая грань по x

        //внутренние грани по x
        for (int j = 0; j < NUM_SPLIT_X - 1; j++) {
            double lambda_left_el, lambda_right_el;
            node_num = array_p[start_el + j].second;
            num_sub = array_p[start_el + j].first;
            double right_V_left_el = CalcFlowFaceYZ(node_num, num_sub, 1); // правый поток левого конечного элемента
            lambda_left_el = lambda(num_sub);

            node_num = array_p[start_el + j + 1].second;
            num_sub = array_p[start_el + j + 1].first;
            double left_V_right_el = CalcFlowFaceYZ(node_num, num_sub, -1); // левый поток правого конечного элемента
            lambda_right_el = lambda(num_sub);

            double res = (right_V_left_el * lambda_right_el - left_V_right_el * lambda_left_el) / (lambda_right_el + lambda_left_el);
            faces_flow_value.push_back(res);
            el_for_face.push_back(start_el + j);
        }
        //внутренние грани по x

        //самая правая грань по x
        node_num = array_p[start_el + NUM_SPLIT_X - 1].second;
        num_sub = array_p[start_el + NUM_SPLIT_X - 1].first;
        faces_flow_value.push_back(CalcFlowFaceYZ(node_num, num_sub, 1));
        el_for_face.push_back(start_el + NUM_SPLIT_X - 1);
        //самая правая грань по x

        start_el += NUM_SPLIT_X;
        int num_face_2_zp = 0;

        //основной цикл
        for (int j = 1; j < NUM_SPLIT_Y; j++) {

            //внутренние грани по y
            for (int k = 0; k < NUM_SPLIT_X; k++) {

                if (IsFictEl(start_el + k) && !IsFictEl(start_el + k - NUM_SPLIT_X)) {
                    node_num = array_p[start_el + k].second;
                    vector<int> _face_2_zp(9);
                    _face_2_zp[0] = node_num[0]; _face_2_zp[1] = node_num[1]; _face_2_zp[2] = node_num[2];
                    _face_2_zp[3] = node_num[9]; _face_2_zp[4] = node_num[10]; _face_2_zp[5] = node_num[11];
                    _face_2_zp[6] = node_num[18]; _face_2_zp[7] = node_num[19]; _face_2_zp[8] = node_num[20];
                    int num_face_2_zp = GetNumFace2ZP(_face_2_zp);
                    faces_flow_value.push_back(-SpecFlow(num_face_2_zp));
                    num_faces_zp.push_back({ faces_flow_value.size() - 1 , 0 });
                    num_face_2_zp++;
                }
                else if (!IsFictEl(start_el + k) && IsFictEl(start_el + k - NUM_SPLIT_X)) {
                    node_num = array_p[start_el + k].second;
                    vector<int> _face_2_zp(9);
                    _face_2_zp[0] = node_num[0]; _face_2_zp[1] = node_num[1]; _face_2_zp[2] = node_num[2];
                    _face_2_zp[3] = node_num[9]; _face_2_zp[4] = node_num[10]; _face_2_zp[5] = node_num[11];
                    _face_2_zp[6] = node_num[18]; _face_2_zp[7] = node_num[19]; _face_2_zp[8] = node_num[20];
                    int num_face_2_zp = GetNumFace2ZP(_face_2_zp);
                    faces_flow_value.push_back(SpecFlow(num_face_2_zp));
                    num_faces_zp.push_back({ faces_flow_value.size() - 1 , 1 });
                    num_face_2_zp++;
                }
                else if (IsFictEl(start_el + k) && IsFictEl(start_el + k - NUM_SPLIT_X)) {
                    faces_flow_value.push_back(0);
                    fict_faces.push_back(faces_flow_value.size() - 1);
                }
                else {
                    double lambda_front_el, lambda_back_el;
                    node_num = array_p[start_el + k - NUM_SPLIT_X].second;
                    num_sub = array_p[start_el + k - NUM_SPLIT_X].first;
                    double back_V_front_el = CalcFlowFaceXZ(node_num, num_sub, 1); // задний поток переднего конечного элемента
                    lambda_front_el = lambda(num_sub);

                    node_num = array_p[start_el + k].second;
                    num_sub = array_p[start_el + k].first;
                    double front_V_back_el = CalcFlowFaceXZ(node_num, num_sub, -1); // передний поток заднего конечного элемента
                    lambda_back_el = lambda(num_sub);

                    double res = (back_V_front_el * lambda_back_el - front_V_back_el * lambda_front_el) / (lambda_back_el + lambda_front_el);
                    faces_flow_value.push_back(res);
                    el_for_face.push_back(start_el + k);
                }
            }
            //внутренние грани по y

            //самая левая грань по x
            node_num = array_p[start_el].second;
            num_sub = array_p[start_el].first;
            faces_flow_value.push_back(-CalcFlowFaceYZ(node_num, num_sub, -1));
            el_for_face.push_back(start_el);
            //самая левая грань по x

            //внутренние грани по x
            for (int k = 0; k < NUM_SPLIT_X - 1; k++) {
                if (!IsFictEl(start_el + k) && IsFictEl(start_el + k + 1)) {
                    node_num = array_p[start_el + k].second;
                    vector<int> _face_2_zp(9);
                    _face_2_zp[0] = node_num[2]; _face_2_zp[1] = node_num[5]; _face_2_zp[2] = node_num[8];
                    _face_2_zp[3] = node_num[11]; _face_2_zp[4] = node_num[14]; _face_2_zp[5] = node_num[17];
                    _face_2_zp[6] = node_num[20]; _face_2_zp[7] = node_num[23]; _face_2_zp[8] = node_num[26];
                    int num_face_2_zp = GetNumFace2ZP(_face_2_zp);
                    faces_flow_value.push_back(-SpecFlow(num_face_2_zp));
                    num_faces_zp.push_back({ faces_flow_value.size() - 1 , 2 });
                    num_face_2_zp++;
                }
                else if (IsFictEl(start_el + k) && !IsFictEl(start_el + k + 1)) {
                    node_num = array_p[start_el + k].second;
                    vector<int> _face_2_zp(9);
                    _face_2_zp[0] = node_num[2]; _face_2_zp[1] = node_num[5]; _face_2_zp[2] = node_num[8];
                    _face_2_zp[3] = node_num[11]; _face_2_zp[4] = node_num[14]; _face_2_zp[5] = node_num[17];
                    _face_2_zp[6] = node_num[20]; _face_2_zp[7] = node_num[23]; _face_2_zp[8] = node_num[26];
                    int num_face_2_zp = GetNumFace2ZP(_face_2_zp);
                    faces_flow_value.push_back(SpecFlow(num_face_2_zp));
                    num_faces_zp.push_back({ faces_flow_value.size() - 1 , 3 });
                    num_face_2_zp++;
                }
                else if (IsFictEl(start_el + k) && IsFictEl(start_el + k + 1)) {
                    faces_flow_value.push_back(0);
                    fict_faces.push_back(faces_flow_value.size() - 1);
                }
                else {
                    double lambda_left_el, lambda_right_el;
                    node_num = array_p[start_el + k].second;
                    num_sub = array_p[start_el + k].first;
                    double right_V_left_el = CalcFlowFaceYZ(node_num, num_sub, 1); // правый поток левого конечного элемента
                    lambda_left_el = lambda(num_sub);

                    node_num = array_p[start_el + k + 1].second;
                    num_sub = array_p[start_el + k + 1].first;
                    double left_V_right_el = CalcFlowFaceYZ(node_num, num_sub, -1); // левый поток правого конечного элемента
                    lambda_right_el = lambda(num_sub);

                    double res = (right_V_left_el * lambda_right_el - left_V_right_el * lambda_left_el) / (lambda_right_el + lambda_left_el);
                    faces_flow_value.push_back(res);
                    el_for_face.push_back(start_el + k);
                }
            }
            //внутренние грани по x

            //самая правая грань по x
            node_num = array_p[start_el + NUM_SPLIT_X - 1].second;
            num_sub = array_p[start_el + NUM_SPLIT_X - 1].first;
            faces_flow_value.push_back(CalcFlowFaceYZ(node_num, num_sub, 1));
            el_for_face.push_back(start_el + NUM_SPLIT_X - 1);
            //самая правая грань по x

            if (j == NUM_SPLIT_Y - 1) continue;
            start_el += NUM_SPLIT_X;
        }
        //основной цикл

        //самая дальняя грань по y
        for (int j = 0; j < NUM_SPLIT_X; j++) {
            node_num = array_p[start_el + j].second;
            num_sub = array_p[start_el + j].first;
            faces_flow_value.push_back(CalcFlowFaceXZ(node_num, num_sub, 1));
            el_for_face.push_back(start_el + j);
        }
        //самая дальняя грань по y

        start_el += NUM_SPLIT_X;
    }
}
#pragma endregion

#pragma region Отчистить старое и создать новый портрет для матрицы
void ClearAndGenPortMatr() { // отчистить старое и создать новый портрет для матрицы
    ia.clear();
    ja.clear();
    aal.clear();
    di.clear();
    b.clear();
    L_sq.clear();
    di_sq.clear();
    list_face.resize(faces_flow_value.size());
    bool face_y = true;
    int num_face_layer = faces_flow_value.size() / NUM_SPLIT_Z; // кол-во ребер в слое (в плоскости xy в одном разбиении по z)
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


    ia.resize(faces_flow_value.size() + 1, 0);
    for (int i = 2; i <= faces_flow_value.size(); i++) {
        ia[i] = ia[i - 1] + list_face[i - 1].size();
    }
    ja.resize(ia[faces_flow_value.size()]);
    auto iter = ja.begin();
    for (int i = 0; i < faces_flow_value.size(); i++) {
        copy(list_face[i].begin(), list_face[i].end(), iter);
        iter += list_face[i].size();
    }
    di.resize(faces_flow_value.size());
    b.resize(faces_flow_value.size());
    L_sq.resize(ia[faces_flow_value.size()]);
    di_sq.resize(faces_flow_value.size());
}
#pragma endregion

#pragma region Генерация матрицы и вектора
void GenMatrAndVec() { // генерация матрицы и вектора

    auto max_it = max_element(faces_flow_value.begin(), faces_flow_value.end(),
        [](double a, double b) { return fabs(a) < fabs(b); });
    double max_V = fabs(*max_it);

    double beta = 1e+8;
    double alpha;
    int num_face_layer = faces_flow_value.size() / NUM_SPLIT_Z; // кол-во ребер в слое (в плоскости xy в одном разбиении по z)

    for (int k = 0; k < NUM_SPLIT_Z; k++) {

        // начальные компоненты вектора
        for (int i = 0; i < NUM_SPLIT_X; i++) { // не работает для разбиений по z
            int S_g_1, S_g_2, S_g_3, S_g_4;

            if (fabs(faces_flow_value[k * num_face_layer + i]) > (max_V * 1e-4))
                alpha = 1. / max_V;
            else
                alpha = 1. / (max_V * 1e-4);

            S_g_1 = (faces_flow_value[k * num_face_layer + i] < 0) ? 1 : -1;
            S_g_2 = (faces_flow_value[k * num_face_layer + i + NUM_SPLIT_X] < 0) ? 1 : -1;
            S_g_3 = (faces_flow_value[k * num_face_layer + i + NUM_SPLIT_X + 1] < 0) ? -1 : 1;
            S_g_4 = (faces_flow_value[k * num_face_layer + i + 2 * NUM_SPLIT_X + 1] < 0) ? -1 : 1;
            di[k * num_face_layer + i] = beta + alpha;
            b[k * num_face_layer + i] -= beta * S_g_1 * (S_g_1 * fabs(faces_flow_value[k * num_face_layer + i]) + S_g_2 * fabs(faces_flow_value[k * num_face_layer + i + NUM_SPLIT_X]) +
                S_g_3 * fabs(faces_flow_value[k * num_face_layer + i + NUM_SPLIT_X + 1]) + S_g_4 * fabs(faces_flow_value[k * num_face_layer + i + 2 * NUM_SPLIT_X + 1]));
        }

        for (int i = 0; i < num_face_layer; i++) {
            // расчет диагонали
            if (fabs(faces_flow_value[k * num_face_layer + i]) > (max_V * 1e-4))
                alpha = 1. / max_V;
            else
                alpha = 1. / (max_V * 1e-4);

            int S_g_1, S_g_2;

            // расчет компонент матрицы
            for (int j = 0; j < list_face[k * num_face_layer + i].size(); j++) {

                if (list_face[k * num_face_layer + i].size() == 1 || (k * num_face_layer + i) == list_face[k * num_face_layer + i][list_face[k * num_face_layer + i].size() - 1] + 1) { // грани по x
                    if (list_face[k * num_face_layer + i].size() == 1) { // самый левый x
                        di[k * num_face_layer + i] = beta + alpha;
                        S_g_1 = (faces_flow_value[k * num_face_layer + i] < 0) ? 1 : -1;
                        S_g_2 = (faces_flow_value[list_face[k * num_face_layer + i][j]] < 0) ? 1 : -1;

                        // расчет компонент вектора правой части
                        int dif_S_g_1, dif_S_g_2;
                        dif_S_g_1 = (faces_flow_value[k * num_face_layer + i + 1] < 0) ? -1 : 1;
                        dif_S_g_2 = (faces_flow_value[list_face[k * num_face_layer + i][j] + 2 * NUM_SPLIT_X + 1] < 0) ? -1 : 1;

                        b[k * num_face_layer + i] -= beta * S_g_1 * (S_g_1 * fabs(faces_flow_value[k * num_face_layer + i]) + S_g_2 * fabs(faces_flow_value[list_face[k * num_face_layer + i][j]]) +
                            dif_S_g_1 * fabs(faces_flow_value[k * num_face_layer + i + 1]) + dif_S_g_2 * fabs(faces_flow_value[list_face[k * num_face_layer + i][j] + 2 * NUM_SPLIT_X + 1]));
                    }
                    else if (list_face[k * num_face_layer + i].size() == 2) { // самый правый x
                        di[k * num_face_layer + i] = beta + alpha;
                        S_g_1 = (faces_flow_value[k * num_face_layer + i] < 0) ? -1 : 1;
                        S_g_2 = (faces_flow_value[list_face[k * num_face_layer + i][j]] < 0) ? 1 : -1;

                        // расчет компонент вектора правой части
                        if (j == 0) {
                            int dif_S_g_1, dif_S_g_2;
                            dif_S_g_1 = (faces_flow_value[k * num_face_layer + i - 1] < 0) ? 1 : -1;
                            dif_S_g_2 = (faces_flow_value[list_face[k * num_face_layer + i][j] + 2 * NUM_SPLIT_X + 1] < 0) ? -1 : 1;

                            b[k * num_face_layer + i] -= beta * S_g_1 * (S_g_1 * fabs(faces_flow_value[k * num_face_layer + i]) + S_g_2 * fabs(faces_flow_value[list_face[k * num_face_layer + i][j]]) +
                                dif_S_g_1 * fabs(faces_flow_value[k * num_face_layer + i - 1]) + dif_S_g_2 * fabs(faces_flow_value[list_face[k * num_face_layer + i][j] + 2 * NUM_SPLIT_X + 1]));
                        }
                    }
                    else {
                        di[k * num_face_layer + i] = 2 * beta + alpha;
                        if (j == 0 || j == 2) {
                            S_g_1 = (faces_flow_value[k * num_face_layer + i] < 0) ? -1 : 1;
                            S_g_2 = (faces_flow_value[list_face[k * num_face_layer + i][j]] < 0) ? 1 : -1;

                            // расчет компонент вектора правой части
                            if (j == 0) {
                                int dif_S_g_1, dif_S_g_2;
                                dif_S_g_1 = (faces_flow_value[k * num_face_layer + i - 1] < 0) ? 1 : -1;
                                dif_S_g_2 = (faces_flow_value[list_face[k * num_face_layer + i][j] + 2 * NUM_SPLIT_X + 1] < 0) ? -1 : 1;

                                b[k * num_face_layer + i] -= beta * S_g_1 * (S_g_1 * fabs(faces_flow_value[k * num_face_layer + i]) + S_g_2 * fabs(faces_flow_value[list_face[k * num_face_layer + i][j]]) +
                                    dif_S_g_1 * fabs(faces_flow_value[k * num_face_layer + i - 1]) + dif_S_g_2 * fabs(faces_flow_value[list_face[k * num_face_layer + i][j] + 2 * NUM_SPLIT_X + 1]));
                            }
                        }
                        else {
                            S_g_1 = (faces_flow_value[k * num_face_layer + i] < 0) ? 1 : -1;
                            S_g_2 = (faces_flow_value[list_face[k * num_face_layer + i][j]] < 0) ? 1 : -1;

                            // расчет компонент вектора правой части
                            int dif_S_g_1, dif_S_g_2;
                            dif_S_g_1 = (faces_flow_value[k * num_face_layer + i + 1] < 0) ? -1 : 1;
                            dif_S_g_2 = (faces_flow_value[list_face[k * num_face_layer + i][j] + 2 * NUM_SPLIT_X + 1] < 0) ? -1 : 1;

                            b[k * num_face_layer + i] -= beta * S_g_1 * (S_g_1 * fabs(faces_flow_value[k * num_face_layer + i]) + S_g_2 * fabs(faces_flow_value[list_face[k * num_face_layer + i][j]]) +
                                dif_S_g_1 * fabs(faces_flow_value[k * num_face_layer + i + 1]) + dif_S_g_2 * fabs(faces_flow_value[list_face[k * num_face_layer + i][j] + 2 * NUM_SPLIT_X + 1]));
                        }
                    }
                }
                else { // грани по y  
                    S_g_1 = (faces_flow_value[k * num_face_layer + i] < 0) ? -1 : 1;
                    if (j == 0 || j == 1) {
                        S_g_2 = (faces_flow_value[list_face[k * num_face_layer + i][j]] < 0) ? 1 : -1;
                    }
                    else {
                        S_g_2 = (faces_flow_value[list_face[k * num_face_layer + i][j]] < 0) ? -1 : 1;
                    }

                    // расчет компонент вектора правой части для граней в конце
                    if (i >= (2 * NUM_SPLIT_X + 1) * NUM_SPLIT_Y && j == 0) {
                        di[k * num_face_layer + i] = beta + alpha;
                        int dif_S_g_1, dif_S_g_2;
                        dif_S_g_1 = (faces_flow_value[list_face[k * num_face_layer + i][j] + NUM_SPLIT_X] < 0) ? 1 : -1;
                        dif_S_g_2 = (faces_flow_value[list_face[k * num_face_layer + i][j] + NUM_SPLIT_X + 1] < 0) ? -1 : 1;

                        b[k * num_face_layer + i] -= beta * S_g_1 * (S_g_1 * fabs(faces_flow_value[k * num_face_layer + i]) + S_g_2 * fabs(faces_flow_value[list_face[k * num_face_layer + i][j]]) +
                            dif_S_g_1 * fabs(faces_flow_value[list_face[k * num_face_layer + i][j] + NUM_SPLIT_X]) + dif_S_g_2 * fabs(faces_flow_value[list_face[k * num_face_layer + i][j] + NUM_SPLIT_X + 1]));
                    }
                    else if (j == 0) { // не в конце
                        di[k * num_face_layer + i] = 2 * beta + alpha;
                        int dif_S_g_1, dif_S_g_2;
                        dif_S_g_1 = (faces_flow_value[list_face[k * num_face_layer + i][j] + NUM_SPLIT_X] < 0) ? 1 : -1;
                        dif_S_g_2 = (faces_flow_value[list_face[k * num_face_layer + i][j] + NUM_SPLIT_X + 1] < 0) ? -1 : 1;

                        b[k * num_face_layer + i] -= beta * S_g_1 * (S_g_1 * fabs(faces_flow_value[k * num_face_layer + i]) + S_g_2 * fabs(faces_flow_value[list_face[k * num_face_layer + i][j]]) +
                            dif_S_g_1 * fabs(faces_flow_value[list_face[k * num_face_layer + i][j] + NUM_SPLIT_X]) + dif_S_g_2 * fabs(faces_flow_value[list_face[k * num_face_layer + i][j] + NUM_SPLIT_X + 1]));

                        int d_S_g_1 = (faces_flow_value[k * num_face_layer + i] < 0) ? 1 : -1;
                        int d_S_g_2 = (faces_flow_value[k * num_face_layer + i + 2 * NUM_SPLIT_X + 1] < 0) ? -1 : 1;
                        dif_S_g_1 = (faces_flow_value[k * num_face_layer + i + NUM_SPLIT_X] < 0) ? 1 : -1;
                        dif_S_g_2 = (faces_flow_value[k * num_face_layer + i + NUM_SPLIT_X + 1] < 0) ? -1 : 1;

                        b[k * num_face_layer + i] -= beta * d_S_g_1 * (d_S_g_1 * fabs(faces_flow_value[k * num_face_layer + i]) + d_S_g_2 * fabs(faces_flow_value[k * num_face_layer + i + 2 * NUM_SPLIT_X + 1]) +
                            dif_S_g_1 * fabs(faces_flow_value[k * num_face_layer + i + NUM_SPLIT_X]) + dif_S_g_2 * fabs(faces_flow_value[k * num_face_layer + i + NUM_SPLIT_X + 1]));
                    }
                }
                double a = beta * S_g_1 * S_g_2;
                aal.push_back(beta * S_g_1 * S_g_2);
            }
            if (IsFindFaceZP(i)) {
                di[k * num_face_layer + i] = beta + alpha;
            }
        }
    }
}
#pragma endregion

#pragma region Учет потоков
void ConsiderKnownFlows(int n) { // учет известных потоков
    b[n] = 0;
    di[n] = 1;
    for (int i = ia[n]; i < ia[n + 1]; i++) {
        int _i = ja[i];
        if (IsFindFaceZP(_i)) {
            aal[i] = 0;
            continue;
        }
        b[_i] -= b[n] * aal[i];
        aal[i] = 0;
    }
    for (int i = n; i < faces_flow_value.size(); i++) {
        int k = 0;
        for (int j = ia[i]; j < ia[i + 1]; j++) {
            if (ja[j] == n) {
                if (IsFindFaceZP(i)) {
                    aal[j] = 0;
                    continue;
                }
                b[i] -= b[n] * aal[j];
                aal[j] = 0;
            }
        }
    }
}

void ConsiderFictFlows(int n) { // учет фиктивных потоков
    b[n] = 0;
    di[n] = 1;
    for (int i = ia[n]; i < ia[n + 1]; i++) {
        int _i = ja[i];
        if (IsFindFictFace(_i)) {
            aal[i] = 0;
            continue;
        }
        b[_i] -= b[n] * aal[i];
        aal[i] = 0;
    }
    for (int i = n; i < faces_flow_value.size(); i++) {
        int k = 0;
        for (int j = ia[i]; j < ia[i + 1]; j++) {
            if (ja[j] == n) {
                if (IsFindFictFace(i)) {
                    aal[j] = 0;
                    continue;
                }
                b[i] -= b[n] * aal[j];
                aal[j] = 0;
            }
        }
    }
}

void ConsiderFlows() { // учет всех краевых
    for (int i = 0; i < num_faces_zp.size(); i++) {
        ConsiderKnownFlows(num_faces_zp[i].first);
    }
    for (int i = 0; i < fict_faces.size(); i++) {
        ConsiderFictFlows(fict_faces[i]);
    }
}
#pragma endregion

#pragma region Расчет небаланса и нахождение сбалансированных потоков
double CalcSumNonBalance() { // рассчитать суммарный небаланс
    double res = 0;
    int start_face = 0;
    int num_face_layer = faces_flow_value.size() / NUM_SPLIT_Z; // кол-во ребер в слое (в плоскости xy в одном разбиении по z)
    for (int i = 0; i < NUM_SPLIT_Z; i++) {
        for (int j = 0; j < NUM_SPLIT_Y; j++) {
            for (int k = 0; k < NUM_SPLIT_X; k++) {
                int curr_el = i * NUM_SPLIT_X * NUM_SPLIT_Y + j * NUM_SPLIT_Y + k;

                if (IsFictEl(curr_el))
                    continue;

                int S_g_1, S_g_2, S_g_3, S_g_4;
                S_g_1 = (faces_flow_value[start_face + k] < 0) ? 1 : -1;
                S_g_2 = (faces_flow_value[start_face + k + NUM_SPLIT_X] < 0) ? 1 : -1;
                S_g_3 = (faces_flow_value[start_face + k + NUM_SPLIT_X + 1] < 0) ? -1 : 1;
                S_g_4 = (faces_flow_value[start_face + k + 2 * NUM_SPLIT_X + 1] < 0) ? -1 : 1;
                double a = S_g_1 * fabs(faces_flow_value[start_face + k]) + S_g_2 * fabs(faces_flow_value[start_face + k + NUM_SPLIT_X]) +
                    S_g_3 * fabs(faces_flow_value[start_face + k + NUM_SPLIT_X + 1]) + S_g_4 * fabs(faces_flow_value[start_face + k + 2 * NUM_SPLIT_X + 1]);
                res += a;
            }
            start_face += 2 * NUM_SPLIT_X + 1;
        }
        start_face += NUM_SPLIT_X;
    }
    return res;
}

void FindBalancedFlows() { // найти сбалансированные потоки
    for (int i = 0; i < faces_flow_value.size(); i++) {
        int S_g = (faces_flow_value[i] < 0) ? -1 : 1;
        faces_flow_value[i] = S_g * fabs(faces_flow_value[i]) + S_g * q_flow[i];
    }
}
#pragma endregion

#pragma region Балансировка потоков
void BalancingFlows() { // балансировка потоков
    GenFacesFlowValue();
    cout << endl << "Стартовый небаланс на всех элементах: " << CalcSumNonBalance() << endl;
    ClearAndGenPortMatr();
    GenMatrAndVec();
    ConsiderFlows();
    q_flow.resize(faces_flow_value.size(), 0);
    vector<double> r(faces_flow_value.size());
    vector<double> z(faces_flow_value.size());
    vector<double> Mult(faces_flow_value.size());
    vector<double> Az(faces_flow_value.size());
    int max_iter = 1000;
    double eps = 1e-15;
    cout << endl;
    MSG::LU_sq_MSG(q_flow, r, z, Az, Mult, faces_flow_value.size(), eps, max_iter);
    FindBalancedFlows();
    cout << endl << "Небаланс на всех элементах после балансировки: " << CalcSumNonBalance() << endl;
}
#pragma endregion