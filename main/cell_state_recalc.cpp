#include "cell_state_recalc.h"
#include "flow_balancing.h"

#pragma region Расчет потоков для каждой фазы
void CalcFlowPh() { // расчитать потоки для каждоый фазы в отдельности
    // порядок заполнения: 
    // сначала грани по y, потом по x (слева - направо; от ближайшей до дальней (если смотреть в сечении xy то снизу вверх)), потом вверх по z
    // (потоки по верхней и нижней z известны и не нужна фактически, а внутренние потоки по z практически нулевые - можно не учитывать (скважина полностью проходит по глубине z))
    faces_flow_ph_value.assign(faces_flow_value.size(), vector<double>(num_ph));
    vector<int> node_num(27);
    int num_sub;
    int num_pos_flow = 0;
    int start_el = 0;
    int num_flow = 0;
    double sum;
    for (int i = 0; i < NUM_SPLIT_Z; i++) {

        //самая передняя грань по y
        for (int j = 0; j < NUM_SPLIT_X; j++) {
            num_sub = array_p[start_el + j].first;
            sum = 0.;
            for (int ph = 0; ph < num_ph; ph++) {
                sum += k_ph(start_el + j, ph) / eta_ph[ph];
            }
            for (int ph = 0; ph < num_ph; ph++) {
                faces_flow_ph_value[num_flow][ph] = faces_flow_value[num_flow] * k_ph(start_el + j, ph) / (eta_ph[ph] * sum);
            }
            num_flow++;
        }
        //самая передняя грань по y

        //самая левая грань по x
        num_sub = array_p[start_el].first;
        sum = 0.;
        for (int ph = 0; ph < num_ph; ph++) {
            sum += k_ph(start_el, ph) / eta_ph[ph];
        }
        for (int ph = 0; ph < num_ph; ph++) {
            faces_flow_ph_value[num_flow][ph] = faces_flow_value[num_flow] * k_ph(start_el, ph) / (eta_ph[ph] * sum);
        }
        num_flow++;
        //самая левая грань по x

        //внутренние грани по x
        for (int j = 0; j < NUM_SPLIT_X - 1; j++) {
            int S_g = (faces_flow_value[num_flow] < 0) ? -1 : 1;
            int curr_el;
            if (S_g < 0)
                curr_el = start_el + j + 1;                           
            else
                curr_el = start_el + j;
            num_sub = array_p[curr_el].first;
            sum = 0.;
            for (int ph = 0; ph < num_ph; ph++) {
                sum += k_ph(curr_el, ph) / eta_ph[ph];
            }
            for (int ph = 0; ph < num_ph; ph++) {
                faces_flow_ph_value[num_flow][ph] = faces_flow_value[num_flow] * k_ph(curr_el, ph) / (eta_ph[ph] * sum);
            }
            num_flow++;
        }
        //внутренние грани по x

        //самая правая грань по x
        num_sub = array_p[start_el + NUM_SPLIT_X - 1].first;
        sum = 0.;
        for (int ph = 0; ph < num_ph; ph++) {
            sum += k_ph(start_el + NUM_SPLIT_X - 1, ph) / eta_ph[ph];
        }
        for (int ph = 0; ph < num_ph; ph++) {
            faces_flow_ph_value[num_flow][ph] = faces_flow_value[num_flow] * k_ph(start_el + NUM_SPLIT_X - 1, ph) / (eta_ph[ph] * sum);
        }
        num_flow++;
        //самая правая грань по x

        start_el += NUM_SPLIT_X;

        //основной цикл
        for (int j = 1; j < NUM_SPLIT_Y; j++) {

            //внутренние грани по y
            for (int k = 0; k < NUM_SPLIT_X; k++) {

                if (IsFictEl(start_el + k) && !IsFictEl(start_el + k - NUM_SPLIT_X)) {
                    if (faces_flow_value[num_flow] < -EPS) {
                        node_num = array_p[start_el + k].second;
                        vector<int> _face_2_zp(9);
                        _face_2_zp[0] = node_num[0]; _face_2_zp[1] = node_num[1]; _face_2_zp[2] = node_num[2];
                        _face_2_zp[3] = node_num[9]; _face_2_zp[4] = node_num[10]; _face_2_zp[5] = node_num[11];
                        _face_2_zp[6] = node_num[18]; _face_2_zp[7] = node_num[19]; _face_2_zp[8] = node_num[20];
                        int num_face_2_zp = GetNumFace2ZP(_face_2_zp);
                        int num_pos_flow = GetNumPosFlowFromAllSetFlows(face_2_zp[num_face_2_zp].second[10]);
                        for (int ph = 0; ph < num_ph; ph++) {
                            faces_flow_ph_value[num_flow][ph] = faces_flow_value[num_flow] * satur_pos_set_flow[num_pos_flow][ph];
                        }                   
                        num_flow++;
                        continue;
                    }
                    num_sub = array_p[start_el + k - NUM_SPLIT_X].first;
                    sum = 0.;
                    for (int ph = 0; ph < num_ph; ph++) {
                        sum += k_ph(start_el + k - NUM_SPLIT_X, ph) / eta_ph[ph];
                    }
                    for (int ph = 0; ph < num_ph; ph++) {
                        faces_flow_ph_value[num_flow][ph] = faces_flow_value[num_flow] * k_ph(start_el + k - NUM_SPLIT_X, ph) / (eta_ph[ph] * sum);
                    }
                    num_flow++;
                }
                else if (!IsFictEl(start_el + k) && IsFictEl(start_el + k - NUM_SPLIT_X)) {
                    if (faces_flow_value[num_flow] > EPS) {
                        node_num = array_p[start_el + k].second;
                        vector<int> _face_2_zp(9);
                        _face_2_zp[0] = node_num[0]; _face_2_zp[1] = node_num[1]; _face_2_zp[2] = node_num[2];
                        _face_2_zp[3] = node_num[9]; _face_2_zp[4] = node_num[10]; _face_2_zp[5] = node_num[11];
                        _face_2_zp[6] = node_num[18]; _face_2_zp[7] = node_num[19]; _face_2_zp[8] = node_num[20];
                        int num_face_2_zp = GetNumFace2ZP(_face_2_zp);
                        int num_pos_flow = GetNumPosFlowFromAllSetFlows(face_2_zp[num_face_2_zp].second[10]);
                        for (int ph = 0; ph < num_ph; ph++) {
                            faces_flow_ph_value[num_flow][ph] = faces_flow_value[num_flow] * satur_pos_set_flow[num_pos_flow][ph];
                        }
                        num_flow++;
                        continue;
                    }
                    num_sub = array_p[start_el + k].first;
                    sum = 0.;
                    for (int ph = 0; ph < num_ph; ph++) {
                        sum += k_ph(start_el + k, ph);
                    }
                    for (int ph = 0; ph < num_ph; ph++) {
                        faces_flow_ph_value[num_flow][ph] = faces_flow_value[num_flow] * k_ph(start_el + k, ph) / (eta_ph[ph] * sum);
                    }
                    num_flow++;
                }
                else if (IsFictEl(start_el + k) && IsFictEl(start_el + k - NUM_SPLIT_X)) {
                    num_flow++;
                }
                else {
                    int S_g = (faces_flow_value[num_flow] < 0) ? -1 : 1;
                    int curr_el;
                    if (S_g < 0)
                        curr_el = start_el + k;
                    else
                        curr_el = start_el + k - NUM_SPLIT_X;
                    sum = 0.;
                    for (int ph = 0; ph < num_ph; ph++) {
                        sum += k_ph(curr_el, ph) / eta_ph[ph];
                    }
                    for (int ph = 0; ph < num_ph; ph++) {
                        faces_flow_ph_value[num_flow][ph] = faces_flow_value[num_flow] * k_ph(curr_el, ph) / (eta_ph[ph] * sum);
                    }
                    num_flow++;
                }
            }
            //внутренние грани по y

            //самая левая грань по x
            num_sub = array_p[start_el].first;
            sum = 0.;
            for (int ph = 0; ph < num_ph; ph++) {
                sum += k_ph(start_el, ph) / eta_ph[ph];
            }
            for (int ph = 0; ph < num_ph; ph++) {
                faces_flow_ph_value[num_flow][ph] = faces_flow_value[num_flow] * k_ph(start_el, ph) / (eta_ph[ph] * sum);
            }
            num_flow++;
            //самая левая грань по x

            //внутренние грани по x
            for (int k = 0; k < NUM_SPLIT_X - 1; k++) {
                if (!IsFictEl(start_el + k) && IsFictEl(start_el + k + 1)) {
                    if (faces_flow_value[num_flow] < -EPS) {
                        node_num = array_p[start_el + k].second;
                        vector<int> _face_2_zp(9);
                        _face_2_zp[0] = node_num[2]; _face_2_zp[1] = node_num[5]; _face_2_zp[2] = node_num[8];
                        _face_2_zp[3] = node_num[11]; _face_2_zp[4] = node_num[14]; _face_2_zp[5] = node_num[17];
                        _face_2_zp[6] = node_num[20]; _face_2_zp[7] = node_num[23]; _face_2_zp[8] = node_num[26];
                        int num_face_2_zp = GetNumFace2ZP(_face_2_zp);
                        int num_pos_flow = GetNumPosFlowFromAllSetFlows(face_2_zp[num_face_2_zp].second[10]);
                        for (int ph = 0; ph < num_ph; ph++) {
                            faces_flow_ph_value[num_flow][ph] = faces_flow_value[num_flow] * satur_pos_set_flow[num_pos_flow][ph];
                        }
                        num_flow++;
                        continue;
                    }
                    num_sub = array_p[start_el + k].first;
                    sum = 0.;
                    for (int ph = 0; ph < num_ph; ph++) {
                        sum += k_ph(start_el + k, ph) / eta_ph[ph];
                    }
                    for (int ph = 0; ph < num_ph; ph++) {
                        faces_flow_ph_value[num_flow][ph] = faces_flow_value[num_flow] * k_ph(start_el + k, ph) / (eta_ph[ph] * sum);
                    }
                    num_flow++;
                }
                else if (IsFictEl(start_el + k) && !IsFictEl(start_el + k + 1)) {
                    if (faces_flow_value[num_flow] > EPS) {
                        node_num = array_p[start_el + k].second;
                        vector<int> _face_2_zp(9);
                        _face_2_zp[0] = node_num[2]; _face_2_zp[1] = node_num[5]; _face_2_zp[2] = node_num[8];
                        _face_2_zp[3] = node_num[11]; _face_2_zp[4] = node_num[14]; _face_2_zp[5] = node_num[17];
                        _face_2_zp[6] = node_num[20]; _face_2_zp[7] = node_num[23]; _face_2_zp[8] = node_num[26];
                        int num_face_2_zp = GetNumFace2ZP(_face_2_zp);
                        int num_pos_flow = GetNumPosFlowFromAllSetFlows(face_2_zp[num_face_2_zp].second[10]);
                        for (int ph = 0; ph < num_ph; ph++) {
                            faces_flow_ph_value[num_flow][ph] = faces_flow_value[num_flow] * satur_pos_set_flow[num_pos_flow][ph];
                        }
                        num_flow++;
                        continue;
                    }
                    num_sub = array_p[start_el + k + 1].first;
                    sum = 0.;
                    for (int ph = 0; ph < num_ph; ph++) {
                        sum += k_ph(start_el + k + 1, ph) / eta_ph[ph];
                    }
                    for (int ph = 0; ph < num_ph; ph++) {
                        faces_flow_ph_value[num_flow][ph] = faces_flow_value[num_flow] * k_ph(start_el + k + 1, ph) / (eta_ph[ph] * sum);
                    }
                    num_flow++;
                }
                else if (IsFictEl(start_el + k) && IsFictEl(start_el + k + 1)) {
                    num_flow++;
                }
                else {
                    int S_g = (faces_flow_value[num_flow] < 0) ? -1 : 1;
                    int curr_el;
                    if (S_g < 0)
                        curr_el = start_el + k + 1;
                    else
                        curr_el = start_el + k;
                    sum = 0.;
                    for (int ph = 0; ph < num_ph; ph++) {
                        sum += k_ph(curr_el, ph) / eta_ph[ph];
                    }
                    for (int ph = 0; ph < num_ph; ph++) {
                        faces_flow_ph_value[num_flow][ph] = faces_flow_value[num_flow] * k_ph(curr_el, ph) / (eta_ph[ph] * sum);
                    }
                    num_flow++;
                }
            }
            //внутренние грани по x

            //самая правая грань по x
            num_sub = array_p[start_el + NUM_SPLIT_X - 1].first;
            sum = 0.;
            for (int ph = 0; ph < num_ph; ph++) {
                sum += k_ph(start_el + NUM_SPLIT_X - 1, ph) / eta_ph[ph];
            }
            for (int ph = 0; ph < num_ph; ph++) {
                faces_flow_ph_value[num_flow][ph] = faces_flow_value[num_flow] * k_ph(start_el + NUM_SPLIT_X - 1, ph) / (eta_ph[ph] * sum);
            }
            num_flow++;
            //самая правая грань по x

            if (j == NUM_SPLIT_Y - 1) continue;
            start_el += NUM_SPLIT_X;
        }
        //основной цикл

        //самая дальняя грань по y
        for (int j = 0; j < NUM_SPLIT_X; j++) {
            num_sub = array_p[start_el + j].first;
            sum = 0.;
            for (int ph = 0; ph < num_ph; ph++) {
                sum += k_ph(start_el + j, ph) / eta_ph[ph];
            }
            for (int ph = 0; ph < num_ph; ph++) {
                faces_flow_ph_value[num_flow][ph] = faces_flow_value[num_flow] * k_ph(start_el + j, ph) / (eta_ph[ph] * sum);
            }
            num_flow++;
        }
        //самая дальняя грань по y

        start_el += NUM_SPLIT_X;
    }
}

void RecalcFlowPh(double min_t) { // регулировка потоков для минимального шага по времени
    for (int ph = 0; ph < num_ph; ph++) {
        int start_face = 0;
        int curr_el = 0;
        double mes_el = 0;
        double h_x, h_y, h_z;
        for (int i = 0; i < NUM_SPLIT_Z; i++) {
            for (int j = 0; j < NUM_SPLIT_Y; j++) {
                for (int k = 0; k < NUM_SPLIT_X; k++) {

                    if (IsFictEl(curr_el)) {
                        curr_el++;
                        continue;
                    }

                    if (max_t[curr_el][ph] < min_t) {                       
                        double sum = 0;
                        double nec_sum;
                        int S_g_1, S_g_2, S_g_3, S_g_4;
                        S_g_1 = (faces_flow_ph_value[start_face + k][ph] < 0) ? 1 : -1;
                        S_g_2 = (faces_flow_ph_value[start_face + k + NUM_SPLIT_X][ph] < 0) ? 1 : -1;
                        S_g_3 = (faces_flow_ph_value[start_face + k + NUM_SPLIT_X + 1][ph] < 0) ? -1 : 1;
                        S_g_4 = (faces_flow_ph_value[start_face + k + 2 * NUM_SPLIT_X + 1][ph] < 0) ? -1 : 1;

                        if (S_g_1 > 0)
                            sum += fabs(faces_flow_ph_value[start_face + k][ph]);
                        if (S_g_2 > 0)
                            sum += fabs(faces_flow_ph_value[start_face + k + NUM_SPLIT_X][ph]);
                        if (S_g_3 > 0)
                            sum += fabs(faces_flow_ph_value[start_face + k + NUM_SPLIT_X + 1][ph]);
                        if (S_g_4 > 0)
                            sum += fabs(faces_flow_ph_value[start_face + k + 2 * NUM_SPLIT_X + 1][ph]);

                        auto node_num = array_p[curr_el].second;

                        h_x = nodes[node_num[2]].x - nodes[node_num[0]].x;
                        h_y = nodes[node_num[6]].y - nodes[node_num[0]].y;
                        h_z = nodes[node_num[18]].z - nodes[node_num[0]].z;

                        int num_sub = array_p[curr_el].first;
                        mes_el = h_x * h_y * h_z;

                        nec_sum = mes_el * F * S_ph[curr_el][ph] / min_t;
                        double diff_sum = fabs(nec_sum - sum);

                        if (ph == 0) {
                            if (S_g_1 > 0) {
                                faces_flow_ph_value[start_face + k][0] -= diff_sum / 2.;
                                faces_flow_ph_value[start_face + k][1] += diff_sum / 2.;
                            }
                            if (S_g_2 > 0) {
                                faces_flow_ph_value[start_face + k + NUM_SPLIT_X][0] -= diff_sum / 2.;
                                faces_flow_ph_value[start_face + k + NUM_SPLIT_X][1] += diff_sum / 2.;
                            }
                            if (S_g_3 > 0) {
                                faces_flow_ph_value[start_face + k + NUM_SPLIT_X + 1][0] -= diff_sum / 2.;
                                faces_flow_ph_value[start_face + k + NUM_SPLIT_X + 1][1] += diff_sum / 2.;
                            }
                            if (S_g_4 > 0) {
                                faces_flow_ph_value[start_face + k + 2 * NUM_SPLIT_X + 1][0] -= diff_sum / 2.;
                                faces_flow_ph_value[start_face + k + 2 * NUM_SPLIT_X + 1][1] += diff_sum / 2.;
                            }
                            max_t[curr_el][ph] = min_t;
                        }
                        else if (ph == 1) {
                            if (S_g_1 > 0) {
                                faces_flow_ph_value[start_face + k][0] += diff_sum / 2.;
                                faces_flow_ph_value[start_face + k][1] -= diff_sum / 2.;
                            }
                            if (S_g_2 > 0) {
                                faces_flow_ph_value[start_face + k + NUM_SPLIT_X][0] += diff_sum / 2.;
                                faces_flow_ph_value[start_face + k + NUM_SPLIT_X][1] -= diff_sum / 2.;
                            }
                            if (S_g_3 > 0) {
                                faces_flow_ph_value[start_face + k + NUM_SPLIT_X + 1][0] += diff_sum / 2.;
                                faces_flow_ph_value[start_face + k + NUM_SPLIT_X + 1][1] -= diff_sum / 2.;
                            }
                            if (S_g_4 > 0) {
                                faces_flow_ph_value[start_face + k + 2 * NUM_SPLIT_X + 1][0] += diff_sum / 2.;
                                faces_flow_ph_value[start_face + k + 2 * NUM_SPLIT_X + 1][1] -= diff_sum / 2.;
                            }
                            max_t[curr_el][ph] = min_t;
                        }
                        else {
                            cout << "Не работает для фаз больше 2" << endl;
                        }                  
                    }
                    curr_el++;               
                }
                start_face += 2 * NUM_SPLIT_X + 1;
            }
            start_face += NUM_SPLIT_X;
        }
    }
}
#pragma endregion

#pragma region Максимальынй шаг по времени
double MaxDeltaT() { //максимальный шаг по времени (учитывая каждую фазу, каждый конечный элемент)
    double h_t = 0;
    double max_h_t = 1e+200;
    max_t.assign(NUM_END_EL, vector<double>(num_ph, 0.0));
    for (int ph = 0; ph < num_ph; ph++) {      
        int start_face = 0;
        int curr_el = 0;
        double mes_el = 0;
        double h_x, h_y, h_z;
        for (int i = 0; i < NUM_SPLIT_Z; i++) {
            for (int j = 0; j < NUM_SPLIT_Y; j++) {
                for (int k = 0; k < NUM_SPLIT_X; k++) {

                    if (IsFictEl(curr_el)) {
                        curr_el++;
                        continue;
                    }

                    double sum = 0;
                    int S_g_1, S_g_2, S_g_3, S_g_4;
                    S_g_1 = (faces_flow_ph_value[start_face + k][ph] < 0) ? 1 : -1;
                    S_g_2 = (faces_flow_ph_value[start_face + k + NUM_SPLIT_X][ph] < 0) ? 1 : -1;
                    S_g_3 = (faces_flow_ph_value[start_face + k + NUM_SPLIT_X + 1][ph] < 0) ? -1 : 1;
                    S_g_4 = (faces_flow_ph_value[start_face + k + 2 * NUM_SPLIT_X + 1][ph] < 0) ? -1 : 1;

                    if (S_g_1 > 0)
                        sum += fabs(faces_flow_ph_value[start_face + k][ph]);
                    if (S_g_2 > 0)
                        sum += fabs(faces_flow_ph_value[start_face + k + NUM_SPLIT_X][ph]);
                    if (S_g_3 > 0)
                        sum += fabs(faces_flow_ph_value[start_face + k + NUM_SPLIT_X + 1][ph]);
                    if (S_g_4 > 0)
                        sum += fabs(faces_flow_ph_value[start_face + k + 2 * NUM_SPLIT_X + 1][ph]);

                    auto node_num = array_p[curr_el].second;

                    h_x = nodes[node_num[2]].x - nodes[node_num[0]].x;
                    h_y = nodes[node_num[6]].y - nodes[node_num[0]].y;
                    h_z = nodes[node_num[18]].z - nodes[node_num[0]].z;

                    int num_sub = array_p[curr_el].first;
                    mes_el = h_x * h_y * h_z;

                    h_t = mes_el * F * S_ph[curr_el][ph] / sum;
                    max_t[curr_el][ph] = h_t;

                    if (h_t < max_h_t)
                        max_h_t = h_t;

                    curr_el++;
                }
                start_face += 2 * NUM_SPLIT_X + 1;
            }
            start_face += NUM_SPLIT_X;
        }
    }  
    return max_h_t;
}

double MaxDeltaTMixture() { //максимальный шаг по времени для смеси
    max_t_mixture.resize(NUM_END_EL);
    int start_face = 0;
    int curr_el = 0;
    double mes_el = 0;
    double h_x, h_y, h_z;
    double h_t = 0;
    double max_h_t = 1e+200;
    for (int i = 0; i < NUM_SPLIT_Z; i++) {
        for (int j = 0; j < NUM_SPLIT_Y; j++) {
            for (int k = 0; k < NUM_SPLIT_X; k++) {

                if (IsFictEl(curr_el)) {
                    curr_el++;
                    continue;
                }

                double sum = 0;
                int S_g_1, S_g_2, S_g_3, S_g_4;
                S_g_1 = (faces_flow_value[start_face + k] < 0) ? 1 : -1;
                S_g_2 = (faces_flow_value[start_face + k + NUM_SPLIT_X] < 0) ? 1 : -1;
                S_g_3 = (faces_flow_value[start_face + k + NUM_SPLIT_X + 1] < 0) ? -1 : 1;
                S_g_4 = (faces_flow_value[start_face + k + 2 * NUM_SPLIT_X + 1] < 0) ? -1 : 1;

                if (S_g_1 > 0)
                    sum += fabs(faces_flow_value[start_face + k]);
                if (S_g_2 > 0)
                    sum += fabs(faces_flow_value[start_face + k + NUM_SPLIT_X]);
                if (S_g_3 > 0)
                    sum += fabs(faces_flow_value[start_face + k + NUM_SPLIT_X + 1]);
                if (S_g_4 > 0)
                    sum += fabs(faces_flow_value[start_face + k + 2 * NUM_SPLIT_X + 1]);

                auto node_num = array_p[curr_el].second;

                h_x = nodes[node_num[2]].x - nodes[node_num[0]].x;
                h_y = nodes[node_num[6]].y - nodes[node_num[0]].y;
                h_z = nodes[node_num[18]].z - nodes[node_num[0]].z;

                int num_sub = array_p[curr_el].first;
                mes_el = h_x * h_y * h_z;

                h_t = mes_el * F / sum;
                max_t_mixture[curr_el] = h_t;

                if (h_t < max_h_t)
                    max_h_t = h_t;
                curr_el++;
            }
            start_face += 2 * NUM_SPLIT_X + 1;
        }
        start_face += NUM_SPLIT_X;
    }
    return max_h_t;
}
#pragma endregion

#pragma region Расчет новых потоков и насыщенностей
void CalcNewFlowAndSatur(double h_t) { // расчитать новые потоки на элементах и значения насыщенностей
    //satur.assign(NUM_SPLIT_Z * NUM_SPLIT_Y * NUM_SPLIT_X, vector<double>(num_ph, 0.0));
    for (int ph = 0; ph < num_ph; ph++) {     
        int start_face = 0;
        int curr_el = 0;
        double mes_el = 0;
        double h_x, h_y, h_z;
        for (int i = 0; i < NUM_SPLIT_Z; i++) {
            for (int j = 0; j < NUM_SPLIT_Y; j++) {
                for (int k = 0; k < NUM_SPLIT_X; k++) {

                    if (IsFictEl(curr_el)) {
                        curr_el++;
                        continue;
                    }

                    double sum = 0;
                    int S_g_1, S_g_2, S_g_3, S_g_4;
                    S_g_1 = (faces_flow_ph_value[start_face + k][ph] < 0) ? 1 : -1;
                    S_g_2 = (faces_flow_ph_value[start_face + k + NUM_SPLIT_X][ph] < 0) ? 1 : -1;
                    S_g_3 = (faces_flow_ph_value[start_face + k + NUM_SPLIT_X + 1][ph] < 0) ? -1 : 1;
                    S_g_4 = (faces_flow_ph_value[start_face + k + 2 * NUM_SPLIT_X + 1][ph] < 0) ? -1 : 1;

                    if (S_g_1 > 0)
                        sum -= fabs(faces_flow_ph_value[start_face + k][ph]) * h_t;
                    else
                        sum += fabs(faces_flow_ph_value[start_face + k][ph]) * h_t;
                    if (S_g_2 > 0)
                        sum -= fabs(faces_flow_ph_value[start_face + k + NUM_SPLIT_X][ph]) * h_t;
                    else
                        sum += fabs(faces_flow_ph_value[start_face + k + NUM_SPLIT_X][ph]) * h_t;
                    if (S_g_3 > 0)
                        sum -= fabs(faces_flow_ph_value[start_face + k + NUM_SPLIT_X + 1][ph]) * h_t;
                    else
                        sum += fabs(faces_flow_ph_value[start_face + k + NUM_SPLIT_X + 1][ph]) * h_t;
                    if (S_g_4 > 0)
                        sum -= fabs(faces_flow_ph_value[start_face + k + 2 * NUM_SPLIT_X + 1][ph]) * h_t;
                    else
                        sum += fabs(faces_flow_ph_value[start_face + k + 2 * NUM_SPLIT_X + 1][ph]) * h_t;

                    auto node_num = array_p[curr_el].second;

                    h_x = nodes[node_num[2]].x - nodes[node_num[0]].x;
                    h_y = nodes[node_num[6]].y - nodes[node_num[0]].y;
                    h_z = nodes[node_num[18]].z - nodes[node_num[0]].z;

                    int num_sub = array_p[curr_el].first;
                    mes_el = h_x * h_y * h_z;

                    sum += mes_el * F * S_ph[curr_el][ph];

                    S_ph[curr_el][ph] = sum / (mes_el * F);

                    curr_el++;
                }
                start_face += 2 * NUM_SPLIT_X + 1;
            }
            start_face += NUM_SPLIT_X;
        }
    }
}
#pragma endregion

#pragma region Пересчет состояния ячеек
void RecalcCellState() { // пересчитать состояние ячеек
    CalcFlowPh();
    double max_h_t = MaxDeltaTMixture();
    double h_t = MaxDeltaT();
    int n = 10;
    double min_t = max_h_t / n;
    RecalcFlowPh(max_h_t);
    for (int i = 1; i <= n; i++) {
        CalcNewFlowAndSatur(i * min_t);
        CalcFlowPh();
        cout << S_ph[36][0] << " " << S_ph[36][1] << endl;
        cout << S_ph[37][0] << " " << S_ph[37][1] << endl;
        cout << S_ph[38][0] << " " << S_ph[38][1] << endl;
        cout << S_ph[39][0] << " " << S_ph[39][1] << endl << endl << endl;
    }
}
#pragma endregion