#include "cell_state_recalc.h"
#include "flow_balancing.h"

#pragma region Расчет потоков для каждоый фазы
void CalcFlowPh() { // расчитать потоки для каждоый фазы в отдельности
    // порядок заполнения: 
    // сначала грани по y, потом по x (слева - направо; от ближайшей до дальней (если смотреть в сечении xy то снизу вверх)), потом вверх по z
    // (потоки по верхней и нижней z известны и не нужна фактически, а внутренние потоки по z практически нулевые - можно не учитывать (скважина полностью проходит по глубине z))
    faces_flow_ph_value.assign(faces_flow_value.size(), vector<double>(num_ph));
    int num_sub;
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

                    h_t = mes_el * F * k_ph(curr_el, ph) / sum;
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
#pragma endregion

#pragma region Расчет новых потоков и насыщенностей
void CalcNewFlowAndSatur(double h_t) { // расчитать новые потоки на элементах и значения насыщенностей
    satur.assign(NUM_SPLIT_Z * NUM_SPLIT_Y * NUM_SPLIT_X, vector<double>(num_ph, 0.0));
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

                    sum += mes_el * F * k_ph(curr_el, ph);

                    satur[curr_el][ph] = sum / (mes_el * F);

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
    double h_t = MaxDeltaT();
    CalcNewFlowAndSatur(h_t);
}
#pragma endregion