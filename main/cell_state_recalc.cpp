#include "cell_state_recalc.h"
#include "flow_balancing.h"
#include "solver_slae.h"

#pragma region ������ ������� ��� ������ ����
void CalcFlowPh() { // ��������� ������ ��� ������� ���� � �����������
    // ������� ����������: 
    // ������� ����� �� y, ����� �� x (����� - �������; �� ��������� �� ������� (���� �������� � ������� xy �� ����� �����)), ����� ����� �� z
    // (������ �� ������� � ������ z �������� � �� ����� ����������, � ���������� ������ �� z ����������� ������� - ����� �� ��������� (�������� ��������� �������� �� ������� z))
    faces_flow_ph_value.assign(faces_flow_value.size(), vector<double>(num_ph));
    vector<int> node_num(27);
    int num_sub;
    int num_pos_flow = 0;
    int start_el = 0;
    int num_flow = 0;
    double sum;
    for (int i = 0; i < NUM_SPLIT_Z; i++) {

        //����� �������� ����� �� y
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
        //����� �������� ����� �� y

        //����� ����� ����� �� x
        num_sub = array_p[start_el].first;
        sum = 0.;
        for (int ph = 0; ph < num_ph; ph++) {
            sum += k_ph(start_el, ph) / eta_ph[ph];
        }
        for (int ph = 0; ph < num_ph; ph++) {
            faces_flow_ph_value[num_flow][ph] = faces_flow_value[num_flow] * k_ph(start_el, ph) / (eta_ph[ph] * sum);
        }
        num_flow++;
        //����� ����� ����� �� x

        //���������� ����� �� x
        for (int j = 0; j < NUM_SPLIT_X - 1; j++) {
            int S_g = (faces_flow_value[num_flow] < -EPS) ? -1 : 1;
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
        //���������� ����� �� x

        //����� ������ ����� �� x
        num_sub = array_p[start_el + NUM_SPLIT_X - 1].first;
        sum = 0.;
        for (int ph = 0; ph < num_ph; ph++) {
            sum += k_ph(start_el + NUM_SPLIT_X - 1, ph) / eta_ph[ph];
        }
        for (int ph = 0; ph < num_ph; ph++) {
            faces_flow_ph_value[num_flow][ph] = faces_flow_value[num_flow] * k_ph(start_el + NUM_SPLIT_X - 1, ph) / (eta_ph[ph] * sum);
        }
        num_flow++;
        //����� ������ ����� �� x

        start_el += NUM_SPLIT_X;

        //�������� ����
        for (int j = 1; j < NUM_SPLIT_Y; j++) {

            //���������� ����� �� y
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
                    else {
                        S_ph[start_el + k][0] = 0;
                        S_ph[start_el + k][1] = 0;
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
                    else {
                        S_ph[start_el + k - NUM_SPLIT_X][0] = 0;
                        S_ph[start_el + k - NUM_SPLIT_X][1] = 0;
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
                else if (IsFictEl(start_el + k) && IsFictEl(start_el + k - NUM_SPLIT_X)) {
                    num_flow++;
                }
                else {
                    int S_g = (faces_flow_value[num_flow] < -EPS) ? -1 : 1;
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
            //���������� ����� �� y

            //����� ����� ����� �� x
            num_sub = array_p[start_el].first;
            sum = 0.;
            for (int ph = 0; ph < num_ph; ph++) {
                sum += k_ph(start_el, ph) / eta_ph[ph];
            }
            for (int ph = 0; ph < num_ph; ph++) {
                faces_flow_ph_value[num_flow][ph] = faces_flow_value[num_flow] * k_ph(start_el, ph) / (eta_ph[ph] * sum);
            }
            num_flow++;
            //����� ����� ����� �� x

            //���������� ����� �� x
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
                    else {
                        S_ph[start_el + k + 1][0] = 0;
                        S_ph[start_el + k + 1][1] = 0;
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
                    else {
                        S_ph[start_el + k][0] = 0;
                        S_ph[start_el + k][1] = 0;
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
                    int S_g = (faces_flow_value[num_flow] < -EPS) ? -1 : 1;
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
            //���������� ����� �� x

            //����� ������ ����� �� x
            num_sub = array_p[start_el + NUM_SPLIT_X - 1].first;
            sum = 0.;
            for (int ph = 0; ph < num_ph; ph++) {
                sum += k_ph(start_el + NUM_SPLIT_X - 1, ph) / eta_ph[ph];
            }
            for (int ph = 0; ph < num_ph; ph++) {
                faces_flow_ph_value[num_flow][ph] = faces_flow_value[num_flow] * k_ph(start_el + NUM_SPLIT_X - 1, ph) / (eta_ph[ph] * sum);
            }
            num_flow++;
            //����� ������ ����� �� x

            if (j == NUM_SPLIT_Y - 1) continue;
            start_el += NUM_SPLIT_X;
        }
        //�������� ����

        //����� ������� ����� �� y
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
        //����� ������� ����� �� y

        start_el += NUM_SPLIT_X;
    }
}
#pragma endregion

#pragma region ������������ ��� �� �������
double MaxDeltaT() { //������������ ��� �� ������� (�������� ������ ����, ������ �������� �������)
    double h_t = 0;
    double max_h_t = 1e+200;
    int num = 0;
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

                    if (h_t < max_h_t) {
                        max_h_t = h_t;
                        num = curr_el;
                    }
                       
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

#pragma region ������ ����� ������� � �������������
void CalcNewFlowAndSatur(double h_t) { // ��������� ����� ������ �� ��������� � �������� �������������
    int start_face = 0;
    int curr_el = 0;
    double mes_el = 0;
    double h_x, h_y, h_z;
    for (int i = 0; i < NUM_SPLIT_Z; i++) {
        for (int j = 0; j < NUM_SPLIT_Y; j++) {
            for (int k = 0; k < NUM_SPLIT_X; k++) {

                vector<double> V_ph(num_ph);
                double sum_V_ph = 0.;

                for (int ph = 0; ph < 2; ph++) {
                    if (IsFictEl(curr_el)) {
                        continue;
                    }

                    double sum = 0;
                    int S_g_1, S_g_2, S_g_3, S_g_4;
                    S_g_1 = (faces_flow_ph_value[start_face + k][ph] < 0) ? 1 : -1;
                    S_g_2 = (faces_flow_ph_value[start_face + k + NUM_SPLIT_X][ph] < 0) ? 1 : -1;
                    S_g_3 = (faces_flow_ph_value[start_face + k + NUM_SPLIT_X + 1][ph] < 0) ? -1 : 1;
                    S_g_4 = (faces_flow_ph_value[start_face + k + 2 * NUM_SPLIT_X + 1][ph] < 0) ? -1 : 1;

                    if (S_g_1 > 0)
                        V_ph[ph] -= fabs(faces_flow_ph_value[start_face + k][ph]) * h_t;
                    else
                        V_ph[ph] += fabs(faces_flow_ph_value[start_face + k][ph]) * h_t;
                    if (S_g_2 > 0)
                        V_ph[ph] -= fabs(faces_flow_ph_value[start_face + k + NUM_SPLIT_X][ph]) * h_t;
                    else
                        V_ph[ph] += fabs(faces_flow_ph_value[start_face + k + NUM_SPLIT_X][ph]) * h_t;
                    if (S_g_3 > 0)
                        V_ph[ph] -= fabs(faces_flow_ph_value[start_face + k + NUM_SPLIT_X + 1][ph]) * h_t;
                    else
                        V_ph[ph] += fabs(faces_flow_ph_value[start_face + k + NUM_SPLIT_X + 1][ph]) * h_t;
                    if (S_g_4 > 0)
                        V_ph[ph] -= fabs(faces_flow_ph_value[start_face + k + 2 * NUM_SPLIT_X + 1][ph]) * h_t;
                    else
                        V_ph[ph] += fabs(faces_flow_ph_value[start_face + k + 2 * NUM_SPLIT_X + 1][ph]) * h_t;

                    

                    auto node_num = array_p[curr_el].second;

                    h_x = nodes[node_num[2]].x - nodes[node_num[0]].x;
                    h_y = nodes[node_num[6]].y - nodes[node_num[0]].y;
                    h_z = nodes[node_num[18]].z - nodes[node_num[0]].z;

                    int num_sub = array_p[curr_el].first;
                    mes_el = h_x * h_y * h_z;

                    V_ph[ph] += mes_el * F * S_ph[curr_el][ph];
                    if (V_ph[ph] < 0) {
                        cout << "������!!!" << endl;
                    }
                    sum_V_ph += V_ph[ph];                
                }

                for (int ph = 0; ph < 2; ph++) {

                    if (IsFictEl(curr_el)) {
                        continue;
                    }

                    S_ph[curr_el][ph] = V_ph[ph] / sum_V_ph;
                }
                curr_el++;
            }
            start_face += 2 * NUM_SPLIT_X + 1;
        }
        start_face += NUM_SPLIT_X;
    }
    
}

vector<double> CalcExtraction(double h_t) { // ���������� ����� ������ ����
    vector<double> extraction_ph_value(num_ph);
    int num_sub;
    int num_pos_flow = 0;
    int start_el = 0;
    int num_flow = 0;
    double sum;
    for (int i = 0; i < NUM_SPLIT_Z; i++) {

        //����� �������� ����� �� y
        for (int j = 0; j < NUM_SPLIT_X; j++) {
            num_flow++;
        }
        //����� �������� ����� �� y

        //����� ����� ����� �� x
        num_flow++;
        //����� ����� ����� �� x

        //���������� ����� �� x
        for (int j = 0; j < NUM_SPLIT_X - 1; j++) {
            num_flow++;
        }
        //���������� ����� �� x

        //����� ������ ����� �� x
        num_flow++;
        //����� ������ ����� �� x

        start_el += NUM_SPLIT_X;

        //�������� ����
        for (int j = 1; j < NUM_SPLIT_Y; j++) {

            //���������� ����� �� y
            for (int k = 0; k < NUM_SPLIT_X; k++) {

                if (IsFictEl(start_el + k) && !IsFictEl(start_el + k - NUM_SPLIT_X)) {
                    if (faces_flow_value[num_flow] > EPS) {                       
                        for (int ph = 0; ph < num_ph; ph++) {
                            extraction_ph_value[ph] += fabs(faces_flow_ph_value[num_flow][ph]) * h_t;
                        }
                        num_flow++;
                        continue;
                    }
                    num_flow++;
                }
                else if (!IsFictEl(start_el + k) && IsFictEl(start_el + k - NUM_SPLIT_X)) {
                    if (faces_flow_value[num_flow] < -EPS) {
                        for (int ph = 0; ph < num_ph; ph++) {
                            extraction_ph_value[ph] += fabs(faces_flow_ph_value[num_flow][ph]) * h_t;
                        }
                        num_flow++;
                        continue;
                    }
                    num_flow++;
                }
                else if (IsFictEl(start_el + k) && IsFictEl(start_el + k - NUM_SPLIT_X)) {
                    num_flow++;
                }
                else {
                    num_flow++;
                }
            }
            //���������� ����� �� y

            //����� ����� ����� �� x
            num_flow++;
            //����� ����� ����� �� x

            //���������� ����� �� x
            for (int k = 0; k < NUM_SPLIT_X - 1; k++) {
                if (!IsFictEl(start_el + k) && IsFictEl(start_el + k + 1)) {
                    if (faces_flow_value[num_flow] > EPS) {
                        for (int ph = 0; ph < num_ph; ph++) {
                            extraction_ph_value[ph] += fabs(faces_flow_ph_value[num_flow][ph]) * h_t;
                        }
                        num_flow++;
                        continue;
                    }
                    num_flow++;
                }
                else if (IsFictEl(start_el + k) && !IsFictEl(start_el + k + 1)) {
                    if (faces_flow_value[num_flow] < -EPS) {
                        for (int ph = 0; ph < num_ph; ph++) {
                            extraction_ph_value[ph] += fabs(faces_flow_ph_value[num_flow][ph]) * h_t;
                        }
                        num_flow++;
                        continue;
                    }
                    num_flow++;
                }
                else if (IsFictEl(start_el + k) && IsFictEl(start_el + k + 1)) {
                    num_flow++;
                }
                else {
                    num_flow++;
                }
            }
            //���������� ����� �� x

            //����� ������ ����� �� x
            num_flow++;
            //����� ������ ����� �� x

            if (j == NUM_SPLIT_Y - 1) continue;
            start_el += NUM_SPLIT_X;
        }
        //�������� ����

        //����� ������� ����� �� y
        for (int j = 0; j < NUM_SPLIT_X; j++) {
            num_flow++;
        }
        //����� ������� ����� �� y

        start_el += NUM_SPLIT_X;
    }
    return extraction_ph_value;
}
#pragma endregion

#pragma region �������� ��������� �����
void RecalcCellState() { // ����������� ��������� �����
    CalcFlowPh();
    double total_time = 0.;
    double time_to_recalc_P = 0.;
    double curr_h_t = MaxDeltaT();
    double extr_oil = 0., extr_water = 0.;

    ofstream f1("extr_oil_neod.txt");
    ofstream f2("perc_oil_neod.txt");
    ofstream f3("time_neod.txt");
    while (total_time < 50 * 86400.) {
    
        total_time += curr_h_t;
        time_to_recalc_P += curr_h_t;
        f3 << total_time / 86400. << endl;

        CalcNewFlowAndSatur(curr_h_t);
        CalcFlowPh();

        auto extr = CalcExtraction(curr_h_t);
        extr_oil += extr[0];
        f1 << extr_oil << endl;

        extr_water += extr[1];
        f2 << extr_oil / (extr_oil + extr_water) * 100. << endl;

        curr_h_t = MaxDeltaT();
    }
    f1.close();
    f2.close();
    f3.close();
}
#pragma endregion