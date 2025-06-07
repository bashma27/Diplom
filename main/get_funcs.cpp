#include "get_funcs.h"

int GetNumSubarea(vector<int>& end_el) { // получить номер подобласти для конечного элемента
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

int GetNumNodes(vector<Coord2>& zone_perf_i, int i_x, int i_y) { // получить номер узла для зоны перфорации
    for (int i = 0; i < NUM_NODES_IN_EDGE_X * NUM_NODES_IN_EDGE_Y - 1; i++) {
        if ((zone_perf_i[i_x].x == nodes[i].x || abs(zone_perf_i[i_x].x - nodes[i].x) < EPS) &&
            (zone_perf_i[i_y].y == nodes[i].y || abs(zone_perf_i[i_y].y - nodes[i].y) < EPS)) return i;
    }
}

int GetNumSubareaBound(int ind, vector<int>& face_2_zp) { // получить номер подобласти для второго краевого (первый индекс 0/1 - грань параллельная оси y/x)
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

int GetNumEndEl(Coord3 point) { // получить номер конечного элемента для точки
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

int GetNumFace2ZP(vector<int>& _face_2_zp) { // получить номер грани зоны перфорации
    for (int i = 0; i < face_2_zp.size(); i++) {
        vector<int> vec;
        vec.assign(face_2_zp[i].second.begin(), face_2_zp[i].second.begin() + 9);
        if (vec == _face_2_zp) {
            return i;
        }
    }
}

int GetNumPosFlowFromAllSetFlows(int num_zone_perf) { // получить номер положительного потока из всех заданных потоков на зонах перфорации
    int num = 0;
    for (int i = 0; i < num_zone_perf; i++) {
        if (set_flow_zp[i] > EPS) {
            num++;
        }   
    }
    return num;
}