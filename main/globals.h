#pragma once
#include "types.h"

extern double EPS;
extern vector<vector<Coord2>> zone_perf; // массив зон перфорации
extern vector<pair<int, vector<double>>> W;
extern int NUM_ZONE_PERF; // количество зон перфорации
extern int NUM_SPLIT_X, NUM_SPLIT_Y, NUM_SPLIT_Z; // суммарное количество разбиений по x, y, z
extern int NUM_NODES_IN_EDGE_X, NUM_NODES_IN_EDGE_Y, NUM_NODES_IN_EDGE_Z, NUM_NODES;
extern int num_ph;
extern vector<double> set_flow_zp; // заданный поток зон перфорации
extern vector<int> ia, ja, choice;
extern vector<double> aal, di, b, q, L_sq, di_sq, normal;
extern vector<vector<double>> k_ph; //массив коэффициентов множителей структурной проницаемости,
extern vector<double> K, eta_ph; // массивы коэффициентов структурной проницаемости, коэффициентов динамической вязкости соответственно
extern vector<Coord3> nodes; // узлы сетки, заданные координатами 
extern unordered_set<int> face_1; // массив с узлами первых краевых
extern vector<pair<int, vector<int>>> face_2_zp; // массив граней с краевыми 2 рода зон перфорации
extern vector<pair<int, vector<int>>> array_p; // массив параллелепипедов
extern vector<vector<double>> G, M; // локальные матрицы масс и жесткости
extern vector<int> ident_fict; // идентификаторы фиктивности (области или конечного элемента)
extern vector<int> i_ident_fict;
extern vector<int> fict_nodes; // массив фиктивных узлов
extern vector<double> faces_flow_value; // массив граней со значением потока
extern vector<pair<int, int>> num_faces_zp; // грани зон перфорации
extern vector<int> fict_faces; // номера фиктивных граней
extern vector<int> fict_el; // номера фиктивных конечных элементов
extern vector<double> q_flow;
extern vector<vector<int>> list_face; 