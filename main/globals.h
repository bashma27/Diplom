#pragma once
#include "types.h"

extern double EPS, F;
extern vector<vector<Coord2>> zone_perf; // массив зон перфорации
extern vector<pair<int, vector<double>>> W;
extern int NUM_ZONE_PERF; // количество зон перфорации
extern int NUM_SPLIT_X, NUM_SPLIT_Y, NUM_SPLIT_Z; // суммарное количество разбиений по x, y, z
extern int NUM_NODES_IN_EDGE_X, NUM_NODES_IN_EDGE_Y, NUM_NODES_IN_EDGE_Z, NUM_NODES, NUM_END_EL;
extern int num_ph; // кол-во фаз
extern vector<double> set_flow_zp; // заданный поток зон перфорации
extern vector<int> ia, ja, choice;
extern vector<double> aal, di, b, q, L_sq, di_sq, normal;
extern vector<vector<double>> S_ph; // значения насыщенности фаз для каждого элемента
extern vector<vector<double>> S_ph_subarea; // значения насыщенности фаз для подоблостей
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
extern vector<vector<double>> faces_flow_ph_value; // массив граней со значением потока фазы
extern vector<pair<int, int>> num_faces_zp; // грани зон перфорации
extern vector<int> fict_faces; // номера фиктивных граней
extern vector<int> fict_el; // номера фиктивных конечных элементов
extern vector<double> q_flow;
extern vector<vector<int>> list_face; 
extern vector<int> el_for_face; // конечный элемент, соответствующий грани
extern vector<vector<double>> satur; // значения насыщенностей для каждого элемента каждой фазы
extern vector<vector<double>> max_t; // максимальные значения времени для перетоков каждого элемента каждой фазы 
extern vector<double> max_t_mixture; // максимальные значения времени для смеси
extern vector<vector<double>> satur_pos_set_flow; // насыщенности каждой фазы для заданных зон перфорации, закачивающих жидкость 
extern vector<vector<pair<int, double>>> coef; // Коэффициеты n, q  для каждой области в виде пары
extern double theta_;