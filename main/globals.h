#pragma once
#include "types.h"

extern double EPS, F;
extern vector<vector<Coord2>> zone_perf; // ������ ��� ����������
extern vector<pair<int, vector<double>>> W;
extern int NUM_ZONE_PERF; // ���������� ��� ����������
extern int NUM_SPLIT_X, NUM_SPLIT_Y, NUM_SPLIT_Z; // ��������� ���������� ��������� �� x, y, z
extern int NUM_NODES_IN_EDGE_X, NUM_NODES_IN_EDGE_Y, NUM_NODES_IN_EDGE_Z, NUM_NODES, NUM_END_EL;
extern int num_ph; // ���-�� ���
extern vector<double> set_flow_zp; // �������� ����� ��� ����������
extern vector<int> ia, ja, choice;
extern vector<double> aal, di, b, q, L_sq, di_sq, normal;
extern vector<vector<double>> S_ph; // �������� ������������ ��� ��� ������� ��������
extern vector<vector<double>> S_ph_subarea; // �������� ������������ ��� ��� �����������
extern vector<double> K, eta_ph; // ������� ������������� ����������� �������������, ������������� ������������ �������� ��������������
extern vector<Coord3> nodes; // ���� �����, �������� ������������ 
extern unordered_set<int> face_1; // ������ � ������ ������ �������
extern vector<pair<int, vector<int>>> face_2_zp; // ������ ������ � �������� 2 ���� ��� ����������
extern vector<pair<int, vector<int>>> array_p; // ������ ����������������
extern vector<vector<double>> G, M; // ��������� ������� ���� � ���������
extern vector<int> ident_fict; // �������������� ����������� (������� ��� ��������� ��������)
extern vector<int> i_ident_fict;
extern vector<int> fict_nodes; // ������ ��������� �����
extern vector<double> faces_flow_value; // ������ ������ �� ��������� ������
extern vector<vector<double>> faces_flow_ph_value; // ������ ������ �� ��������� ������ ����
extern vector<pair<int, int>> num_faces_zp; // ����� ��� ����������
extern vector<int> fict_faces; // ������ ��������� ������
extern vector<int> fict_el; // ������ ��������� �������� ���������
extern vector<double> q_flow;
extern vector<vector<int>> list_face; 
extern vector<int> el_for_face; // �������� �������, ��������������� �����
extern vector<vector<double>> satur; // �������� ������������� ��� ������� �������� ������ ����
extern vector<vector<double>> max_t; // ������������ �������� ������� ��� ��������� ������� �������� ������ ���� 
extern vector<double> max_t_mixture; // ������������ �������� ������� ��� �����
extern vector<vector<double>> satur_pos_set_flow; // ������������ ������ ���� ��� �������� ��� ����������, ������������ �������� 
extern vector<vector<pair<int, double>>> coef; // ����������� n, q  ��� ������ ������� � ���� ����
extern double theta_;