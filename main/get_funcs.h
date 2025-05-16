#pragma once
#include "globals.h"

int GetNumSubarea(vector<int>& end_el); // �������� ����� ���������� ��� ��������� ��������
int GetNumNodes(vector<Coord2>& zone_perf_i, int i_x, int i_y); // �������� ����� ���� ��� ���� ����������
int GetNumSubareaBound(int ind, vector<int>& face_2_zp); // �������� ����� ���������� ��� ������� �������� (������ ������ 0/1 - ����� ������������ ��� y/x)
int GetNumEndEl(Coord3 point); // �������� ����� ��������� �������� ��� �����
int GetNumFace2ZP(vector<int>& _face_2_zp); // �������� ����� ����� ���� ����������
