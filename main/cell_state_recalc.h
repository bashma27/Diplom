#pragma once
#include "globals.h"
#include "get_funcs.h"
#include "fem.h"

#pragma region ������ ������� ��� ������� ����
void CalcFlowPh(); // ��������� ������ ��� ������� ���� � �����������
#pragma endregion

#pragma region ������������ ��� �� �������
double MaxDeltaT(); //������������ ��� �� ������� (�������� ������ ����, ������ �������� �������)
#pragma endregion

#pragma region ������ ����� ������� � �������������
void CalcNewFlowAndSatur(double h_t); // ��������� ����� ������ �� ��������� � �������� �������������
#pragma endregion

#pragma region �������� ��������� �����
void RecalcCellState(); // ����������� ��������� �����
#pragma endregion