#pragma once
#include "globals.h"
#include "get_funcs.h"
#include "fem.h"

#pragma region ������ ������� ��� ������ ����
void CalcFlowPh(); // ��������� ������ ��� ������� ���� � �����������
void RecalcFlowPh(double min_t); // ����������� ������� ��� ������������ ���� �� �������
#pragma endregion

#pragma region ������������ ��� �� �������
double MaxDeltaT(); //������������ ��� �� ������� (�������� ������ ����, ������ �������� �������)
double MaxDeltaTMixture(); //������������ ��� �� ������� ��� �����
#pragma endregion

#pragma region ������ ����� ������� � �������������
void CalcNewFlowAndSatur(double h_t); // ��������� ����� ������ �� ��������� � �������� �������������
#pragma endregion

#pragma region ������� ���������� ������ ��� � �������� ������
void CalcSumVPhAndPoreV(double max_h_t); // ��������� ��������� ����� ��� � ������� �����
#pragma endregion

#pragma region �������� ��������� �����
void RecalcCellState(); // ����������� ��������� �����
#pragma endregion