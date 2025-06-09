#pragma once
#include "fem.h"

#pragma region ������������ ��� ������������ �������� �������
void RelativeErrorNorm(); // ������������� ����� �����������
vector<double> grad_u(double x, double y, double z); // �������� ������������� �������
vector<double> grad_uh(double xi, double eta, double zeta, int num_end_el); // �������� �������, ��������� ��������
void EnergyNorm(); // �������������� �����
#pragma endregion

#pragma region ��������� � ������������� ��������
double AnalitP(double x, double y, double z);
void VecAnalitP();
#pragma endregion