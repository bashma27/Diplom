#pragma once
#include "globals.h"
#include "fem.h"

#pragma region ������ �������
double CalcFlowF(vector<int>& node_num, int num_sub); // ������ ������, ������������ ����������
double CalcFlowFaceYZ(vector<int>& node_num, int num_sub, int n); // ������ ������ ����� ����� yz
double CalcFlowFaceXZ(vector<int>& node_num, int num_sub, int n); // ������ ������ ����� ����� xz
double CalcFlowFaceXY(vector<int>& node_num, int num_sub, int n); // ������ ������ ����� ����� xy
#pragma endregion

#pragma region ��������� ������ �� ��������� ������
void GenFacesFlowValue(); // ��������� ������ �� ��������� ������
#pragma endregion

#pragma region ��������� ������ � ������� ����� ������� ��� �������
void ClearAndGenPortMatr(); // ��������� ������ � ������� ����� ������� ��� �������
#pragma endregion

#pragma region ��������� ������� � �������
void GenMatrAndVec(); // ��������� ������� � ������� 
#pragma endregion

#pragma region ���� �������
void ConsiderKnownFlows(int n); // ���� ��������� �������
void ConsiderFictFlows(int n); // ���� ��������� �������
void ConsiderFlows(); // ���� ���� �������
#pragma endregion

#pragma region ������ ��������� � ���������� ���������������� �������
double CalcSumNonBalance(); // ���������� ��������� ��������
void FindBalancedFlows(); // ����� ���������������� ������
#pragma endregion

#pragma region ������������ �������
void BalancingFlows(); // ������������ �������
#pragma endregion
