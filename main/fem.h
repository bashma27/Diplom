#pragma once
#include "globals.h"

#pragma region ��������������� �������
double BasicFunc(int i, double input_point, double point_1, double point_2, double point_3); // ������������ �������� �������
double DerivativeLocalBasicFuncs(int i, double psi); // ����������� �� �������� ��������
bool IsFictitious(int num); // ����������� ����������� ����
void GenFictNodes(); // ��������� ��������� �����
void GenArrayFictEndEl(); // ��������� ������� ��������� �������� ���������
bool IsFictEl(int num_end_el); // ����������� ����������� ��������� ��������
bool IsFindFaceZP(int num_face); // ���������� ����� � ������� ������ ��� ����������
bool IsFindFictFace(int num_face); // ���������� ����� � ������� ��������� ������
#pragma endregion

#pragma region ������� ������� �������, �� ������� � ������� ������ ����� f
double u_g(double x, double y, double z); // ������� ������� ������� ����
double k_ph(int num_el, int num_ph); // ��������� ����������� �������������
double lambda(int num_sub); // ������ ������ ��� ��������� ������������
double theta(int num_face_2_zp); // ���������� ����� �� ��������� ������
double SpecFlow(int num_face_2_zp); // �������� ����� �� ������ ����� ��� ���� ����������
double f(int num_sub, double x, double y); // ������� ������ ����� ���. ���������
#pragma endregion

#pragma region ��������� ������� �������� ��������� � �������� �������
void GenArrayParallelepipeds(); // ��������� ������� ����������������
void GenPortraitMatr(); // ��������� �������� �������
#pragma endregion

#pragma region ������� ���������� � ����. ������� � ������ ������ �����
void AddLocalMatr(vector<int>& node_num, vector<vector<double>>& loc_matr); // ���������� ������� � ����������
void AddLocalVec(vector<int>& node_num, vector<double>& loc_vec); // ���������� ������� � ����������
void AddLocalVecBound(int num_face_2_zp, vector<double>& vec); // ���������� ������� �� ������� �������� � ����������
#pragma endregion

#pragma region ������ � �������� ��������� (��������� � ����)
void GenFirstBoundCondit(); // �������� ������� � ������� �������� ���������
void GenSecBoundCondit(); // �������� ������� � ������� �������� ���������
void ConsiderBoundConditFirstType(int n, vector<vector<pair<int, int>>>& columnToRows); // ���� ������� ������� ������� ����
void ConsiderBoundConditSecType(int num_face_2_zp); // ���� ������� ������� ������� ����
void ConsiderFictitiousNodes(); // ���� ��������� �����
void ConsiderBoundCondit(); // ���� ���� �������
#pragma endregion

#pragma region ������ � ���������� ��������� � ���������
void CreateLocalMatrStiff(int num_sub, vector<int>& node_num); // ��������� ������� ���������
void CreateLocalVec(int num_sub, vector<int>& node_num); // ��������� ������ ������ �����
#pragma endregion

#pragma region ���������� ���������� ������� � �������
void BuildMatrVec(); // ��������� ���������� ������� � ������
#pragma endregion

#pragma region ������������ ��� ������������ �������� �������
void Test(); // ����
#pragma endregion

#pragma region ������� �������� ������� � ����� �����
double GetResUInPoint(Coord3 point); // �������� �������� ���������� ������� � �����
#pragma endregion

#pragma region ����� ������� ���������� ������� ���� � ����
void OutputResult(); // ����� ���������� � ����
#pragma endregion