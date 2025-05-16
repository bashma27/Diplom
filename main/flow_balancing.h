#pragma once
#include "globals.h"
#include "fem.h"

#pragma region Расчет потоков
double CalcFlowF(vector<int>& node_num, int num_sub); // расчет потока, создаваемого источником
double CalcFlowFaceYZ(vector<int>& node_num, int num_sub, int n); // расчет потока через грань yz
double CalcFlowFaceXZ(vector<int>& node_num, int num_sub, int n); // расчет потока через грань xz
double CalcFlowFaceXY(vector<int>& node_num, int num_sub, int n); // расчет потока через грань xy
#pragma endregion

#pragma region Генерация граней со значением потока
void GenFacesFlowValue(); // генерация граней со значением потока
#pragma endregion

#pragma region Отчистить старое и создать новый портрет для матрицы
void ClearAndGenPortMatr(); // отчистить старое и создать новый портрет для матрицы
#pragma endregion

#pragma region Генерация матрицы и вектора
void GenMatrAndVec(); // генерация матрицы и вектора 
#pragma endregion

#pragma region Учет потоков
void ConsiderKnownFlows(int n); // учет известных потоков
void ConsiderFictFlows(int n); // учет фиктивных потоков
void ConsiderFlows(); // учет всех краевых
#pragma endregion

#pragma region Расчет небаланса и нахождение сбалансированных потоков
double CalcSumNonBalance(); // рассчитать суммарный небаланс
void FindBalancedFlows(); // найти сбалансированные потоки
#pragma endregion

#pragma region Балансировка потоков
void BalancingFlows(); // балансировка потоков
#pragma endregion
