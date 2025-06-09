#pragma once
#include "fem.h"

#pragma region Тестирование для аналитически заданных функций
void RelativeErrorNorm(); // относительная норма погрешности
vector<double> grad_u(double x, double y, double z); // градиент аналитической функции
vector<double> grad_uh(double xi, double eta, double zeta, int num_end_el); // градиент функции, полученой численно
void EnergyNorm(); // энергетическая норма
#pragma endregion

#pragma region Сравнение с аналитическим решением
double AnalitP(double x, double y, double z);
void VecAnalitP();
#pragma endregion