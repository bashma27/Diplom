#pragma once
#include "globals.h"
#include "get_funcs.h"
#include "fem.h"

#pragma region Расчет потоков для каждоый фазы
void CalcFlowPh(); // расчитать потоки для каждоый фазы в отдельности
#pragma endregion

#pragma region Максимальынй шаг по времени
double MaxDeltaT(); //максимальный шаг по времени (учитывая каждую фазу, каждый конечный элемент)
#pragma endregion

#pragma region Расчет новых потоков и насыщенностей
void CalcNewFlowAndSatur(double h_t); // расчитать новые потоки на элементах и значения насыщенностей
#pragma endregion

#pragma region Пересчет состояния ячеек
void RecalcCellState(); // пересчитать состояние ячеек
#pragma endregion