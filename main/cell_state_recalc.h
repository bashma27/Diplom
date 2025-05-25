#pragma once
#include "globals.h"
#include "get_funcs.h"
#include "fem.h"

#pragma region Расчет потоков для каждой фазы
void CalcFlowPh(); // расчитать потоки для каждоый фазы в отдельности
void RecalcFlowPh(double min_t); // регулировка потоков для минимального шага по времени
#pragma endregion

#pragma region Максимальынй шаг по времени
double MaxDeltaT(); //максимальный шаг по времени (учитывая каждую фазу, каждый конечный элемент)
double MaxDeltaTMixture(); //максимальный шаг по времени для смеси
#pragma endregion

#pragma region Расчет новых потоков и насыщенностей
void CalcNewFlowAndSatur(double h_t); // расчитать новые потоки на элементах и значения насыщенностей
#pragma endregion

#pragma region Подсчет суммарного объема фаз и порового объема
void CalcSumVPhAndPoreV(double max_h_t); // расчитать суммарный объем фаз и поровый объем
#pragma endregion

#pragma region Пересчет состояния ячеек
void RecalcCellState(); // пересчитать состояние ячеек
#pragma endregion