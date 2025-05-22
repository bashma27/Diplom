#pragma once
#include "globals.h"

#pragma region Вспомогательные функции
double BasicFunc(int i, double input_point, double point_1, double point_2, double point_3); // используемые базисные функции
double DerivativeLocalBasicFuncs(int i, double psi); // производная по базисным функциям
bool IsFictitious(int num); // определение фиктивности узла
void GenFictNodes(); // генерация фиктивных узлов
void GenArrayFictEndEl(); // генерация массива фиктивных конечных элементов
bool IsFictEl(int num_end_el); // определение фиктивности конечного элемента
bool IsFindFaceZP(int num_face); // нахождение грани в массиве граней зон перфорации
bool IsFindFictFace(int num_face); // нахождение грани в массиве фиктивных граней
#pragma endregion

#pragma region Функции краевых условий, их функций и вектора правой части f
double u_g(double x, double y, double z); // краевое условие первого рода
double k_ph(int num_el, int num_ph); // множитель структурной проницаемости
double lambda(int num_sub); // расчет лямбды как итогового коэфиициента
double theta(int num_face_2_zp); // вычисление тетты по заданному потоку
double SpecFlow(int num_face_2_zp); // заданный поток на нужной грани для зоны перфорации
double f(int num_sub, double x, double y); // функция правой части диф. уравнения
#pragma endregion

#pragma region Генерация массива конечных элементов и портрета матрицы
void GenArrayParallelepipeds(); // генерация массива параллелепипедов
void GenPortraitMatr(); // генерация портрета матрицы
#pragma endregion

#pragma region Функции добавления в глоб. матрицу и вектор правой части
void AddLocalMatr(vector<int>& node_num, vector<vector<double>>& loc_matr); // добавление матрицы в глобальную
void AddLocalVec(vector<int>& node_num, vector<double>& loc_vec); // добавление вектора в глобальный
void AddLocalVecBound(int num_face_2_zp, vector<double>& vec); // добавление вектора из второго краевого в глобальный
#pragma endregion

#pragma region Работа с краевыми условиями (генерация и учет)
void GenFirstBoundCondit(); // создание массива с первыми краевыми условиями
void GenSecBoundCondit(); // создание массива с вторыми краевыми условиями
void ConsiderBoundConditFirstType(int n, vector<vector<pair<int, int>>>& columnToRows); // учет краевых условий первого типа
void ConsiderBoundConditSecType(int num_face_2_zp); // учет краевых условий второго типа
void ConsiderFictitiousNodes(); // учет фиктивных узлов
void ConsiderBoundCondit(); // учет всех краевых
#pragma endregion

#pragma region Работа с локальными матрицами и векторами
void CreateLocalMatrStiff(int num_sub, vector<int>& node_num); // локальная матрица жесткости
void CreateLocalVec(int num_sub, vector<int>& node_num); // локальный вектор правой части
#pragma endregion

#pragma region Построение глобальной матрицы и вектора
void BuildMatrVec(); // построить глобальную матрицу и вектор
#pragma endregion

#pragma region Тестирование для аналитически заданных функций
void Test(); // тест
#pragma endregion

#pragma region Возврат значения функции в любой точке
double GetResUInPoint(Coord3 point); // получить значение полученной функции в точке
#pragma endregion

#pragma region Вывод вектора результата решения СЛАУ в файл
void OutputResult(); // вывод результата в файл
#pragma endregion