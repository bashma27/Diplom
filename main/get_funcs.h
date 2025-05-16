#pragma once
#include "globals.h"

int GetNumSubarea(vector<int>& end_el); // получить номер подобласти для конечного элемента
int GetNumNodes(vector<Coord2>& zone_perf_i, int i_x, int i_y); // получить номер узла для зоны перфорации
int GetNumSubareaBound(int ind, vector<int>& face_2_zp); // получить номер подобласти для второго краевого (первый индекс 0/1 - грань параллельная оси y/x)
int GetNumEndEl(Coord3 point); // получить номер конечного элемента для точки
int GetNumFace2ZP(vector<int>& _face_2_zp); // получить номер грани зоны перфорации
