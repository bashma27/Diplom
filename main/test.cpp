#include "test.h"

#pragma region Тестирование для аналитически заданных функций
void RelativeErrorNorm() { // относительная норма погрешности
    vector<double> q_u(NUM_NODES, 0);
    for (int i = 0; i < NUM_NODES; i++) {
        q_u[i] = u_g(nodes[i].x, nodes[i].y, nodes[i].z);
    }
    double norm_vec_err = 0, norm_vec_q_u = 0; // норма вектора погрешности и q_u
    for (int i = 0; i < NUM_NODES; i++) {
        if (IsFictitious(i)) continue;
        norm_vec_err += (q[i] - q_u[i]) * (q[i] - q_u[i]);
        norm_vec_q_u += (q_u[i]) * (q_u[i]);
    }
    cout << endl;
    cout << "Относительная норма вектора погрешности полученного решения:" << endl;
    cout << sqrt(norm_vec_err) / sqrt(norm_vec_q_u) << endl << endl;
}

const double gauss_points[4] = {
    0.5 * (-0.8611363116 + 1.),
    0.5 * (-0.3399810436 + 1.),
    0.5 * (0.3399810436 + 1.),
    0.5 * (0.8611363116 + 1.)
};

const double gauss_weights[4] = {
    0.5 * 0.3478548451,
    0.5 * 0.6521451549,
    0.5 * 0.6521451549,
    0.5 * 0.3478548451
};

vector<double> grad_u(double x, double y, double z) { // градиент аналитической функции
    return {-theta_ * 1. / lambda(0) * x / (x * x + y * y), -theta_ * 1. / lambda(0) * y / (x * x + y * y), 0.};
}

vector<double> grad_uh(double xi, double eta, double zeta, int num_end_el) { // градиент функции, полученой численно
    double gx = 0., gy = 0., gz = 0.;
    double h_x, h_y, h_z;
    h_x = nodes[array_p[num_end_el].second[2]].x - nodes[array_p[num_end_el].second[0]].x;
    h_y = nodes[array_p[num_end_el].second[6]].y - nodes[array_p[num_end_el].second[0]].y;
    h_z = nodes[array_p[num_end_el].second[18]].z - nodes[array_p[num_end_el].second[0]].z;

    for (int i = 0; i < array_p[num_end_el].second.size(); i++) {
        
        if (IsFictEl(i)) {
            continue;
        }

        int ind = array_p[num_end_el].second[i];

        gx += q[ind] * DerivativeLocalBasicFuncs(i % 3, xi) * LocalBasicFuncs((i / 3) % 3, eta) * LocalBasicFuncs(i / 9, zeta);
        gy += q[ind] * LocalBasicFuncs(i % 3, xi) * DerivativeLocalBasicFuncs((i / 3) % 3, eta) * LocalBasicFuncs(i / 9, zeta);
        gz += q[ind] * LocalBasicFuncs(i % 3, xi) * LocalBasicFuncs((i / 3) % 3, eta) * DerivativeLocalBasicFuncs(i / 9, zeta);
    }
    
    return { gx / h_x, gy / h_y, gz / h_z };
}

void EnergyNorm() { // энергетическая норма
    double norm = 0.;

    for (int num_end_el = 0; num_end_el < NUM_SPLIT_X * NUM_SPLIT_Y * NUM_SPLIT_Z; num_end_el++) {

        if (IsFictEl(num_end_el)) {
            continue;
        }

        double sum = 0.;
        double h_x, h_y, h_z;
        h_x = nodes[array_p[num_end_el].second[2]].x - nodes[array_p[num_end_el].second[0]].x;
        h_y = nodes[array_p[num_end_el].second[6]].y - nodes[array_p[num_end_el].second[0]].y;
        h_z = nodes[array_p[num_end_el].second[18]].z - nodes[array_p[num_end_el].second[0]].z;

        for (int i = 0; i < 4; i++) {
            double xi = gauss_points[i];
            double wi = gauss_weights[i];

            for (int j = 0; j < 4; j++) {
                double eta = gauss_points[j];
                double wj = gauss_weights[j];

                for (int k = 0; k < 4; k++) {
                    double zeta = gauss_points[k];
                    double wk = gauss_weights[k];

                    // Переводим в физические координаты
                    double x = nodes[array_p[num_end_el].second[0]].x + xi * h_x;
                    double y = nodes[array_p[num_end_el].second[0]].y + eta * h_y;
                    double z = nodes[array_p[num_end_el].second[0]].z + zeta * h_z;

                    auto grad_exact = grad_u(x, y, z);
                    auto grad_approx = grad_uh(xi, eta, zeta, num_end_el);

                    double ex = grad_exact[0] - grad_approx[0];
                    double ey = grad_exact[1] - grad_approx[1];
                    double ez = grad_exact[2] - grad_approx[2];

                    double integrand = ex * ex + ey * ey + ez * ez;
                    sum += wi * wj * wk * integrand;
                }
            }
        }
        norm += lambda(array_p[num_end_el].first) * sum * h_x * h_y * h_z;
    }
    cout << "Энергетическая норма ошибки: " << sqrt(norm) << endl;
}
#pragma endregion

#pragma region Сравнение с аналитическим решением
double AnalitP(double x, double y, double z) {   
    double res = theta_ * 1. / lambda(0) * log(100. / sqrt(x * x + y * y)) + 130;
    return res;
}

void VecAnalitP() {
    ofstream f1("relative_error.txt");
    ofstream f2("r.txt");
    double r = 100.;
    int start = NUM_NODES_IN_EDGE_X * (2. * coef[1][0].first + 1);
    int end = NUM_NODES_IN_EDGE_X * (2. * coef[1][0].first + 1) + (2. * coef[0][0].first + 1);
    for (int num_nodes = start; num_nodes < end; num_nodes++) {
        double analyt_value_P = AnalitP(nodes[num_nodes].x, nodes[num_nodes].y, nodes[num_nodes].z);
        double num_value_P = q[num_nodes];
        double r = -nodes[num_nodes].x - 1;
        double relative_error = fabs(num_value_P - analyt_value_P) / fabs(analyt_value_P);
        f1 << relative_error << endl;
        f2 << r << endl;
    }
    f1.close();
    f2.close();
}

//void VecAnalitP() {
//    double norm = 0.;
//    for (int num_end_el = NUM_SPLIT_X * coef[1][0].first; num_end_el < NUM_SPLIT_X * coef[1][0].first + coef[0][0].first; num_end_el++) {
//        double sum = 0.;
//        double h_x, h_y, h_z;
//        h_x = nodes[array_p[num_end_el].second[2]].x - nodes[array_p[num_end_el].second[0]].x;
//        h_y = nodes[array_p[num_end_el].second[6]].y - nodes[array_p[num_end_el].second[0]].y;
//        h_z = nodes[array_p[num_end_el].second[18]].z - nodes[array_p[num_end_el].second[0]].z;
//
//        for (int i = 0; i < 4; i++) {
//            double xi = gauss_points[i];
//            double wi = gauss_weights[i];
//
//            for (int j = 0; j < 4; j++) {
//                double eta = gauss_points[j];
//                double wj = gauss_weights[j];
//
//                for (int k = 0; k < 4; k++) {
//                    double zeta = gauss_points[k];
//                    double wk = gauss_weights[k];
//
//                    // Переводим в физические координаты
//                    double x = nodes[array_p[num_end_el].second[0]].x + xi * h_x;
//                    double y = nodes[array_p[num_end_el].second[0]].y + eta * h_y;
//                    double z = nodes[array_p[num_end_el].second[0]].z + zeta * h_z;
//
//                    auto grad_exact = grad_u(x, y, z);
//                    auto grad_approx = grad_uh(xi, eta, zeta, num_end_el);
//
//                    double ex = grad_exact[0] - grad_approx[0];
//                    double ey = grad_exact[1] - grad_approx[1];
//                    double ez = grad_exact[2] - grad_approx[2];
//
//                    double integrand = ex * ex + ey * ey + ez * ez;
//                    sum += wi * wj * wk * integrand;
//                }
//            }
//        }
//        norm += lambda(array_p[num_end_el].first) * sum * h_x * h_y * h_z;
//    }
//    cout << "Энергетическая норма ошибки (аналитическое решение): " << sqrt(norm) << endl;
//}
#pragma endregion