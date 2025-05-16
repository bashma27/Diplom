#include "solver_slae.h"

namespace MSG {
    double ScalarMult(vector<double>& first_vec, vector<double>& second_vec, int n) {
        double res = 0;
        for (int i = 0; i < n; i++) {
            res += first_vec[i] * second_vec[i];
        }
        return res;
    }

    void MultMartVec(vector<double>& input_vec, vector<double>& res_vec, int n) {
        for (int i = 0; i < n; i++) {
            res_vec[i] = di[i] * input_vec[i];
            for (int k = ia[i]; k < ia[i + 1]; k++) {
                int j = ja[k];
                res_vec[i] += aal[k] * input_vec[j];
                res_vec[j] += aal[k] * input_vec[i];
            }
        }
    }

    void LU_sq(int n) {
        for (int i = 0; i < n; i++) {
            double s_d = 0;
            for (int k = ia[i]; k < ia[i + 1]; k++) {
                double s_l = 0;
                int j = ja[k];
                int j0 = ia[j];
                int j1 = ia[j + 1];
                int ki = ia[i];
                int kj = j0;
                for (; ki < k && kj < j1;) {
                    int jl = ja[ki];
                    int ju = ja[kj];
                    if (jl == ju) {
                        s_l += L_sq[kj] * L_sq[ki];
                        ki++, kj++;
                    }
                    else if (jl < ju) ki++;
                    else kj++;
                }
                L_sq[k] = (aal[k] - s_l) / di_sq[j];
                s_d += L_sq[k] * L_sq[k];
            }
            di_sq[i] = sqrt(di[i] - s_d);
        }
    }

    void Ly_f(vector<double>& f, int n) {
        for (int i = 0; i < n; i++) {
            double y = 0;
            for (int k = ia[i]; k < ia[i + 1]; k++) {
                int j = ja[k];
                y += L_sq[k] * f[j];
            }
            f[i] = (f[i] - y) / di_sq[i];
        }
    }

    void Ly_f_transp(vector<double>& f, int n) {
        for (int i = n - 1; i >= 0; i--) {
            f[i] /= di_sq[i];
            for (int k = ia[i]; k < ia[i + 1]; k++) {
                int j = ja[k];
                f[j] -= L_sq[k] * f[i];
            }
        }
    }

    void PrMinusR(vector<double>& r, int n) {
        for (int i = 0; i < n; i++) {
            r[i] = b[i] - r[i];
        }
    }

    void FirstEqualSecond(vector<double>& z, vector<double>& r, int n) {
        for (int i = 0; i < n; i++) {
            z[i] = r[i];
        }
    }

    double CulcResid(vector<double>& r, double norma_f, int n) {
        return sqrt(ScalarMult(r, r, n)) / norma_f;
    }

    void CulcX(vector<double>& x, vector<double>& z, double alpha, int n) {
        for (int i = 0; i < n; i++) {
            x[i] += alpha * z[i];
        }
    }

    void CulcR(vector<double>& r, vector<double>& Az, double alpha, int n) {
        for (int i = 0; i < n; i++) {
            r[i] -= alpha * Az[i];
        }
    }

    void CulcZ(vector<double>& Mult, vector<double>& z, double beta, int n) {
        for (int i = 0; i < n; i++) {
            z[i] = Mult[i] + beta * z[i];
        }
    }

    double CulcAlpha(vector<double>& Az, vector<double>& z, vector<double>& Mult, vector<double>& r, int n) {
        return ScalarMult(Mult, r, n) / ScalarMult(Az, z, n);
    }

    double CulcBeta(double numerator, double denominator) {
        return numerator / denominator;
    }

    void CulcMult(vector<double>& Mult, int n) {
        Ly_f(Mult, n);
        Ly_f_transp(Mult, n);
    }

    void LU_sq_MSG(vector<double>& x, vector<double>& r, vector<double>& z, vector<double>& Az, vector<double>& Mult, int n, double eps, int max_iter) {
        double resid, norma_f;
        int counter = 1;
        LU_sq(n);
        MultMartVec(x, r, n);
        PrMinusR(r, n);
        norma_f = sqrt(ScalarMult(b, b, n));
        resid = CulcResid(r, norma_f, n);
        cout << "Start resid: " << resid << endl;
        FirstEqualSecond(Mult, r, n);
        CulcMult(Mult, n);
        FirstEqualSecond(z, Mult, n);
        for (; counter < max_iter && resid > eps; counter++) {
            double alpha, beta, numerator, denominator;
            resid = 0;
            MultMartVec(z, Az, n);
            alpha = CulcAlpha(Az, z, Mult, r, n);
            denominator = ScalarMult(Mult, r, n);
            CulcX(x, z, alpha, n);
            CulcR(r, Az, alpha, n);
            FirstEqualSecond(Mult, r, n);
            CulcMult(Mult, n);
            numerator = ScalarMult(Mult, r, n);
            beta = CulcBeta(numerator, denominator);
            CulcZ(Mult, z, beta, n);
            resid = sqrt(ScalarMult(r, r, n)) / norma_f;
            cout << setprecision(15) << resid << " " << counter << endl;
        }
        MultMartVec(x, r, n);
        PrMinusR(r, n);
        resid = CulcResid(r, norma_f, n);
        cout << "End resid: " << resid << endl;
    }
}