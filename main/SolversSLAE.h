#pragma once

using namespace std;

// пространство имен для решателя СЛАУ МСГ
namespace MSG {
	double ScalarMult(vector<double>& first_vec, vector<double>& second_vec, int n);
	void MultMartVec(vector<double>& input_vec, vector<double>& res_vec, int n);
	void LU_sq(int n);
	void Ly_f(vector<double>& f, int n);
	void Ly_f_transp(vector<double>& f, int n);
	void PrMinusR(vector<double>& r, int n);
	void FirstEqualSecond(vector<double>& z, vector<double>& r, int n);
	double CulcResid(vector<double>& r, double norma_f, int n);
	void CulcX(vector<double>& x, vector<double>& z, double alpha, int n);
	void CulcR(vector<double>& r, vector<double>& Az, double alpha, int n);
	void CulcZ(vector<double>& Mult, vector<double>& z, double beta, int n);
	double CulcAlpha(vector<double>& Az, vector<double>& z, vector<double>& Mult, vector<double>& r, int n);
	double CulcBeta(double numerator, double denominator);
	void CulcMult(vector<double>& Mult, int n);
	void LU_sq_MSG(vector<double>& x, vector<double>& r, vector<double>& z, vector<double>& Az, vector<double>& Mult, int n, double eps, int max_iter);
}