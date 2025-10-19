#include <iostream>
#include <vector>
#include <cmath>
using namespace std;


class NewtonPolynomial {
private:
    vector<double> x;
    vector<double> coefficients;

public:
    NewtonPolynomial(const vector<double> &x_, const vector<double> &coeff_)
        : x(x_), coefficients(coeff_) {}

    double operator()(double x0) const {
        double result = 0.0;
        double term = 1.0;
        for (size_t i = 0; i < coefficients.size(); i++) {
            result += coefficients[i] * term;
            term *= (x0 - x[i]);
        }
        return result;
    }
};



class HermitePolynomial {
private:
    NewtonPolynomial p;            
    vector<double> x;              
    vector<double> coefficients;   

public:
    HermitePolynomial(const NewtonPolynomial &poly,
                      const vector<double> &x_,
                      const vector<double> &coeff_)
        : p(poly), x(x_), coefficients(coeff_) {}

    double operator()(double x0) const {
        return p(x0);
    }

    double diff(double x0) const {
        double result = 0.0;
        double term = 1.0;
        double term_deriv = 0.0;

        for (size_t i = 0; i < coefficients.size(); i++) {
            result += coefficients[i] * term_deriv;
            term_deriv = term_deriv * (x0 - x[i]) + term;
            term *= (x0 - x[i]);
        }
        return result;
    }
};


auto hermiteInterp(const vector<double> &x,
                   const vector<double> &y,
                   const vector<double> &z)
{
    int n = x.size();
    int m = 2 * n;  
    vector<double> X(m);
    vector<double> Y(m);


    for (int i = 0; i < n; i++) {
        X[2 * i] = x[i];
        X[2 * i + 1] = x[i];
        Y[2 * i] = y[i];
        Y[2 * i + 1] = y[i];
    }

    vector<double> coeff = Y;
    for (int k = 1; k < m; k++) {
        for (int i = m - 1; i >= k; i--) {
            if (fabs(X[i] - X[i - k]) < 1e-14) {
                coeff[i] = z[i / 2];
            } else {
                coeff[i] = (coeff[i] - coeff[i - 1]) / (X[i] - X[i - k]);
            }
        }
    }

    NewtonPolynomial newtonP(X, coeff);

    return HermitePolynomial(newtonP, X, coeff);

}


// 下面是评测系统使用的测试函数
class Func{
public:
    double operator () (const double &x) const{
        double tt = (x == 0) ? 0 : pow(x, 1.1) * log(x);
        return sin(M_PI * x) + sin(M_PI * 4 * (x - 0.2)) + tt;
    }
    double diff (const double &x) const{
        double tt = (x == 0) ? 0 : 1.1 * pow(x,0.1) * log(x);
        return M_PI * cos(M_PI * x) + 4 * M_PI * cos(M_PI * 4 * (x - 0.2)) + (tt + pow(x,0.1));
    }
} f;


/*
auto p = hermiteInterp(x, y, z);
// 你需要确保评测系统可以用这种方式使用你的返回对象：通过小括号运算，获取你的插值多项式在 x0 处的取值。
double f0 = p(x0);
// 通过 diff 函数，获取你的插值多项式在 x0 处的导数值。
double df0 = p.diff(x0);
*/

int main() {
    vector<double> x = {0.0, 0.5, 1.0};
    vector<double> y = {f(0.0), f(0.5), f(1.0)};
    vector<double> z = {f.diff(0.0), f.diff(0.5), f.diff(1.0)};

    auto p = hermiteInterp(x, y, z);

    cout << fixed;
    for (double x0 = 0.0; x0 <= 1.0; x0 += 0.2) {
        double f_val = f(x0);
        double p_val = p(x0);
        double f_diff = f.diff(x0);
        double p_diff = p.diff(x0);
        cout << "x = " << x0 
             << "  f(x) = " << f_val 
             << "  p(x) = " << p_val 
             << "  |误差| = " << fabs(f_val - p_val) << endl;
        cout << "       f'(x) = " << f_diff 
             << "  p'(x) = " << p_diff 
             << "  |误差| = " << fabs(f_diff - p_diff) << endl;
    }

    return 0;
}