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
    NewtonPolynomial p;            // 调用你已有的牛顿多项式
    vector<double> x;              // 自己保存一份节点
    vector<double> coefficients;   // 自己保存一份系数

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
    int m = 2 * n;  // 每个点重复一次
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

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

int main() {
    // 使用和之前一样的函数 f(x)
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

    // 📌 1. 定义插值节点
    vector<double> x = {0.0, 0.5, 1.0};
    vector<double> y, z;
    for (double xi : x) {
        y.push_back(f(xi));
        z.push_back(f.diff(xi));
    }

    // 📌 2. 构造 Hermite 插值多项式
    auto p = hermiteInterp(x, y, z);

    // 📌 3. 插值点误差检测
    cout << fixed << setprecision(8);
    cout << "=== 插值点误差检测 ===" << endl;
    for (size_t i = 0; i < x.size(); i++) {
        double xi = x[i];
        double fx = f(xi);
        double dfx = f.diff(xi);
        double px = p(xi);
        double dpx = p.diff(xi);
        double err_f = fabs(fx - px);
        double err_df = fabs(dfx - dpx);
        cout << "x = " << setw(8) << xi
             << "  f(x) = " << setw(12) << fx
             << "  p(x) = " << setw(12) << px
             << "  |误差| = " << setw(12) << err_f << endl;
        cout << "          f'(x) = " << setw(12) << dfx
             << "  p'(x) = " << setw(12) << dpx
             << "  |误差| = " << setw(12) << err_df << endl;
    }

    return 0;
}