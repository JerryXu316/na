#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

class Func{
public:
    double operator () (const double &x) const{
        return 1.0 / (1.0 + 25.0 * x * x);
    }
} g;



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

auto newtonInterp(const vector<double> &x, const vector<double> &y) { 
    int n = x.size(); 
    vector<double> coefficients = y; 
    for (int k = 1; k < n; k++) {
        for (int i = n - 1; i >= k; i--) {
            coefficients[i] = (coefficients[i] - coefficients[i - 1]) / (x[i] - x[i - k]);
        }
    }
    return NewtonPolynomial(x, coefficients);

} 

auto interpRungeFunction(int n) {
    vector<double> x(n + 1);
    vector<double> y(n + 1);

    for(int k = 0; k < n + 1; k++){
        x[k] = cos((k + 1.0) / (n + 1) * M_PI);
        y[k] = g(x[k]);
    }

    auto p = newtonInterp(x, y);

    return p;
    
}

/*
auto p = interpRungeFunction(n);
// 你需要确保评测系统可以用这种方式使用你的返回对象：通过小括号运算，获取你的插值多项式在 x0 处的取值。
double p0 = p(x0);
*/
// 测试主函数
int main() {
    int n = 10; // 插值节点数
    auto p = interpRungeFunction(n);

    cout << fixed << setprecision(6);
    cout << " x\t真实值 f(x)\t插值值 p(x)\t误差\n";
    cout << "--------------------------------------------\n";

    for (double x0 = -1.0; x0 <= 1.0; x0 += 0.2) {
        double fx = g(x0);
        double px = p(x0);
        double err = fabs(fx - px);
        cout << setw(5) << x0 << "\t" << setw(10) << fx
             << "\t" << setw(10) << px << "\t" << setw(10) << err << "\n";
    }

    return 0;
}