#include <iostream>
#include <vector>
using namespace std;

// 在函数外定义类
class NewtonPolynomial {
private:
    vector<double> x;              // 节点
    vector<double> coefficients;   // 差商系数

public:
    NewtonPolynomial(const vector<double> &x_, const vector<double> &coef_)
        : x(x_), coefficients(coef_) {}

    // 让对象能用 p(x0) 方式调用
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

// 不能改函数格式
auto newtonInterp(const vector<double> &x,
                  const vector<double> &y)
{
    int n = x.size();
    vector<double> coef = y; // 初始为y值

    // 计算差商
    for (int i = 1; i < n; i++) {
        for (int j = n - 1; j >= i; j--) {
            coef[j] = (coef[j] - coef[j - 1]) / (x[j] - x[j - i]);
        }
    }

    // 返回多项式对象
    return NewtonPolynomial(x, coef);
}

int main() {
    vector<double> x = {1.0, 2.0, 3.0};
    vector<double> y = {2.0, 3.0, 5.0};

    auto p = newtonInterp(x, y);  // 使用固定函数签名

    double x0 = 2.5;
    cout << "P(" << x0 << ") = " << p(x0) << endl;

    return 0;
}
