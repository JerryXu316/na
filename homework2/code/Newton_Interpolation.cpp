#include <vector>
using namespace std;

class NewtonPolynomial {
    vector<double> coeff; // 差商系数 a0, a1, ...
    vector<double> x;     // 节点

public:
    NewtonPolynomial(const vector<double>& xs, const vector<double>& ys) {
        int n = xs.size();
        x = xs;
        coeff = ys; // 初始差商表的第一列是 y 值

        // 计算差商，原地更新
        for (int j = 1; j < n; j++) {
            for (int i = n - 1; i >= j; i--) {
                coeff[i] = (coeff[i] - coeff[i - 1]) / (xs[i] - xs[i - j]);
            }
        }
        // coeff[0..n-1] 就是 a0..an
    }

    double operator()(double x0) const {
        double res = 0;
        int n = coeff.size();
        // Newton 公式展开
        for (int i = n - 1; i >= 0; i--) {
            res = res * (x0 - x[i]) + coeff[i];
        }
        return res;
    }
};

auto newtonInterp(const vector<double> &x, const vector<double> &y) {
    return NewtonPolynomial(x, y);
}


#include <iostream>

int main() {
    vector<double> x = {0, 1, 2};
    vector<double> y = {1, 3, 2};
    auto p = newtonInterp(x, y);

    cout << p(0) << endl; // 应该等于 1
    cout << p(1) << endl; // 应该等于 3
    cout << p(2) << endl; // 应该等于 2
    cout << p(1.5) << endl; // 插值结果
}