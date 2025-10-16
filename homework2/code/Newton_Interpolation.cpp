#include <iostream>
#include <vector> 
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



auto newtonInterp(const vector<double> &x, const vector<double> &y) { 
    int n = x.size(); 
    vector<double> coefficients = y; 
    for(int i = 0; i < n - 1; i++){ 
        for(int j = 0; j < n - i - 1; j++){ 
            coefficients[n - j - 1] = (coefficients[n - j - 1] - coefficients[n - j - 2]) / (x[n - j - 1] - x[n - j - i - 2]);
        } 
    } 

    return NewtonPolynomial(x, coefficients);

} 
        
        
/*
auto p = newtonInterp(x, y);
// 你需要确保评测系统可以用这种方式使用你的返回对象：通过小括号运算，获取你的插值多项式在 x0 处的取值。
double f0 = p(x0);
*/

int main() {
    vector<double> x = {1.0, 2.0, 3.0};
    vector<double> y = {2.0, 3.0, 5.0};

    try {
        auto p = newtonInterp(x, y);
        double x0 = 2.5;
        cout << "P(" << x0 << ") = " << p(x0) << endl; // 4.25
    } catch (const exception &e) {
        cerr << "Error: " << e.what() << endl;
    }

    return 0;
}