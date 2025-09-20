#include <iostream>
#include <cmath>
using namespace std;


class Function {
public:
    virtual double operator () (double x) = 0;
    virtual double d(double x) = 0;
};

class NewtonSolver {
protected:
    Function &F;
public:
    NewtonSolver(Function & F) : F(F) {}
    double solve(double x0) {
        const int MAX_ITERATION = 1000;
        const double EPS = 1e-9;

        double x = x0;

        for (int i = 0; i < MAX_ITERATION ; i++){
            double fx = F(x);
            double dfx = F.d(x);

            if (fabs(dfx) < 1e-12){
                cout << "导数太小，牛顿法失败！" << endl;
                return x;
            }

            double x_next = x - fx / dfx;

            if (fabs(x_next - x) < EPS){
                return x_next;
            } 

            x = x_next;
        } 
        
        cout << "迭代未收敛，返回最后结果。" << endl;
        return x;

    }
};


// 你可以使用下面的代码在自己电脑上测试正确性，但提交时请勿包含以下代码
class Func1 : public Function {
public:
    double operator () (double x) {
        return x - tan(x);
    }
    double d(double x) {
        return 1 - 1 / pow(cos(x), 2);
    }
};

int main() {
    Func1 func;
    NewtonSolver solver(func);
    double x = solver.solve(4.5);
    cout << x << endl;
    return 0;
}
