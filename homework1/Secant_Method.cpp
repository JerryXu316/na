#include <iostream>
#include <cmath>


class Function {
public:
    virtual double operator () (double x) = 0;
};

class SecantSolver {
protected:
    Function &F;
public:
    SecantSolver(Function & F) : F(F) {}
    double solve(double x0, double x1) {
        // 在这里实现你的割线法，将答案返回
        double a = x1;
        double b = x0;

        double u = F(a);
        double v = F(b);


        double s = 0;

        for (int i = 0; i < 28; ++i){
            s = (a - b) / (u - v);
            b = a;
            v = u;
            a = a - u * s;
            if (std::fabs(a - b) < 1e-10){
                return a;
            }
            u = F(a);
            if (std::fabs(u) < 1e-12){
                return a;
            }
        }

        return a;
    }
};


class FuncSqrt2 : public Function {
public:
    double operator()(double x) {
        return x * x - 2.0;
    }
};

int main() {
    FuncSqrt2 f;
    SecantSolver solver(f);
    double root = solver.solve(1.0, 2.0); // 初值选 [1,2]
    std::cout.precision(15);
    std::cout << "Root ≈ " << root << std::endl;
    std::cout << "Error = " << fabs(root - sqrt(2.0)) << std::endl;
    return 0;
}
