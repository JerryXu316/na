#include <vector>
#include <string>
using namespace std;

#ifndef SPLINE_WRAPPER
#define SPLINE_WRAPPER
#include <memory>
// SplineWrapper 用于将您返回的自定义 class 包装为标准接口，可以接收的返回类型包括：普通类对象、派生类对象、lambda 函数。
class SplineWrapper {
    struct Concept { virtual ~Concept() = default; virtual double eval(double) const = 0; };

    template<class T>
    struct Model final : Concept {
        T impl;
        explicit Model(T x) : impl(std::move(x)) {}
        double eval(double x) const override { return impl(x); }
    };

    std::shared_ptr<const Concept> self;

public:
    SplineWrapper() = default;

    template<class T>
    SplineWrapper(T&& x)
    : self(std::make_shared<Model<std::decay_t<T>>>(std::forward<T>(x))) {}

    double operator()(double x) const { return self->eval(x); }
    explicit operator bool() const { return static_cast<bool>(self); }
};
#endif



// 使用追赶法（Thomas算法）求解三对角线性方程组 Ax = b
vector<double> tridiagonalSolver(const vector<double>& diagonal,
                                    const vector<double>& upperDiagonal,
                                    const vector<double>& lowerDiagonal,
                                    const vector<double>& b) {
    int n = b.size();
    vector<double> x(n, 0.0);
    vector<double> c(n, 0.0);
    vector<double> d(n, 0.0);

    // 初始化c和d
    c[0] = upperDiagonal[0] / diagonal[0];
    d[0] = b[0] / diagonal[0];

    for (int i = 1; i < n; ++i) {
        double m = diagonal[i] - lowerDiagonal[i - 1] * c[i - 1];
        c[i] = upperDiagonal[i] / m;
        d[i] = (b[i] - lowerDiagonal[i - 1] * d[i - 1]) / m;
    }

    // 回代求解x
    x[n - 1] = d[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = d[i] - c[i] * x[i + 1];
    }

    return x;
}


// 循环三对角矩阵求解（适用于周期样条）
vector<double> cyclicTridiagonalSolver(const vector<double>& diag,
                                       const vector<double>& upper,
                                       const vector<double>& lower,
                                       const vector<double>& b) {
    int n = b.size();

    // Step 1: 构造 u, v
    vector<double> u(n, 0.0), v(n, 0.0);
    u[0] = lower.back();    // A(0, n-1)
    u.back() = upper.back(); // A(n-1, 0)
    v[0] = 1.0;
    v.back() = 1.0;

    // Step 2: 去掉首尾连接，形成普通三对角矩阵
    vector<double> diag_mod = diag;
    diag_mod[0] -= u[0];      // 调整主对角以保证稳定性
    diag_mod.back() -= u.back();

    // Step 3: 两次求解 T*y=b 和 T*z=u
    vector<double> y = tridiagonalSolver(diag_mod, upper, lower, b);
    vector<double> z = tridiagonalSolver(diag_mod, upper, lower, u);

    // Step 4: Sherman–Morrison 修正
    double vTy = y.front() + y.back();
    double vTz = z.front() + z.back();
    double factor = vTy / (1.0 + vTz);

    vector<double> x(n);
    for (int i = 0; i < n; ++i)
        x[i] = y[i] - factor * z[i];

    return x;
}


SplineWrapper getCubicSpline(
    const vector<double> &t, 
    const vector<double> &f, 
    const string &bdType, 
    const vector<double> &bdVal = {}) {
    
    // 请将您自定义的样条对象（需要支持小括号运算）直接作为此函数的返回值，SplineWrapper 会自动将您的返回对象包装成标准接口
    int n = t.size();
    vector<double> u(n - 2, 0.0);
    vector<double> r(n - 2, 0.0);
    for(int i = 0; i < n - 2; ++i){
        u[i] = (t[i + 1] - t[i]) / (t[i + 2] - t[i]);
        r[i] = (t[i + 2] - t[i + 1]) / (t[i + 2] - t[i]);
    }

    // 一阶差商f[xi,xi+1]计算
    vector<double> first_difference(n - 1, 0.0);
    for(int i = 0; i < n - 1; ++i){
        first_difference[i] = (f[i + 1] - f[i]) / (t[i + 1] - t[i]);
    }
    
    // 二阶差商f[xi,xi+1,xi+2]计算
    vector<double> second_difference(n - 2, 0.0);
    for(int i = 0; i < n - 2; ++i){
        second_difference[i] = (first_difference[i + 1] - first_difference[i]) / (t[i + 2] - t[i]);
    }

    // 插值点的二阶导数
    vector<double> M(n, 0.0);

    // 分情况计算一二阶导数
    if (bdType == "Complete") {
        // 这是课本上已经证明了的情况，我们只需求解三对角阵的线性方程组即可
        // 求解二阶导数Mi
        vector<double> M_b(n, 0.0);
        M_b[0] = 6 * (first_difference[0] - bdVal[0]) / (t[1] - t[0]);
        M_b[n - 1] = 6 * (bdVal[1] - first_difference[n - 2]) / (t[n - 1] - t[n - 2]);
        for(int i = 1; i < n - 1; ++i){
            M_b[i] = 6 * second_difference[i - 1];
        }

        //对角线元素
        vector<double> M_diag(n, 2.0);

        //上对角线元素        
        vector<double> M_upperdiag(n - 1, 0.0);
        M_upperdiag[0] = 1.0;
        for(int i = 1; i < n - 1; ++i){
            M_upperdiag[i] = r[i - 1];
        }


        //下对角线元素
        vector<double> M_lowerdiag(n - 1, 0.0);
        M_lowerdiag[n - 2] = 1.0;
        for(int i = 0; i < n - 2; ++i){
            M_upperdiag[i] = u[i];
        }

        M = tridiagonalSolver(M_diag, M_upperdiag, M_lowerdiag, M_b);
    } else if (bdType == "Natural") {
        //这种情况就是上面的简化，因为首尾都知道了
        vector<double> M_b(n, 0.0);
        M_b[0] = 0.0;
        M_b[n - 1] = 0.0;
        for(int i = 1; i < n - 1; ++i){
            M_b[i] = 6 * second_difference[i - 1];
        }

        //对角线元素
        vector<double> M_diag(n, 2.0);
        M_diag[0] = 1.0;
        M_diag[n - 1] = 1.0;

        //上对角线元素
        vector<double> M_upperdiag(n - 1, 0.0);
        M_upperdiag[0] = 0.0;
        for(int i = 1; i < n - 1; ++i){
            M_upperdiag[i] = r[i - 1];
        }


        //下对角线元素        
        vector<double> M_lowerdiag(n - 1, 0.0);
        M_lowerdiag[n - 2] = 0.0;
        for(int i = 0; i < n - 2; ++i){
            M_upperdiag[i] = u[i];
        }

        M = tridiagonalSolver(M_diag, M_upperdiag, M_lowerdiag, M_b);


    } else if (bdType == "Periodic") {
        //这里的情况略有不同，但是我们还是一样的处理方法，再代入消去xn(这里变量命名暂且与课本理论部分一致)
        //我们有关系式2 * (x2 + xn - x1 - xn-1) * M1 + (x2 - x1) * M2 + (xn - xn-1) * Mn-1 = 6f[x1, x2] - 6f[xn-1, xn]
        //这时矩阵并不是一个严格的三对角阵，而是在右上角与左下角分别有一个元素，我们考虑使用Sherman–Morrison 修正法进行求解
        vector<double> M_b(n - 1, 0.0);
        M_b[0] = 6 * first_difference[0] - 6 * first_difference[n - 2];
        for(int i = 1; i < n - 1; ++i){
            M_b[i] = 6 * second_difference[i - 1];
        }

        //对角线元素
        vector<double> M_diag(n - 1, 2.0);
        M_diag[0] = 2.0 * (t[1] + t[n - 1] - t[0] - t[n - 2]);


        //上对角线元素
        vector<double> M_upperdiag(n - 1, 0.0);
        M_upperdiag[0] = t[1] - t[0];
        for(int i = 1; i < n - 1; ++i){
            M_upperdiag[i] = r[i - 1];
        }


        //下对角线元素        
        vector<double> M_lowerdiag(n - 1, 0.0);
        M_lowerdiag[n - 2] = t[n - 1] - t[n - 2];
        for(int i = 0; i < n - 2; ++i){
            M_upperdiag[i] = u[i];
        }


        vector<double> x(n - 1, 0.0);
        x = cyclicTridiagonalSolver(M_diag, M_upperdiag, M_lowerdiag, M_b);
        for(int i = 0; i < n - 1; ++i){
            M[i] = x[i];
        }
        M[n - 1] = x[0];


    } else if (bdType == "Second-Derivatives") {
        //这与Natural的情况是一致的，只要把0.0改为对应的值
        vector<double> M_b(n, 0.0);
        M_b[0] = bdVal[0];
        M_b[n - 1] = bdVal[1];
        for(int i = 1; i < n - 1; ++i){
            M_b[i] = 6 * second_difference[i - 1];
        }

        //对角线元素        
        vector<double> M_diag(n, 2.0);
        M_diag[0] = 1.0;
        M_diag[n - 1] = 1.0;


        vector<double> M_upperdiag(n - 1, 0.0);
        M_upperdiag[0] = 0.0;
        for(int i = 1; i < n - 1; ++i){
            M_upperdiag[i] = r[i - 1];
        }

        vector<double> M_lowerdiag(n - 1, 0.0);
        M_lowerdiag[n - 2] = 0.0;
        for(int i = 0; i < n - 2; ++i){
            M_upperdiag[i] = u[i];
        }

        M = tridiagonalSolver(M_diag, M_upperdiag, M_lowerdiag, M_b);
    } else if (bdType == "Not-A-Knot") {


    }





}