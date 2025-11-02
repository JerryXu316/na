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

class PPform {
private:
    vector<double> x; // 节点
    vector<double> f; // 函数值
    vector<double> M; // 二阶导数

    vector<double> a, b, c, d; // 每段系数

    void computeCoefficients() {
        size_t n = x.size();
        a.resize(n - 1);
        b.resize(n - 1);
        c.resize(n - 1);
        d.resize(n - 1);

        for (size_t i = 0; i < n - 1; ++i) {
            double h = x[i + 1] - x[i];
            a[i] = f[i];
            b[i] = (f[i + 1] - f[i]) / h - h * (2 * M[i] + M[i + 1]) / 6.0;
            c[i] = M[i] / 2.0;
            d[i] = (M[i + 1] - M[i]) / (6.0 * h);
        }
    }

    // 找到x0所在的区间
    size_t findIndex(double x0) const {
        for (size_t i = 0; i < x.size() - 1; ++i)
            if (x0 >= x[i] && x0 <= x[i + 1])
                return i;

        return x.size() - 2; // 最后一个区间
    }

public:
    PPform(const vector<double> &x_, const vector<double> &f_, const vector<double> &M_)
        : x(x_), f(f_), M(M_) {
        computeCoefficients();
    }

    // 计算样条插值值
    double operator()(double x0) const {
        size_t i = findIndex(x0);
        double dx = x0 - x[i];
        return a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
    }
};



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



// Not-a-Knot 样条的线性系统求解
vector<double> solveDenseLinearSystem(vector<vector<double>> A, vector<double> b) {
    int n = (int)b.size();
    const double EPS = 1e-12;
    for (int i = 0; i < n; ++i) {
        // pivot
        int piv = i;
        for (int j = i+1; j < n; ++j)
            if (fabs(A[j][i]) > fabs(A[piv][i])) piv = j;
        if (fabs(A[piv][i]) < EPS) {
            // singular or nearly singular; continue (will produce unstable result)
            continue;
        }
        if (piv != i) {
            swap(A[piv], A[i]);
            swap(b[piv], b[i]);
        }
        // normalize and eliminate
        double diag = A[i][i];
        for (int k = i; k < n; ++k) A[i][k] /= diag;
        b[i] /= diag;
        for (int r = 0; r < n; ++r) {
            if (r == i) continue;
            double factor = A[r][i];
            if (fabs(factor) < EPS) continue;
            for (int c = i; c < n; ++c) A[r][c] -= factor * A[i][c];
            b[r] -= factor * b[i];
        }
    }
    return b; 
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
            M_lowerdiag[i] = u[i];
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
            M_lowerdiag[i] = u[i];
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
            M_lowerdiag[i] = u[i];
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
            M_lowerdiag[i] = u[i];
        }

        M = tridiagonalSolver(M_diag, M_upperdiag, M_lowerdiag, M_b);
    }else if (bdType == "Not-A-Knot") {
        // 构造 h_i = t[i+1]-t[i]
        vector<double> h(n-1);
        for (int i = 0; i < n-1; ++i) h[i] = t[i+1] - t[i];

        // 构造密矩阵 A (n x n) 和右端 b (size n)
        vector<vector<double>> A(n, vector<double>(n, 0.0));
        vector<double> bb(n, 0.0);

        // 首行 (Not-a-knot at x1): -h1*M0 + (h0+h1)*M1 - h0*M2 = 0
        // 对应列 0,1,2
        A[0][0] = -h[1];
        A[0][1] = h[0] + h[1];
        A[0][2] = -h[0];
        bb[0] = 0.0;

        // 中间标准三对角方程 i = 1 .. n-2:
        // h_{i-1} M_{i-1} + 2*(h_{i-1}+h_i) M_i + h_i M_{i+1} = 6*(first_diff[i] - first_diff[i-1])
        for (int i = 1; i <= n-2; ++i) {
            A[i][i-1] = h[i-1];
            A[i][i]   = 2.0 * (h[i-1] + h[i]);
            A[i][i+1] = h[i];
            // RHS: 6*( (f[i+1]-f[i])/h[i] - (f[i]-f[i-1])/h[i-1] )
            bb[i] = 6.0 * ( first_difference[i] - first_difference[i-1] );
        }

        // 末行 (Not-a-knot at x_{n-2}):
        // h_{n-3} * M_{n-1} - (h_{n-3}+h_{n-2}) * M_{n-2} + h_{n-2} * M_{n-3} = 0
        // 将其放在行 n-1，列 n-3,n-2,n-1
        A[n-1][n-3] = h[n-2];
        A[n-1][n-2] = -(h[n-3] + h[n-2]);
        A[n-1][n-1] = h[n-3];
        bb[n-1] = 0.0;

        // 求解密线性系统
        vector<double> Msol = solveDenseLinearSystem(A, bb);

        // 复制到 M
        for (int i = 0; i < n; ++i) M[i] = Msol[i];
    }
    return SplineWrapper(PPform(t, f, M));
}

#include <iostream>
#include <vector>
#include <cmath>

// 测试函数，用于检查样条插值是否正确
void testCubicSpline() {
    // 测试数据
    vector<double> t = {0, 1, 2, 3, 4}; // 节点
    vector<double> f = {0, 0.5, 2, 4.5, 8}; // 函数值
    vector<double> bdVal = {0, 0}; // 边界条件值，例如 Natural 边界条件

    // 获取样条插值函数
    auto spline = getCubicSpline(t, f, "Natural", bdVal);

    // 测试点
    vector<double> testPoints = {0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4};
    vector<double> expectedValues = {0, 0.25, 0.5, 1.25, 2, 3.75, 4.5, 6.25, 8};

    // 进行测试
    cout << "Testing Cubic Spline..." << endl;
    for (size_t i = 0; i < testPoints.size(); ++i) {
        double x = testPoints[i];
        double actualValue = spline(x);
        double expectedValue = expectedValues[i];
        double error = fabs(actualValue - expectedValue);

        if (error < 1e-5) {
            cout << "Test point " << x << " passed. Expected: " << expectedValue << ", Got: " << actualValue << endl;
        } else {
            cout << "Test point " << x << " failed. Expected: " << expectedValue << ", Got: " << actualValue << endl;
        }
    }
}

int main() {
    testCubicSpline();
    return 0;
}