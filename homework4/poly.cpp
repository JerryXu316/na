#include <vector>
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// 返回的对象类型
class PolyMSP {
private:
    int p;                                  // 多项式阶数
    VectorXd c;                             // 系数 c_0 .. c_p
    vector<double> x_data, y_data;          // 原始数据，用来计算残差

public:
    // 构造函数：不再需要存储 basis 对象，只需要系数和阶数
    PolyMSP(int p,
            const VectorXd &c,
            const vector<double> &x_data,
            const vector<double> &y_data)
        : p(p), c(c), x_data(x_data), y_data(y_data) {}

    // 1) 小括号运算：计算 f(x)
    // 优化：使用递推公式在线计算 P_k(x) 的值，避免 O(p^2) 的重复计算，降为 O(p)
    double operator()(double x) const {
        if (p == 0) return c(0);
        
        double Pkm2 = 1.0; // P_{k-2}, 对应 P_0
        double Pkm1 = x;   // P_{k-1}, 对应 P_1
        
        double val = c(0) * Pkm2;
        if (p >= 1) val += c(1) * Pkm1;

        for (int k = 2; k <= p; ++k) {
            // Legendre 递推: k*P_k = (2k-1)*x*P_{k-1} - (k-1)*P_{k-2}
            double Pk = ((2.0 * k - 1.0) * x * Pkm1 - (k - 1.0) * Pkm2) / (double)k;
            val += c(k) * Pk;
            
            // 更新前两项
            Pkm2 = Pkm1;
            Pkm1 = Pk;
        }
        return val;
    }

    // 2) 返回残差 sum (y_i - f(x_i))^2
    double getRes() const {
        double r = 0.0;
        for (size_t i = 0; i < x_data.size(); ++i) {
            double diff = y_data[i] - (*this)(x_data[i]);
            r += diff * diff;
        }
        return r;
    }
};

auto getPolyMSP(int p, 
                const vector<double> &x, 
                const vector<double> &y) {
    int n = x.size();

    // 创建一个 n 行 p + 1 列的 double 矩阵 A
    MatrixXd A(n, p + 1);

    // 优化：在填充矩阵时直接使用递推公式，复杂度从 O(np^2) 降为 O(np)
    for (int i = 0; i < n; ++i) {
        double xi = x[i];
        
        // P_0(x) = 1
        A(i, 0) = 1.0;
        
        if (p >= 1) {
            // P_1(x) = x
            A(i, 1) = xi;
        }

        // P_k(x) 递推
        for (int j = 2; j <= p; ++j) {
            // 利用矩阵中已经计算好的前两列 A(i, j-1) 和 A(i, j-2)
            A(i, j) = ((2.0 * j - 1.0) * xi * A(i, j - 1) - (j - 1.0) * A(i, j - 2)) / (double)j;
        }
    }

    // y 转换成 Eigen 向量
    Map<const VectorXd> yv(y.data(), n);

    // 对矩阵 A 进行 Householder QR 分解
    HouseholderQR<MatrixXd> qr(A);
    
    // 优化：不要显式计算 Q (n*n 矩阵)，直接计算 Q^T * y
    // Eigen 会智能处理 householderQ() * vector，复杂度 O(np) 且不消耗巨大内存
    VectorXd qTy = qr.householderQ().transpose() * yv;

    // 取前 p+1 个元素作为 b1 (对应数学公式中的 R1 * c = Q1^T * y)
    VectorXd b1 = qTy.head(p + 1);

    // 取 R 的上三角前 (p+1)x(p+1) 并求解
    // qr.matrixQR() 存储了 R (上三角) 和 Householder 向量 (下三角)
    VectorXd c = qr.matrixQR().topRows(p + 1).triangularView<Upper>().solve(b1);

    return PolyMSP(p, c, x, y);
}

// 本地编译测试
void compileTest() {
    vector<double> x = {-0.2, -0.5, 0, 0.5, 0.2};
    auto f = getPolyMSP(1, x, x);
    cout << "f(0.1) = " << f(0.1) << ", Residual = " << f.getRes() << endl;
}

int main() {
    compileTest();
    return 0;
}