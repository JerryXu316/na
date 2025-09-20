#include <iostream>
using namespace std;

int main() {
    int a = 5;
    double b = 2.7;

    // 隐式转换：int -> double
    double sum = a + b;
    cout << "a + b = " << sum << endl;

    // 显式转换：double -> int
    int c = (int)b;
    cout << "b 转换为 int 后: " << c << endl;

    // 函数风格
    int d = int(b);
    cout << "函数风格转换 b -> int: " << d << endl;

    return 0;
}
