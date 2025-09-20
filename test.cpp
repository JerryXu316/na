#include <iostream>
#include <string>

using namespace std;

struct Kunkun{
    string hobby1 = "唱";
    string hobby2 = "跳";
    string hobby3 = "rap";
    string hobby4 = "篮球";

    string song1 = "只因你太美";
    string song2 = "Deadman";

    void 被打了(){
        cout << "坤坤：哎呦，你干嘛~" << endl;
    }

    void 见到熟人(){
        cout << "坤坤：哎呦，真的是你呀~" << endl;
    }

    void 展示爱好(){
        cout << "坤坤的爱好是："
             << hobby1 << "、"
             << hobby2 << "、"
             << hobby3 << "、"
             << hobby4 << endl;
    }

    void 唱歌(){
        int song;
        cout << "请告诉我你想要唱哪首歌？" << endl;
        cout << "1. 只因你太美" << endl;
        cout << "2. Deadman" << endl;
        cout << "3. 退出" << endl;
        cout << "输入你的选择：";
        cin >> song;

        switch (song){
            case 1:     
                cout << "只因你太美!" << endl;
                break;
            case 2:
                cout << "I am a deadman!" << endl;
                break;
            case 3:
                cout << "哎呦你干嘛" << endl;
                break;
            default:
                cout << "无效输入，请重新选择！" << endl;
        }
    }
};

int main(){
    Kunkun k;
    int choice;

    cout << "=== 欢迎进入坤坤交互系统 ===" << endl;

    while (true) {
        cout << "\n请选择操作：" << endl;
        cout << "1. 展示爱好" << endl;
        cout << "2. 打坤坤" << endl;
        cout << "3. 遇到熟人" << endl;
        cout << "4. 唱歌" << endl;
        cout << "0. 退出" << endl;
        cout << "输入你的选择：";
        cin >> choice;

        switch (choice){
            case 1:
                k.展示爱好();
                break;
            case 2:
                k.被打了();
                break;
            case 3:
                k.见到熟人();
                break;
            case 4:
                k.唱歌();
                break;
            case 0:
                cout << "系统已退出，再见！" << endl;
                return 0;
            default:
                cout << "无效输入，请重新选择！" << endl;
        }
    }

    return 0;
}