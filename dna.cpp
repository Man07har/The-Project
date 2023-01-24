#include <iostream>
using namespace std;

int main() {
    char t[] = "AAGGAG";
    char p[] = "AG";
    int n = sizeof(t) / sizeof(t[0]) - 1;
    int m = sizeof(p) / sizeof(p[0]) - 1;
    int window_index[n];
    int count = 0;
    int num_window = 0;

    while (count <= n - m) {
        if (t[count] == p[0]) {
            if (t[count + m - 1] == p[m - 1]) {
                window_index[num_window] = count;
                num_window++;
            }
        }
        count++;
    }

    cout << "The pattern was found at starting indexes: ";
    for (int i = 0; i < num_window; i++) {
        cout << window_index[i] << " ";
    }
    cout << endl;
    cout<<"\n  number of count was: "<<count<<" "<<endl;
    return 0;
}