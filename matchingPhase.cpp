#include <iostream>
using namespace std;

int main() {
    int count = 0;
    int num_match = 0;
    int num_window = 100; 
    int window_index[num_window];
    int m = 5; 
    char p[m]; 
    char t[1000]; 
    int match_index[num_window];
    
    while (count < num_window) {
        int s = window_index[count];
        int c = 1;
        while (c <= m - 2) {
            if (p[c] != t[s + c]) {
                break;
            }
            c++;
        }
        if (c == m - 1) {
            match_index[num_match] = s;
            num_match++;
        }
        count++;
    }
    return 0;
}
