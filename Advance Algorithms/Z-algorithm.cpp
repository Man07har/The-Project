#include <iostream>
#include <vector>

using namespace std;

vector<int> z_algorithm(const string& str) {
    int n = str.length();
    vector<int> Z(n, 0);

    int l = 0, r = 0;
    for (int i = 1; i < n; ++i) {
        if (i <= r) {
            Z[i] = min(r - i + 1, Z[i - l]);
        }

        while (i + Z[i] < n && str[Z[i]] == str[i + Z[i]]) {
            Z[i]++;
        }

        if (i + Z[i] - 1 > r) {
            l = i;
            r = i + Z[i] - 1;
        }
    }

    return Z;
}

int main() {
    string text = "abababab";
    string pattern = "aba";

    // Concatenate pattern and text with a special character
    string concatStr = pattern + "$" + text;

    vector<int> Z = z_algorithm(concatStr);

    // Search for occurrences of the pattern in the Z array
    for (int i = 0; i < Z.size(); ++i) {
        if (Z[i] == pattern.length()) {
            cout << "Pattern found at index " << i - pattern.length() - 1 << endl;
        }
    }

    return 0;
}
