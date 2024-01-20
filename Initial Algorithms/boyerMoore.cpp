#include <iostream>
#include <vector>
#include <algorithm>

const int CHAR_SET_SIZE = 256;

// Function to pre-process the pattern and create a bad character heuristic
void badCharHeuristic(const std::string& pattern, std::vector<int>& badChar) {
    int m = pattern.length();
    for (int i = 0; i < CHAR_SET_SIZE; i++)
        badChar[i] = -1;

    for (int i = 0; i < m; i++)
        badChar[static_cast<int>(pattern[i])] = i;
}

// Boyer-Moore algorithm for string searching
void boyerMoore(const std::string& text, const std::string& pattern) {
    int m = pattern.length();
    int n = text.length();

    std::vector<int> badChar(CHAR_SET_SIZE, 0);
    badCharHeuristic(pattern, badChar);

    int s = 0; // s is the shift of the pattern with respect to the text

    while (s <= (n - m)) {
        int j = m - 1;

        // Keep reducing the index of the pattern until there is a mismatch
        while (j >= 0 && pattern[j] == text[s + j])
            j--;

        if (j < 0) {
            // Pattern found, print the index
            std::cout << "Pattern found at index " << s << std::endl;

            // Shift the pattern to find the next occurrence
            s += (s + m < n) ? m - badChar[text[s + m]] : 1;
        } else {
            // Shift the pattern based on the bad character heuristic
            s += std::max(1, j - badChar[text[s + j]]);
        }
    }
}

int main() {
    std::string text = "ABAAABCD";
    std::string pattern = "ABC";
    
    boyerMoore(text, pattern);

    return 0;
}
