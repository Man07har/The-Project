#include <iostream>
using namespace std;

int main() {
    string t;
    string p;
   cout<<"Enter the string: "<<endl;
   getline(cin,t);
   cout<<"Enter the pattern: "<<endl;
   getline(cin,p);
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

    count=0;
    int num_match=0;
    int match_index[num_window];
    while(count<num_window){
        int s=window_index[count];
       int c=1;
       while(c<=m-2){
        if(p[c]!=t[s+c]){
            break;
        }
        c=c+1;
       }
    if(c=m-1){
        match_index[num_match] =s;
        num_match=num_match+1;
    }
    count=count+1;
    }
    return 0;
}