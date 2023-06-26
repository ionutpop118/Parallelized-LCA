#include <bits/stdc++.h>
#include <fstream>
#include <cstring>
#include <iostream>

using namespace std;

// Makes only the input file
int main(int argc, char **argv)
{
    if (argc != 3) {
        cout << "Invalid usage: <N> <Q>\n";
        return 0;
    }

    int N, Q;

    N = atoi(argv[1]);
    Q = atoi(argv[2]);
    ofstream fout("../io_files/lca.in");

    mt19937 mt_rand(chrono::system_clock::now().time_since_epoch().count());
    uniform_int_distribution<> nodeRange(N / 5, N);
    
    fout << N << " " << Q << "\n";
    vector <int> v(N + 1);

    v[1] = 1;
    for (int i = 2; i <= N; ++i) {
        uniform_int_distribution<> fatherRange(i / 2, i - 1);
        v[i] = fatherRange(mt_rand);
    }

    for (int i = 2; i < N; ++i)
        fout << v[i] << " ";
    fout << v[N] << "\n";

    for (int i = 1; i <= Q; ++i) {
        int x = nodeRange(mt_rand);
        int y = nodeRange(mt_rand);
        fout << x << " " << y << "\n";
    }

    return 0;
}