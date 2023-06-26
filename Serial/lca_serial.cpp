#include <fstream>
#include <iostream>
#include <vector>

#include <chrono>

using namespace std;
using namespace std::chrono;

string in_file = "../io_files/inputs/large.in";
string out_file = "../io_files/outputs/serial/large.out";

/// Depth-first search of the graph
void build_tree(int node, vector <bool> &viz, vector < vector <int> > &g, vector < vector <int> > &children) {
	viz[node] = 1;
	for (unsigned int i = 0; i < g[node].size(); ++i) {
		if (viz[g[node][i]] == 0) {
			children[node].push_back(g[node][i]);
			build_tree(g[node][i], viz, g, children);
		}
	}
}

/// Transforms the graph into a tree
void get_tree(int root, int &n, vector < vector <int> > &g, vector < vector <int> > &children) {
	vector <bool> viz(n + 1, 0);
	build_tree(root, viz, g, children);
}

void build_lg(int &tour_length, vector <int> &lg) {
	 for (int i = 2; i <= tour_length; ++i) {
        lg[i] = 1 + lg[(i >> 1)]; // Vector of base 2 logarithms
    }
}

/// Depth-first search of the tree
void dfs(int node, int current_level, vector< vector<int> > &children,
         vector <int> &euler, vector <int> &level, vector<int> &first) {
    euler.push_back(node);
    level.push_back(current_level);
    first[node] = euler.size() - 1; // current tour length

    for (unsigned int i = 0; i < children[node].size(); ++i)
    {
        dfs(children[node][i], current_level + 1, children, euler, level, first);
        euler.push_back(node);
        level.push_back(current_level);
    }
}
/// Generates the Euler tour
void generate_euler_tour(vector <vector<int> > &children,
                         vector <int> &euler, vector <int> &level, vector<int> &first) {
    dfs(1, 0, children, euler, level, first); // node 1 is root in this scenario
}

/// Builds the rmq
void build_rmq(vector < vector<int> > &rmq, vector <int> &lg,
               int &tour_length, vector <int> &level) {
    for (int i = 1; i <= tour_length; ++i) {
        rmq[0][i] = i;
    }

    for (int i = 1; i <= lg[tour_length]; ++i) {
        for (int j = 1; j <= tour_length; ++j) {
            if ((1 << i) <= j) {
                rmq[i][j] = rmq[i - 1][j];
                if (level[rmq[i][j]] > level[rmq[i - 1][j - (1 << (i - 1))]]) {
                    rmq[i][j] = rmq[i - 1][j - (1 << (i - 1))];
                }
            }
        }
    }
}
/// Calculates and returns the LCA of 2 nodes x and y
int get_LCA(int x, int y, vector <vector<int> > &rmq, vector <int> &lg,
            vector <int> &euler, vector<int> &first, vector <int> &level) {
    x = first[x];
    y = first[y];
    if (x > y) {
        swap(x, y);
    }
    int mid = lg[y - x + 1];
    int ans = rmq[mid][y];
    if (level[ans] > level[rmq[mid][x + (1 << mid) - 1]]) {
        ans = rmq[mid][x + (1 << mid) - 1];
    }
    return euler[ans];
}

int main() {
    ifstream fin(in_file);
    ofstream fout(out_file);
    int n; // Tree size
    int q; // Number of queries
    // Input
    fin >> n >> q;

    vector< vector<int> > children(n + 1); // List of children for each node
    vector< vector<int> > g(n + 1);
    for (int i = 2; i <= n; ++i) {
        int x; // Father of node i
        fin >> x;
        //children[x].push_back(i);
        g[i].push_back(x);
        g[x].push_back(i);
    }

    get_tree(1, n, g, children);

    vector <int> euler; // Euler tour of a tree
    euler.push_back(0); // index from 1

    vector <int> level; // Level of a node in the Euler tour
    level.push_back(0); // index from 1

    vector <int> first(n + 1); // Level of a node in the Euler tour

    generate_euler_tour(children, euler, level, first);
    int tour_length = euler.size() - 1; // Length of the Euler tour

    vector <int> lg (tour_length + 1, 0);
    build_lg(tour_length, lg);

    vector <vector <int> > rmq(lg[tour_length] + 1, vector <int> (tour_length + 1, 0)); // RMQ
    build_rmq(rmq, lg, tour_length, level);

    vector <pair<int, int>> queries(q + 1, std::make_pair(0, 0));
    vector <int> outputs(q + 1, 0);

    for (int i = 1; i <= q; ++i) {
        fin >> queries[i].first >> queries[i].second;
    }

    for (int i = 1; i <= q; ++i) {
        outputs[i] = get_LCA(queries[i].first, queries[i].second, rmq, lg, euler, first, level);
    }

    for (int i = 1; i <= q; ++i) {
        fout << outputs[i] << "\n";
    }

    
    return 0;
}