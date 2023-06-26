#include <fstream>
#include <iostream>
#include <vector>
#include <thread>
#include <chrono>

#include <mpi.h>
#include <omp.h>

using namespace std;

const int MASTER = 0;

vector <int> euler; // Euler tour of a tree
vector <int> level; // Level of a node in the Euler tour
vector <int> first; // First occurence of a node in the Euler tour

vector <vector<int>> rmq; 
vector <int> lg; // Vector of logarithms

vector <vector<int>> queries;
vector <vector<int>> query_lines;
vector <int> outputs;

string in_file = "../io_files/inputs/large.in";
string out_file = "../io_files/outputs/distributed/large.out";

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
void build_rmq_line(int i, int tour_length) {
    int j;
    #pragma omp parallel for private(j)
    for (j = (1 << i); j <= tour_length; ++j) {
        if ((1 << i) <= j) {
            rmq[i][j] = rmq[i - 1][j];
            if (level[rmq[i][j]] > level[rmq[i - 1][j - (1 << (i - 1))]]) {
                rmq[i][j] = rmq[i - 1][j - (1 << (i - 1))];
            }
        }
    }
}

int main(int argc, char **argv) {
    int n; // Tree size
    int q; // Number of queries
    
    omp_set_num_threads(NUM_THREADS);
	int rank, proc_count;

    // initializam procesele
    MPI_Init(&argc, &argv);
    // fiecare proces are un rank
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // rank
    // stabilim numarul de procese
    MPI_Comm_size(MPI_COMM_WORLD, &proc_count);

    if (rank == MASTER) {
        ifstream fin(in_file);
        fin >> n >> q;
        fin.close();
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&q, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int tour_length;
    first.assign(n + 1, 0);
    queries.assign(2, vector<int>(q + 1, 0));

    if (rank == MASTER) {
        ifstream fin(in_file);
        fin >> n >> q;
        vector<vector<int>> children(n + 1); // List of children for each node
		vector<vector<int>> g(n + 1);

		for (int i = 2; i <= n; ++i) {
			int x; // Father of node i
			fin >> x;
			g[i].push_back(x);
			g[x].push_back(i);
		}

        for (int i = 1; i <= q; ++i) {
			fin >> queries[0][i] >> queries[1][i];
		}

		get_tree(1, n, g, children);

		euler.push_back(0); // index from 1
		level.push_back(0); // index from 1

		generate_euler_tour(children, euler, level, first);
		tour_length = euler.size() - 1; // Length of the Euler tour

		lg.assign(tour_length + 1, 0);
		build_lg(tour_length, lg);
    }
	
	// Send euler tour data
	MPI_Bcast(&tour_length, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

	if (rank != MASTER) {
		euler.resize(tour_length + 1);
		level.resize(tour_length + 1);
		lg.resize(tour_length + 1);
	}
    MPI_Bcast(&first[0], n + 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(&euler[0], tour_length + 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(&level[0], tour_length + 1, MPI_INT, MASTER, MPI_COMM_WORLD);

	// Send logarithms and prepare RMQ
	MPI_Bcast(&lg[0], tour_length + 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    rmq.assign(lg[tour_length] + 1, vector <int> (tour_length + 1, 0));

    // Send queries to workers
	MPI_Bcast(&queries[0][0], q + 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(&queries[1][0], q + 1, MPI_INT, MASTER, MPI_COMM_WORLD);

    outputs.assign(q + 1, 0);
    // OTHERS will now separate the queries based on the line they need 
    if (rank != MASTER) {
        query_lines.assign(lg[tour_length] + 1, vector<int>());

        int start_idx = 1 + (rank - 1) * (double) q / (proc_count - 1);
        int end_idx = 1 + (rank) * (double) q / (proc_count - 1);

        if (end_idx > q + 1) {
            end_idx = q + 1;
        }


        for (int i = start_idx; i < end_idx; ++i) {
            int x = first[queries[0][i]];
            int y = first[queries[1][i]];

            if (x > y) {
                swap(x, y);
            }

            int mid = lg[y - x + 1];
            query_lines[mid].push_back(i);
        }

    }

    // From now on, MASTER will build the RMQ, while OTHERS will answer the queries
    if (rank == MASTER) {
        for (int i = 1; i <= tour_length; ++i) {
			rmq[0][i] = i;
		}
    }
    
    // MASTER will build the i-th line, OTHERS will answer queries from the i-th line
	for (int i = 1; i <= lg[tour_length]; ++i) {

        if (rank == MASTER) {
		    build_rmq_line(i, tour_length);
        }

        MPI_Bcast(&rmq[i][0], tour_length+1, MPI_INT, MASTER, MPI_COMM_WORLD);

        if (rank != MASTER) {
            int j, len = query_lines[i].size();
            #pragma omp parallel for private(j)
            for (j = 0; j < len; ++j) {
                int idx = query_lines[i][j];
                int x = first[queries[0][idx]];
                int y = first[queries[1][idx]];

                if (x > y) {
                    swap(x, y);
                }

                int mid = lg[y - x + 1];
                int ans1 = rmq[mid][y];
                int ans2 = rmq[mid][x + (1 << mid) - 1];

                if (level[ans1] > level[ans2]) {
                    ans1 = ans2;
                }

                outputs[idx] = euler[ans1];
            }
        }
	}

    // Send segment back to master
    if (rank != MASTER) {
        int start_idx = 1 + (rank - 1) * (double) q / (proc_count - 1);
        int end_idx = 1 + (rank) * (double) q / (proc_count - 1);
        MPI_Send(&outputs[start_idx], end_idx-start_idx, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
    } else {
        for (int j = 1; j < proc_count; ++j) {
            MPI_Status status;
            vector<int> recv_buff(q/(proc_count - 1) + 10, 0);

            MPI_Recv(&recv_buff[0], q / (proc_count - 1) + 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

            int from = status.MPI_SOURCE;
            int start_segm = 1 + (from - 1) * (double) q / (proc_count - 1);
            int end_segm = 1 + (from) * (double) q / (proc_count - 1);
            if (end_segm > q + 1) {
                end_segm = q + 1;
            }

            // Copy received segment in master memory
            int i;
            #pragma omp parallel for private(i)
            for (i = start_segm; i < end_segm; ++i) {
                outputs[i] = recv_buff[i-start_segm];
            }
	    }
    }

	// Print queries
	if (rank == MASTER) {
        ofstream fout(out_file);
		for (int i = 1; i <= q; ++i) {
			fout << outputs[i] << "\n";
		}
	}

	MPI_Finalize();

    return 0;
}
