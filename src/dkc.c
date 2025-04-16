/**
 * File: pivoter.c
 * Author: Cui Donghang
 * Last Update: 2022/04/24
 * 
 * This file implements all function declared in pivoter.h.
 */

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <mpi.h>
#include <omp.h>
#include <time.h>

#include "dkc.h"

double comb[2000][2000];

/**
 * Read graph data from a specified file.
 * The first line of the file should contain two number: max{v} in V and |E|, 
 * and each of the next |E| lines should contain two vertices connected by an edge.
 * Note that "max{v} in V" is not always identical to |V|!
 * Return the pointer of the initialized struct Graph.
 * 
 */

Graph *read_graph(const char *filename) {

	freopen(filename, "r", stdin);
	int n, m, u, v;
	scanf("%d%d", &n, &m);

	Graph *g = init_graph(n, (unsigned int)(m << 1));

	for (int i = 1; i <= m; ++i) {
		scanf("%d%d", &u, &v);
		++u; ++v;
		add_edge(g, u, v);
		add_edge(g, v, u);
	}

	return g;

}

/**
 * Initialize the combination number from C(1, 1) to C(n, n).
 */

void comb_init(int n) {
	for (int i = 0; i <= n; ++i)
		comb[i][0] = 1.;

	for (int i = 1; i <= n; ++i)
		for (int j = 1; j <= i; ++j) {
			comb[i][j] = comb[i - 1][j - 1] + comb[i - 1][j];
		}
}

Graph *degree_orientation(Graph *g) {
	Graph *d = init_graph(g->v_num, g->e_num >> 1);

	for (int u = 1; u <= g->v_num; ++u)
		for (unsigned int i = g->first[u]; ~i; i = g->e[i].nxt) {
			int v = g->e[i].dst;
			if (g->in_deg[u] < g->in_deg[v] || (g->in_deg[u] == g->in_deg[v] && u < v))
				add_edge(d, u, v);
		}
	return d;
}

/**
 * An implementation of core/degeneracy decomposition in O(m) time.
 * Orient the undirected graph for lowest out degree.
 * Return the pointer of the oriented (directed) graph.
 */
Graph *core_decomposition(Graph *g) {
	Graph *d = init_graph(g->v_num, g->e_num >> 1);
	int *bin = (int *)calloc(g->max_in_deg + 2, sizeof(int));
	int *rank = (int *)calloc(g->v_num + 1, sizeof(int));
	int *rank_pos = (int *)calloc(g->v_num + 1, sizeof(int));
	int *deg = (int *)calloc(g->v_num + 1, sizeof(int));
	memcpy(deg, g->in_deg, sizeof(int) * (g->v_num + 1));


	for (int i = 1; i <= g->v_num; ++i)
		bin[deg[i] + 1]++;

	for (int i = 1; i <= g->max_in_deg; ++i) {
		bin[i] += bin[i - 1];
	}
	for (int i = 1; i <= g->v_num; ++i) {
		rank[i] = ++bin[deg[i]];
		rank_pos[rank[i]] = i;
	}

	for (int i = 1; i <= g->v_num; ++i) {
		int u = rank_pos[i];
		for (unsigned int j = g->first[u]; ~j; j = g->e[j].nxt) {
			int v = g->e[j].dst;
			if (deg[v] <= deg[u]) continue; 
			int w = rank_pos[++bin[--deg[v]]];
			rank[w] = rank[v];
			rank[v] = bin[deg[v]];
			rank_pos[rank[v]] = v;
			rank_pos[rank[w]] = w;
		}
	}

	for (int u = 1; u <= g->v_num; ++u) {
		for (unsigned int i = g->first[u]; ~i; i = g->e[i].nxt) {
			int v = g->e[i].dst;
			if (rank[u] < rank[v])
				add_edge(d, u, v);
		}
	}
		

	free(bin);
	free(rank);
	free(rank_pos);
	free(deg);

	return d;
}

/**
 * Init the set S and S_p with 1 ... vertex_num.
 * Note that S[0] is the size of S.
 */
void S_init(int *S, int *S_p, int vertex_num) {
	S[0] = vertex_num;
	for (int i = 1; i <= vertex_num; ++i) {
		S[i] = i;
		S_p[i] = i;
	}
}

/**
 * Insert a vertex into the set S.
 * Note that S[0] is the size of S.
 */
void S_insert(int *S, int *S_p, int vertex) {
	int v_p = S_p[vertex];
	if (v_p <= S[0]) return;
	S[v_p] = S[++S[0]];
	S[S[0]] = vertex;
	S_p[vertex] = S[0];
	S_p[S[v_p]] = v_p;
}

/**
 * Delete a vertex from the set S.
 * Note that S[0] is the size of S.
 */
void S_delete(int *S, int *S_p, int vertex) {
	int v_p = S_p[vertex];
	if (v_p > S[0]) return;
	S[v_p] = S[S[0]];
	S[S[0]] = vertex;
	S_p[vertex] = S[0]--;
	S_p[S[v_p]] = v_p;
}

/**
 * This is a depth-first implementation of Pivoter.
 * Params:
 * * S: the vertex set of a pivoter call
 * * S_p: the index of vertices in S
 * * S_nbr: the neighbor information of vertices in S
 * * S_cand: the candidates from S for the next call
 * * C: the array that stores the clique count
 * * k: the max k value
 * * c_cnt: the number of vertices in the current clique
 * * p_cnt: the number of pivot vertices in the current clique
 */
void pivoter_dfs(int *S, int *S_p, int **S_nbr, int **S_cand, double *C, int c_cnt, int p_cnt, int k) {

	// S is empty or the number of hold vertices reaches k
	if (!S[0] || c_cnt - p_cnt == k) {
		for (int i = 0; i <= min(p_cnt, k - c_cnt + p_cnt); ++i) {
# 			pragma omp atomic
			C[i + (c_cnt - p_cnt)] += comb[p_cnt][i];
		}
		return;
	}


	// Choose a pivot
	int pivot = 0, max_nbr_cnt = -1;
	for (int i = 1; i <= S[0]; ++i) {
		int u = S[i], nbr_cnt = 0;
		for (int j = 1; j <= S_nbr[u][0]; ++j) {
			int v = S_nbr[u][j];
			if (S_p[v] <= S[0])
				++nbr_cnt;
			else
				break;
		}
		if (nbr_cnt > max_nbr_cnt) {
			max_nbr_cnt = nbr_cnt;
			pivot = u;
		}
	}

	int *cand = S_cand[c_cnt];

	memcpy(cand, S, sizeof(int) * (S[0] + 1));

	for (int i = 1; i <= S_nbr[pivot][0]; ++i) {
		int u = S_nbr[pivot][i];
		if (S_p[u] <= S[0])
			cand[S_p[u]] = -1;
		else
			break;
	}

	cand[0] = 0;
	for (int i = 1; i <= S[0]; ++i) {
		if (cand[i] == -1) continue;
		cand[++cand[0]] = cand[i];
	}



	// search N(S, candidate)

	int S0_bak = S[0];

	for (int i = 1; i <= cand[0]; ++i) {
		int u = cand[i], S_sz = S[0];
		S[0] = 0;

		for (int j = 1; j <= S_nbr[u][0]; ++j) {
			int v = S_nbr[u][j];
			if (S_p[v] <= S_sz)
				S_insert(S, S_p, v);
			else if (S_p[v] > S0_bak)
				break;

		}

		for (int j = 1; j <= S[0]; ++j) {
			int v = S[j], cnt = 0;
			for (int k = 1; k <= S_nbr[v][0]; ++k) {
				int v_nbr = S_nbr[v][k];
				if (S_p[v_nbr] <= S[0]) {
					S_nbr[v][k] = S_nbr[v][++cnt];
					S_nbr[v][cnt] = v_nbr;
				}
				if (S_p[v_nbr] > S0_bak)
					break;
			}
		}

		pivoter_dfs(S, S_p, S_nbr, S_cand, C, c_cnt + 1, p_cnt + (int)(u == pivot), k);

		S[0] = S_sz;

		S_delete(S, S_p, u);

	}

	S[0] = S0_bak;

}

/**
 * Push a subgraph into buffer
 * Params:
 * * S_nbr_buf: the pointer of buffer
 * * S_nbr: the neighbor list of a subgraph
 * * num_vertice: the number of vertice in the subgraph
 * * buffer_size: the maximum buffer size
 */
void buffer_push(int **S_nbr_buf, int **S_nbr, int num_vertice, int buffer_size) {

	int flag = 0;
	while (!flag) {
#		pragma omp critical(S_nbr_buf_lock)
		{
			// Add current task to buffer

			if (S_nbr_buf[0][0] >= buffer_size) flag = 0;
			else {
				int *buf = S_nbr_buf[++S_nbr_buf[0][0]];
				int buf_sz = 0;

				buf[0] = num_vertice;
				for (int v = 1; v <= num_vertice; ++v) {
					memcpy(buf + buf_sz + 1, S_nbr[v], sizeof(int) * (S_nbr[v][0] + 1));
					buf_sz += S_nbr[v][0] + 1;
				}
				flag = 1;
			}
		}
	}
} 

/**
 * Pop a subgraph out of buffer
 * Params:
 * * S_nbr_buf: the pointer of buffer
 * * msg_buf: the formatted subgraph
 * * num_running_producers: pop condition
 * * num_tasks: pop condition 
 */
int buffer_pop(int **S_nbr_buf, int *msg_buf, int *num_running_producers, int num_tasks) {
	int msg_buf_sz = 0;

	for (msg_buf[0] = 0; msg_buf[0] < num_tasks; ) {
		int flag = 0;

		while (!flag) {
#			pragma omp critical(S_nbr_buf_lock)
			{
				if (S_nbr_buf[0][0] == 0) {
					if (*num_running_producers == 0) flag = 2;
					else flag = 0;
				} else {
					for (; S_nbr_buf[0][0] > 0 && msg_buf[0] <= num_tasks; --S_nbr_buf[0][0], ++msg_buf[0]) {
						int *buf = S_nbr_buf[S_nbr_buf[0][0]];

						msg_buf[++msg_buf_sz] = buf[0];

						for (int u = 1, i = 1; u <= buf[0]; ++u, i += buf[i] + 1) {
							memcpy(msg_buf + msg_buf_sz + 1, buf + i, sizeof(int) * (buf[i] + 1));
							msg_buf_sz += buf[i] + 1;
						}

					}
					flag = 1;
				}
			}
		}

		if (flag == 2) break;
	}
	return msg_buf_sz;
}

/**
 * Divide k-clique task into multiple subtasks with orientation method
 * Params:
 * * S: the vertex set of a recursion call
 * * S_p: the index of vertices in S
 * * S_nbr: the neighbor information of vertices in S
 * * S_cand: the candidates from S for the next call
 * * new_S_nbr: the new S_nbr array of the discretized graph
 * * new_v: the new vertex id of the discretized graph
 * * S_nbr_buf: the pointer of buffer
 * * level: division level
 * * buffer_size: the size of buffer
 */
void orientation_dfs(int *S, int *S_p, int **S_nbr, int **S_cand, int **new_S_nbr, int *new_v, int **S_nbr_buf, int level, int buffer_size) {
	if (level == 0) {
		
		for (int i = 1; i <= S[0]; ++i) {
			new_v[S[i]] = ++new_v[0];
			new_S_nbr[new_v[0]][0] = 0;
		}

		for (int i = 1; i <= S[0]; ++i) {
			int v = S[i];
			int nv = new_v[v];
			for (int j = 1; j <= S_nbr[v][0]; ++j) {
				int v_nbr = S_nbr[v][j];
				int nv_nbr = new_v[v_nbr];
				if (nv_nbr) {
					new_S_nbr[nv][++new_S_nbr[nv][0]] = nv_nbr;
					new_S_nbr[nv_nbr][++new_S_nbr[nv_nbr][0]] = nv;
				}
			}
		}


		buffer_push(S_nbr_buf, new_S_nbr, new_v[0], buffer_size);

		new_v[0] = 0;
		for (int i = 1; i <= S[0]; ++i)
			new_v[S[i]] = 0;
		return;
	}

	int *cand = S_cand[level];

	memcpy(cand, S, sizeof(int) * (S[0] + 1));

	int S0_bak = S[0];

	for (int i = 1; i <= cand[0]; ++i) {

		int u = cand[i], S_sz = S[0];
		S[0] = 0;

		for (int j = 1; j <= S_nbr[u][0]; ++j) {
			int v = S_nbr[u][j];
			if (S_p[v] <= S_sz)
				S_insert(S, S_p, v);
			else if (S_p[v] > S0_bak)
				break;

		}

		for (int j = 1; j <= S[0]; ++j) {
			int v = S[j], cnt = 0;
			for (int k = 1; k <= S_nbr[v][0]; ++k) {
				int v_nbr = S_nbr[v][k];
				if (S_p[v_nbr] <= S[0]) {
					S_nbr[v][k] = S_nbr[v][++cnt];
					S_nbr[v][cnt] = v_nbr;
				}
				if (S_p[v_nbr] > S0_bak)
					break;
			}
		}

		orientation_dfs(S, S_p, S_nbr, S_cand, new_S_nbr, new_v, S_nbr_buf, level - 1, buffer_size);

		S[0] = S_sz;

	}

	S[0] = S0_bak;


}

/**
 * Start the Distributed K-clique Counting algorithm.
 * Params:
 * * input_file: the filename of the input graph data
 * * master_num_threads: the number of threads running on master node
 * * worker_num_threads: the number of threads running on worker node
 * * buffer_size: the size of sending buffer
 * * num_tasks: the number of tasks sent from master node to worker node each time
 * * k: the maximum value of k
 */

void dkc(const char *input_file, int master_num_threads, int worker_num_threads, int buffer_size, int num_tasks, int task_level, int k) {

	int my_rank, comm_sz;

	int time_read, time_orient, time_calc;

	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	Graph *d = NULL;
	int deg;

	// Init input data
	if (my_rank == MASTER_NODE) {
		time_read = time(NULL);
		Graph *g = read_graph(input_file);
		time_read = time(NULL) - time_read;
		fprintf(stderr, "Running time of reading graph: %ds\n", time_read);
		time_orient = time(NULL);
		d = core_decomposition(g);
		//d = degree_orientation(g);
		time_orient = time(NULL) - time_orient;
		fprintf(stderr, "Running time of orientation: %ds\n", time_orient);
		fprintf(stderr, "max degree = %d\ndegeneracy = %d\n", g->max_out_deg, d->max_out_deg);
		free_graph(g);
		deg = d->max_out_deg;
	}

	// Broadcast the degeneracy to every worker node
	MPI_Bcast(&deg, 1, MPI_INT, 0, MPI_COMM_WORLD);

	double *C = (double *)calloc(deg + 2, sizeof(double));
	double *C_sum = (double *)calloc(deg + 2, sizeof(double));

	// Init combination
	comb_init(deg + 1);

// ************ Master node process ************
	if (my_rank == MASTER_NODE) {

		time_calc = time(NULL);

		int curr_v = 1;
		int num_threads = max(master_num_threads + 1, 2);
		int num_running_producers = num_threads - 1;

		int **S_nbr_buf = (int **)calloc(buffer_size + 1, sizeof(int *));
		for (int i = 0; i <= buffer_size; ++i)
			S_nbr_buf[i] = (int *)calloc((deg + 1) * (deg + 1), sizeof(int));

#		pragma omp parallel num_threads(num_threads)
		{
			int t = omp_get_thread_num();

			if (t > 0) { // Producer

				int n = d->v_num;

				int *new_v = (int *)calloc(n + 1, sizeof(int));
				int **S_nbr = (int **)calloc(deg + 1, sizeof(int *));
				for (int i = 0; i <= deg; ++i)
					S_nbr[i] = (int *)calloc(deg + 1, sizeof(int));

				int **new_S_nbr = (int **)calloc(deg + 1, sizeof(int *));
				for (int i = 0; i <= deg; ++i)
					new_S_nbr[i] = (int *)calloc(deg + 1, sizeof(int));
				

				int *S = (int *)calloc(deg + 1, sizeof(int));
				int *S_p = (int *)calloc(deg + 1, sizeof(int));
				int **S_cand = (int **)calloc(deg + 2, sizeof(int *));
				for (int i = 0; i <= deg + 1; ++i)
					S_cand[i] = (int *)calloc(deg + 1, sizeof(int));

				

				while (1) {

					int u, flag = 0;

#					pragma omp critical(curr_v_lock)
					{
						if (curr_v > n) flag = 1;
						else {
							u = curr_v++;
							flag = 0;
						}
					}

					if (flag) break;

					for (int i = d->first[u]; ~i; i = d->e[i].nxt) {
						int v = d->e[i].dst;
						new_v[v] = ++new_v[0];
						S_nbr[new_v[v]][0] = 0;
					}

					for (int i = d->first[u]; ~i; i = d->e[i].nxt) {
						int v = d->e[i].dst;
						int nv = new_v[v];
						for (int j = d->first[v]; ~j; j = d->e[j].nxt) {
							int v_nbr = d->e[j].dst;
							int nv_nbr = new_v[v_nbr];
							if (nv_nbr) {
								S_nbr[nv][++S_nbr[nv][0]] = nv_nbr;
								//S_nbr[nv_nbr][++S_nbr[nv_nbr][0]] = nv;
							}
						}
					}

					S_init(S, S_p, new_v[0]);

					// Clear new id of vertex
					new_v[0] = 0;
					for (int i = d->first[u]; ~i; i = d->e[i].nxt) {
						int v = d->e[i].dst;
						new_v[v] = 0;
					}

					if (comm_sz == 1) {
						// Count cliques
						pivoter_dfs(S, S_p, S_nbr, S_cand, C, 1, 0, k);
					}
					else {
						// Generate tasks
						orientation_dfs(S, S_p, S_nbr, S_cand, new_S_nbr, new_v, S_nbr_buf, task_level - 1, buffer_size);
					}

				}

#				pragma omp atomic
				--num_running_producers;

				free(S);
				free(S_p);
				free(new_v);
				for (int i = 0; i <= deg; ++i) {
					free(S_nbr[i]);
					free(new_S_nbr[i]);
				}
				free(S_nbr);
				free(new_S_nbr);
				for (int i = 0; i <= deg + 1; ++i) free(S_cand[i]);
				free(S_cand);


			}
			else { // Thread for sending message

				if (comm_sz > 1) {
					// Send tasks to worker nodes

					int *msg_buf = (int *)calloc((deg + 1) * (deg + 1) * (num_tasks + 1), sizeof(int));
					
					while (1) {

						int msg_buf_sz = buffer_pop(S_nbr_buf, msg_buf, &num_running_producers, num_tasks);

						int rk;
						if (msg_buf[0]) {
							MPI_Recv(&rk, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							MPI_Send(msg_buf, msg_buf_sz + 1, MPI_INT, rk, 1, MPI_COMM_WORLD);
						}
						else {
							for (int i = 0; i < worker_num_threads * (comm_sz - 1); ++i) {
								int tmp = -1;
								MPI_Recv(&rk, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								MPI_Send(&tmp, 1, MPI_INT, rk, 1, MPI_COMM_WORLD);
							}
							break;
						}
					
					}

					free(msg_buf);

				}

			}
		}
		
		for (int i = 0; i <= buffer_size; ++i) free(S_nbr_buf[i]);
		free(S_nbr_buf);
		free_graph(d);
	} 

// ************ Worker node process ************
	else 
	// Start sub-threads
#	pragma omp parallel num_threads(worker_num_threads)
	{

		int *S = (int *)calloc(deg + 1, sizeof(int));

		int *S_p = (int *)calloc(deg + 1, sizeof(int));
					
		int **S_nbr = (int **)calloc(deg + 1, sizeof(int *));
		for (int i = 0; i <= deg; ++i)
			S_nbr[i] = (int *)calloc(deg + 1, sizeof(int));

		int **S_cand = (int **)calloc(deg + 2, sizeof(int *));
		for (int i = 0; i <= deg + 1; ++i)
			S_cand[i] = (int *)calloc(deg + 1, sizeof(int));

		int *msg_buf = (int *)calloc((deg + 1) * (deg + 1) * (num_tasks + 1), sizeof(int));
				
		while(1) {
#			pragma omp critical(recv_msg)
			{
				// Receive data from master node
				MPI_Send(&my_rank, 1, MPI_INT, MASTER_NODE, 0, MPI_COMM_WORLD);
				MPI_Recv(msg_buf, (deg + 1) * (deg + 1) * (num_tasks + 1), MPI_INT, MASTER_NODE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}

			// Termination signal received
			if (msg_buf[0] == -1) break;

			for (int i = 1, msg_buf_i = 1; i <= msg_buf[0]; ++i) {
				int n = msg_buf[msg_buf_i++];
				for (int u = 1; u <= n; ++u) {
					memcpy(S_nbr[u], msg_buf + msg_buf_i, sizeof(int) * (msg_buf[msg_buf_i] + 1));
					msg_buf_i += msg_buf[msg_buf_i] + 1;
				}

				S_init(S, S_p, n);
				pivoter_dfs(S, S_p, S_nbr, S_cand, C, task_level, 0, k);
			}

			// for (int i = 1, S_nbr_buf_i = 1; i <= S_buf[t][0]; ++i) {

			// 	int n = S_buf[t][i];

			// 	for (int u = 1; u <= n; ++u) {
			// 		memcpy(S_nbr[t][u], S_nbr_buf[t] + S_nbr_buf_i, sizeof(int) * (S_nbr_buf[t][S_nbr_buf_i] + 1));
			// 		S_nbr_buf_i += S_nbr_buf[t][S_nbr_buf_i] + 1;
			// 	}

			// 	S_init(S[t], S_p[t], n);
			// 	pivoter_dfs(S[t], S_p[t], S_nbr[t], S_cand[t], C, 1, 0, k);

			// }
	

		}

		free(S);
		free(S_p);
		for (int i = 0; i <= deg; ++i) free(S_nbr[i]);
		free(S_nbr);
		for (int i = 0; i <= deg + 1; ++i) free(S_cand[i]);
		free(S_cand);
		free(msg_buf);

	}

	// Reduce all count
	MPI_Reduce(C, C_sum, deg + 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	// Print all k-clique counts
	if (my_rank == MASTER_NODE) {
		time_calc = time(NULL) - time_calc;
		fprintf(stderr, "Running time of counting cliques: %ds\n", time_calc);
		fprintf(stderr, "\nResult: \n");
		for (int i = 1; i <= deg + 1; ++i)
			if (C_sum[i] > 1e-8) {
				fprintf(stderr, "C(%d) = %.0lf\n", i, C_sum[i]);
			}
	}

	// Free all allocated memory
	free(C);
	free(C_sum);
}
