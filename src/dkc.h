/**
 * File: pivoter.h
 * Author: Cui Donghang
 * Last Update: 2022/04/24
 * 
 * This file declares all methods for the Distributed Pivoter.
 */

#ifndef PIVOTER_H
#define PIVOTER_H

#include "graph.h"

#define MASTER_NODE 0

int read_num();

Graph *read_graph(const char *filename);

void comb_init(int n);

Graph *degree_orientation(Graph *g);
Graph *core_decomposition(Graph *g);

void S_init(int *S, int *S_p, int vertex_num);
void S_insert(int *S, int *S_p, int vertex);
void S_delete(int *S, int *S_p, int vertex);

void buffer_push(int **S_nbr_buf, int **S_nbr, int num_vertice,  int buffer_size);
int buffer_pop(int **S_nbr_buf, int *msg_buf, int *num_running_producers, int num_tasks);

void orientation_dfs(int *S, int *S_p, int **S_nbr, int **S_cand, int **new_S_nbr, int *new_v, int **S_nbr_buf, int level, int buffer_size);
void pivoter_dfs(int *S, int *S_p, int **S_nbr, int **S_cand, double *C, int c_cnt, int p_cnt, int k);
void dkc(const char *input_file, int master_num_threads, int worker_num_threads, int buffer_size, int num_tasks, int task_level, int k);

#endif // PIVOTER_H
