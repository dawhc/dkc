/**
 * File: graph.h
 * Author: Cui Donghang
 * Last Update: 2022/04/24
 * 
 * This file declares struct Graph and its methods.
 * Graph uses a kind of data structure called Linked Star to store directed/undirected graph.
 */

#ifndef GRAPH_H
#define GRAPH_H

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

typedef struct {
	int src, dst;
        unsigned int nxt;
} Edge; 

typedef struct {
	unsigned int e_num, e_cnt;
	int v_num;
	unsigned int *first;
	Edge *e;

	int *in_deg, *out_deg;
	int max_in_deg, max_out_deg;
} Graph;


void clear_graph(Graph *g);

Graph *init_graph(int vertex_num, unsigned int edge_num);

void free_graph(Graph *g);

void add_edge(Graph *g, int src, int dst);

#endif // GRAPH_H
