/**
 * File: graph.c
 * Author: Cui Donghang
 * Last Update: 2022/04/24
 * This file implements all methods for struct Graph declared in graph.h.
 */
#include <string.h>
#include <stdlib.h>
#include "graph.h"

/**
 * Clear all data in graph.
 * Note that This function will not free this graph.
 */
void clear_graph(Graph *g) {
	g->e_cnt = g->max_in_deg = g->max_out_deg = 0;
	memset(g->first, -1, sizeof(unsigned int) * (g->v_num + 1));
	memset(g->in_deg, 0, sizeof(int) * (g->v_num + 1));
	memset(g->out_deg, 0, sizeof(int) * (g->v_num + 1));
}

/**
 * Initialize a graph with specified |V| and |E|.
 * Params:
 * * vertex_num: the max vertex id in V
 * * edge_num: the number |E|
 * Return the pointer of the initialized graph.
 * Note that the max vertex id in V is not always identical to |V|!
 */
Graph *init_graph(int vertex_num, unsigned int edge_num) {
	Graph *g = (Graph *)malloc(sizeof(Graph));
	g->e_num = edge_num;
	g->v_num = vertex_num;
	g->first = (unsigned int *)malloc(sizeof(unsigned int) * (vertex_num + 1));
	g->e = (Edge *)malloc(sizeof(Edge) * (edge_num + 10));
	g->in_deg = (int *)malloc(sizeof(int) * (vertex_num + 1));
	g->out_deg = (int *)malloc(sizeof(int) * (vertex_num + 1));
	clear_graph(g);
	return g;
}

/**
 * Free all the memory used by a graph.
 */
void free_graph(Graph *g) {
	free(g->first);
	free(g->e);
	free(g->in_deg);
	free(g->out_deg);
	free(g);
}

/**
 * add a DIRECTED edge to a graph with weight.
 * * Params:
 * * g: the pointer of a graph
 * * src: source vertex id of the DIRECTED edge
 * * dst: destination vertex id of the DIRECTED edge
 * * w: weight of the DIRECTED edge
 */
void add_edge(Graph *g, int src, int dst) {
	g->e[++g->e_cnt] = (Edge) {src, dst, g->first[src]};
	g->first[src] = g->e_cnt;	
	g->in_deg[dst]++; 
	g->out_deg[src]++;
	g->max_in_deg = max(g->in_deg[dst], g->max_in_deg);
	g->max_out_deg = max(g->out_deg[src], g->max_out_deg);
}

