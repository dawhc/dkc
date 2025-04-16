#include <algorithm>
#include <cstdio>

struct Edge {
	int u, v;
	bool operator < (const Edge &e) const {
		return u == e.u ? v < e.v : u < e.u;
	}
} e[100000000];

inline int max(int a, int b) {
	return a > b ? a : b;
}

int main(int argc, char **argv) {
	freopen(argv[1], "r", stdin);
	freopen(argv[2], "w", stdout);
	int n, m; scanf("%d%d", &n, &m);
	int max_v = 0;
	for (int i = 1; i <= m; ++i) {
		int u, v; scanf("%d%d", &u, &v);
		if (u > v) std::swap(u, v);
		e[i].u = u;
		e[i].v = v;
	}
	std::sort(e + 1, e + m + 1);
	e[0].u = e[0].v = -1;
	int num_edges = 0;
	for (int i = 1; i <= m; ++i) {
		if (e[i].u == e[i - 1].u && e[i].v == e[i - 1].v) continue;
		e[++num_edges].u = e[i].u;
		e[num_edges].v = e[i].v;
		max_v = max(max_v, max(e[i].u, e[i].v));
	}
	printf("%d\t%d\n", max_v + 1, num_edges);
	for (int i = 1; i <= num_edges; ++i)
		printf("%d\t%d\n", e[i].u, e[i].v);
	return 0;
}
