### Distributed K-clique Counting

---

本科毕业设计《分布式k-团计数算法的设计与实现》源代码

基于Pivoter算法实现，具体参阅：Jain, Shweta, and C. Seshadhri. "The power of pivoting for exact clique counting." Proceedings of the 13th international conference on web search and data mining. 2020.

编译和运行方式：

```bash
make
bash ./run_dkc.sh <params...>
```

清除编译内容：

```bash
make clean
```



查看帮助：

```bash
bash ./run_dkc.sh -h
```

指定节点数量和数据集，在本地环境下使用多进程运行：

```bash
bash ./run_dkc.sh -n <num_nodes> -d <dataset>
```

指定hostfile和数据集，在多机集群上运行，其中的hostfile参照MPI的hostfile格式：

```bash
bash ./run_dkc.sh -f <hostfile> -d <dataset>
```

指定hostfile、数据集、主节点并行数量（可选，默认为4）从节点并行数量（可选，默认为2）、任务划分等级（可选，默认为1）、k-团的最大k值（可选，默认为0x7fffffff）、任务组打包任务数量（可选，默认为64）、任务缓冲队列大小（可选，默认为1024）：

```bash
bash ./run_dkc.sh \
	-f <hostfile> \
	-d <dataset> \
	-m <num_threads_master_node> \
    -w <num_threads_worker_node> \
    -l <division_level> \
    -k <max_k> \
    -t <taskgroup_size> \
    -b <buffer_size>
```

命令中的`<dataset>`格式如下，其中第一行为顶点和边的数量，第二行开始每一行均为**无向边**的两个顶点编号：

```
<number of vertice> <number of edges>
<vertex 1> <vertex 2>
<vertex 3> <vertex 4>
...
```

在`data/`目录下提供了两个格式化的数据集供测试使用。
