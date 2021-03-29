#include "loi_maximal_clique.hpp"
LoiMaximalClique::LoiMaximalClique()
{
    v_num = 0;
    e_num = 0;
    align_malloc((void **)&pool_edges, 32, sizeof(int) * PACK_NODE_POOL_SIZE);
    align_malloc((void **)&pool_mc, 32, sizeof(int) * PACK_NODE_POOL_SIZE);
}

LoiMaximalClique::~LoiMaximalClique()
{
    free(pool_edges);
    free(pool_mc);
}

void LoiMaximalClique::build(const EdgeVector &_e_v)
{
    edge_vec.reserve(_e_v.size());
    for (auto &e : _e_v)
        if (e.first != e.second)
        {
            edge_vec.push_back(e);
            v_num = std::max(v_num, e.first);
            v_num = std::max(v_num, e.second);
        }

    std::sort(edge_vec.begin(), edge_vec.end(), edge_idpair_cmp);
    edge_vec.erase(std::unique(edge_vec.begin(), edge_vec.end()), edge_vec.end());
    v_num++;
    e_num = (long long)edge_vec.size();

    graph.resize(v_num);

    int cur_node_idx = 0;
    int prev_u = -1;
    for (auto &e : edge_vec)
    {
        if (e.first != prev_u)
        {
            prev_u = e.first;
            graph[e.first].start = cur_node_idx;
        }
        max_deg = std::max(max_deg, ++graph[e.first].deg);
        pool_edges[cur_node_idx++] = e.second;
        if (e.first > e.second)
        {
            graph[e.first].offset++;
        }
    }
    printf("v_num=%d e_num=%lld max_deg=%d\n", v_num, e_num, max_deg);
}

int LoiMaximalClique::build_matrix(const QVertex &a)
{
    int triangles = 0;
    root_start = a.start;
    root_offset = a.offset;
    root_deg = a.deg;
    root_vector_size = get_vector_size(a.deg);
    memset(matrix, 0, a.deg * root_vector_size * sizeof(Bitmap));
    for (int b_idx = 0; b_idx < a.deg; b_idx++)
    {
        const int b_id = pool_edges[a.start + b_idx];
        const QVertex &b = graph[b_id];
        triangles += mark_intersect(pool_edges + a.start, a.deg,
                                    pool_edges + b.start, b.deg, &matrix[b_idx * root_vector_size]);
    }
    return triangles;
}

int LoiMaximalClique::maximal_clique_bk()
{

    align_malloc((void **)&pool_sets, 32, sizeof(int) * PACK_NODE_POOL_SIZE);
    int mem_idx = 0;
    max_pool_sets_idx = 0;
    maximum_clique_size = 0;
    pool_mc_idx = 0;
    mc_num = 0;
    intersect_call_time = 0;
    big_intersect_call_time = 0;

    for (int i = 0; i < v_num; ++i)
        pool_sets[mem_idx++] = i;
    QVertex P(0, v_num), X(mem_idx, 0);
    std::vector<int> R;
    R.reserve(2048);

    BronKerbosch(R, P, X, mem_idx);
    printf("max_pool_sets_idx=%d\n", max_pool_sets_idx);
    printf("maximum_clique_size=%d\n", maximum_clique_size);
    printf("intersect_call_ratio=%.3f%%(%d/%d)\n",
           big_intersect_call_time * 100.0 / intersect_call_time,
           big_intersect_call_time, intersect_call_time);

    free(pool_sets);
    return mc_num;
}

int LoiMaximalClique::maximal_clique_degen()
{
    align_malloc((void **)&pool_sets, 32, sizeof(int) * PACK_NODE_POOL_SIZE);
    align_malloc((void **)&temp_set, 32, sizeof(int) * v_num);
    max_pool_sets_idx = 0;
    maximum_clique_size = 0;
    pool_mc_idx = 0;
    mc_num = 0;
    intersect_call_time = 0;
    big_intersect_call_time = 0;

    // std::unordered_set<int> visited, unvisited;
    // for (int i = 0; i < v_num; ++i) unvisited.insert(i);
    bool *visited = new bool[v_num];
    memset(visited, 0, v_num * sizeof(bool));
    std::vector<int> dorder = degeneracy_order();

    std::vector<int> R;
    R.reserve(2048);
    R.push_back(-1);

    for (auto v : dorder)
    {
        // unvisited.erase(v);
        R[0] = v;
        QVertex P(0, 0);
        for (int i = 0; i < graph[v].deg; ++i)
        {
            int u = pool_edges[graph[v].start + i];
            // if (unvisited.find(u) != unvisited.end()) {
            if (visited[u] == 0)
            {
                pool_sets[P.start + P.deg] = u;
                P.deg++;
            }
        }
        QVertex X(P.deg, 0);
        for (int i = 0; i < graph[v].deg; ++i)
        {
            int u = pool_edges[graph[v].start + i];
            // if (visited.find(u) != visited.end()) {
            if (visited[u] == 1)
            {
                pool_sets[X.start + X.deg] = u;
                X.deg++;
            }
        }

        Tomita(R, P, X);
        // visited.insert(v);
        visited[v] = 1;
    }
    R.pop_back();

    printf("max_pool_sets_idx=%d\n", max_pool_sets_idx);
    printf("maximum_clique_size=%d\n", maximum_clique_size);
    // printf("intersect_call_ratio=%.3f%%(%d/%d)\n",
    //     big_intersect_call_time * 100.0 / intersect_call_time,
    //     big_intersect_call_time, intersect_call_time);

    free(pool_sets);
    free(temp_set);
    delete[] visited;

    return mc_num;
}
int LoiMaximalClique::maximal_clique_loi()
{

#define LOG(x) std::cout << x << std::endl
    int buffer_vector_size = get_vector_size(max_deg);
    matrix = (Bitmap *)calloc(max_deg * buffer_vector_size, sizeof(Bitmap));
    clique_pool = (Bitmap *)calloc((max_deg + 1) * buffer_vector_size, sizeof(Bitmap));
    index_pool = (int *)calloc(max_deg * buffer_vector_size, sizeof(int));
    id_vec = (int *)calloc(max_deg + 1, sizeof(int));
    index_vec = (int *)calloc(max_deg, sizeof(int));
    index_vec[0] = -1;
    mc_num = 0;
    u_id = 0;
    for (QVertex &u : graph)
    {
        id_vec[0] = u_id++;
        build_matrix(u);
        memset(clique_pool, 0xff, root_vector_size * sizeof(Bitmap));
        for (int v_index = u.offset; v_index < u.deg; v_index++)
        {
            dfs_clz(v_index, 1);
        }
    }

    free(matrix);
    free(clique_pool);
    free(id_vec);
    free(index_vec);
    free(index_pool);

    printf("max_pool_sets_idx=%d\n", max_pool_sets_idx);
    printf("maximum_clique_size=%d\n", maximum_clique_size);
    return mc_num;
}

void LoiMaximalClique::dfs_avx2(int v_index, int depth)
{
    // check if v is part of a maximal clique
    index_vec[depth] = v_index;
    id_vec[depth] = pool_edges[root_start + v_index];
    intersect(get_bitmap(v_index), get_clique(depth - 1), get_clique(depth), root_vector_size);
    if (all_zero(get_clique(depth), root_vector_size))
    {
        mc_num++;
        if (pool_mc_idx + depth + 1 < PACK_NODE_POOL_SIZE)
        {
            memcpy(pool_mc + pool_mc_idx, id_vec, (depth + 1) * sizeof(int));
            pool_mc_idx += depth + 1;
            pool_mc[pool_mc_idx++] = -1;
        }

        // max_pool_sets_idx = std::max(max_pool_sets_idx, depth + 1);
        maximum_clique_size = std::max(maximum_clique_size, depth + 1);
        // if (root_start == 0)
        // {
        //     LOG("id: " << id_vec[0] << " offset: " << root_offset << " deg: " << root_deg << " vector_size: " << root_vector_size);
        //     LOG("pool_edges: " << to_string(pool_edges + root_start, root_deg));
        //     LOG("matrix: \n"
        //         << matrix_to_string());
        //     LOG("clique_pool: " << bitmap_to_string(get_clique(depth), root_vector_size));
        //     LOG("index_vec: " << to_string(index_vec, depth + 1));
        //     LOG("id_vec: " << to_string(id_vec, depth + 1) << "\n");
        // }
        return;
    }
    int num = expandToIndex(get_clique(depth), get_index_by_depth(depth), v_index + 1, root_deg);
    int *next = get_index_by_depth(depth);
    for (int i = 0; i < num; i++)
    {
        dfs_avx2(next[i], depth + 1);
    }
}

void LoiMaximalClique::dfs_clz(int v_index, int depth) {
    // check if v is part of a maximal clique
    index_vec[depth] = v_index;
    id_vec[depth] = pool_edges[root_start + v_index];
    bool found = intersect_allzero(get_bitmap(v_index), get_clique(depth - 1), get_clique(depth), root_vector_size);
    if (found)
    {
        mc_num++;
        if (pool_mc_idx + depth + 1 < PACK_NODE_POOL_SIZE)
        {
            memcpy(pool_mc + pool_mc_idx, id_vec, (depth + 1) * sizeof(int));
            pool_mc_idx += depth + 1;
            pool_mc[pool_mc_idx++] = -1;
        }
        // if (root_start == 0)
        // {
        //     LOG("id: " << id_vec[0] << " offset: " << root_offset << " deg: " << root_deg << " vector_size: " << root_vector_size);
        //     LOG("pool_edges: " << to_string(pool_edges + root_start, root_deg));
        //     LOG("matrix: \n"
        //         << matrix_to_string());
        //     LOG("clique_pool: \n" << bitmap_to_string(get_clique(depth), root_vector_size));
        //     LOG("index_vec: " << to_string(index_vec, depth + 1));
        //     LOG("id_vec: " << to_string(id_vec, depth + 1) << "\n");
        // }
        return;
    }
    if(v_index + 1 == root_deg) {
        return;
    }
    Bitmap *bitvec = get_clique(depth);
    int v_filter = v_index + 1;
    int v_start = v_filter / (sizeof(Bitmap) * 8);
    Bitmap bitset = bitvec[v_start];
    bitset = bitset & (BITMASK >> (v_filter % 64));
    while (bitset)
    {
        int r = __builtin_clzll(bitset);
        int new_index = r + 64 * v_start;
        bitset ^= 1LLU << (63 - r);
        dfs_clz(new_index, depth + 1);
    }
    // COUNT LEADING ZERO
    for (int i = v_start + 1; i < root_vector_size; i++)
    {
        bitset = bitvec[i];
        while (bitset)
        {
            int r = __builtin_clzll(bitset);
            int new_index = r + 64 * i;
            bitset ^= 1LLU << (63 - r);
            dfs_clz(new_index, depth+1);
        }
    }
}
void LoiMaximalClique::report_mc_num()
{
    int counter = 0;
    for (;;)
    {
        std::this_thread::sleep_for(std::chrono::seconds(REPORT_ELAPSE));
        counter++;
        if (counter > MAX_REPORT_TIME)
        {
            break;
        }
        std::cout << counter << " seconds: " << mc_num << " vertex processesd: " << u_id << std::endl;
    }
}

void LoiMaximalClique::BronKerbosch(std::vector<int> &R, QVertex P, QVertex X, int mem_idx)
{
    if (P.deg == 0)
    { // report R as a maximal clique.
        if (X.deg == 0)
        {
            std::sort(R.begin(), R.end());
            memcpy(pool_mc + pool_mc_idx, R.data(), R.size() * sizeof(int));
            pool_mc_idx += R.size();
            pool_mc[pool_mc_idx++] = -1; // split different maximal cliques by -1.
            mc_num++;
            max_pool_sets_idx = std::max(max_pool_sets_idx, mem_idx);
            maximum_clique_size = std::max(maximum_clique_size, (int)R.size());
        }
        return;
    }

    R.push_back(-1); // enlarge R by a placeholder(-1).
    for (int i = 0; i < P.deg; ++i)
    {
        int v = pool_sets[P.start + i];
        R.back() = v;
        int newPstart = mem_idx;
        int newPdeg = intersect(pool_sets + P.start + i, P.deg - i,
                                pool_edges + graph[v].start, graph[v].deg, pool_sets + newPstart);
        QVertex newP(newPstart, newPdeg);
        int newXstart = newPstart + newPdeg;
        int newXdeg1 = intersect(pool_sets + X.start, X.deg,
                                 pool_edges + graph[v].start, graph[v].deg, pool_sets + newXstart);
        int newXdeg2 = intersect(pool_sets + P.start, i,
                                 pool_edges + graph[v].start, graph[v].deg, pool_sets + newXstart + newXdeg1);
        QVertex newX(newXstart, newXdeg1 + newXdeg2);
        // for (int j = 1; j < newX.deg; ++j)
        //     assert(pool_sets[newX.start + j] > pool_sets[newX.start + j - 1]);
        BronKerbosch(R, newP, newX, newX.start + newX.deg);
    }
    R.pop_back();
}

int LoiMaximalClique::maximal_clique_pivot()
{
    align_malloc((void **)&pool_sets, 32, sizeof(int) * PACK_NODE_POOL_SIZE);
    align_malloc((void **)&temp_set, 32, sizeof(int) * v_num);
    max_pool_sets_idx = 0;
    maximum_clique_size = 0;
    pool_mc_idx = 0;
    mc_num = 0;
    intersect_call_time = 0;
    big_intersect_call_time = 0;

    for (int i = 0; i < v_num; ++i)
        pool_sets[i] = i;
    QVertex P(0, v_num), X(v_num, 0);
    std::vector<int> R;
    R.reserve(2048);

    Tomita(R, P, X);

    printf("max_pool_sets_idx=%d\n", max_pool_sets_idx);
    printf("maximum_clique_size=%d\n", maximum_clique_size);
    printf("intersect_call_ratio=%.3f%%(%d/%d)\n",
           big_intersect_call_time * 100.0 / intersect_call_time,
           big_intersect_call_time, intersect_call_time);

    free(pool_sets);
    free(temp_set);

    return mc_num;
}

void LoiMaximalClique::Tomita(std::vector<int> &R, QVertex P, QVertex X)
{
    if (P.deg == 0)
    { // report R as a maximal clique.
        if (X.deg == 0)
        {
            std::sort(R.begin(), R.end());
            memcpy(pool_mc + pool_mc_idx, R.data(), R.size() * sizeof(int));
            pool_mc_idx += R.size();
            pool_mc[pool_mc_idx++] = -1; // split different maximal cliques by -1.
            mc_num++;
            max_pool_sets_idx = std::max(max_pool_sets_idx, X.start + X.deg);
            maximum_clique_size = std::max(maximum_clique_size, (int)R.size());
        }
        return;
    }

    // heuristic 1: choose a pivot u from P or X to maximize |P and N(u)|.
    // int u = -1, max_inter_cnt = -1;
    // for (int i = 0; i < P.deg; ++i) {
    //     int cu = pool_sets[P.start + i];
    //     int inter_cnt = intersect_count(pool_sets + P.start, P.deg,
    //             pool_edges + graph[cu].start, graph[cu].deg);
    //     if (max_inter_cnt < inter_cnt) {
    //         max_inter_cnt = inter_cnt;
    //         u = cu;
    //     }
    // }
    // for (int i = 0; i < X.deg; ++i) {
    //     int cu = pool_sets[X.start + i];
    //     int inter_cnt = intersect_count(pool_sets + X.start, X.deg,
    //             pool_edges + graph[cu].start, graph[cu].deg);
    //     if (max_inter_cnt < inter_cnt) {
    //         max_inter_cnt = inter_cnt;
    //         u = cu;
    //     }
    // }

    // heuristic 2: choose a pivot u from P or X to maximize |N(u)|.
    // int u = -1, max_deg = -1;
    // for (int i = 0; i < P.deg; ++i) {
    //     int cu = pool_sets[P.start + i];
    //     if (max_deg < graph[cu].deg) {
    //         max_deg = graph[cu].deg;
    //         u = cu;
    //     }
    // }
    // for (int i = 0; i < X.deg; ++i) {
    //     int cu = pool_sets[X.start + i];
    //     if (max_deg < graph[cu].deg) {
    //         max_deg = graph[cu].deg;
    //         u = cu;
    //     }
    // }

    // heuristic 3: choose the first vertex from P as the pivot.
    int u;
    if (X.deg > 0)
        u = pool_sets[X.start];
    else
        u = pool_sets[P.start];

    // only need to enumerate verices in P - N(u).
    int *N_u_ptr = pool_edges + graph[u].start;
    int *N_u_end = pool_edges + graph[u].start + graph[u].deg;
    int Pbound = 0;

    R.push_back(-1); // enlarge R by a placeholder(-1).
    for (int i = 0; i < P.deg; ++i)
    {
        int v = pool_sets[P.start + i];
        while (N_u_ptr != N_u_end && *N_u_ptr < v)
            N_u_ptr++;
        if (N_u_ptr != N_u_end && *N_u_ptr == v)
        { // skip vertices in N(u).
            N_u_ptr++;
            continue;
        }
        R.back() = v;

        int newPstart = X.start + X.deg;

        // intersect_cnt += 2;
        // gettimeofday(&time_start, NULL);

#if SIMD_STATE == 2
        int newPdeg = intersect_scalar2x(pool_sets + P.start + Pbound, P.deg - Pbound,
                                         pool_edges + graph[v].start, graph[v].deg, pool_sets + newPstart);
#elif SIMD_STATE == 4
        int newPdeg = intersect_simd4x(pool_sets + P.start + Pbound, P.deg - Pbound,
                                       pool_edges + graph[v].start, graph[v].deg, pool_sets + newPstart);
#else
        int newPdeg = intersect(pool_sets + P.start + Pbound, P.deg - Pbound,
                                pool_edges + graph[v].start, graph[v].deg, pool_sets + newPstart);
#endif
        QVertex newP(newPstart, newPdeg);
        int newXstart = newPstart + newPdeg;

        int mergeXdeg = merge(pool_sets + P.start, Pbound,
                              pool_sets + X.start, X.deg, temp_set);
#if SIMD_STATE == 2
        int newXdeg = intersect_scalar2x(temp_set, mergeXdeg,
                                         pool_edges + graph[v].start, graph[v].deg, pool_sets + newXstart);
#elif SIMD_STATE == 4
        int newXdeg = intersect_simd4x(temp_set, mergeXdeg,
                                       pool_edges + graph[v].start, graph[v].deg, pool_sets + newXstart);
#else
        int newXdeg = intersect(temp_set, mergeXdeg,
                                pool_edges + graph[v].start, graph[v].deg, pool_sets + newXstart);
#endif
        // gettimeofday(&time_end, NULL);
        // intersect_time += (time_end.tv_sec - time_start.tv_sec) * 1000.0 + (time_end.tv_usec - time_start.tv_usec) / 1000.0;

        QVertex newX(newXstart, newXdeg);

        Tomita(R, newP, newX);

        // move v from P to X.
        for (int j = i; j > Pbound; --j)
            pool_sets[P.start + j] = pool_sets[P.start + j - 1];
        pool_sets[P.start + Pbound] = v;
        Pbound++;
    }
    R.pop_back();
}

std::vector<int> LoiMaximalClique::degeneracy_order()
{
    std::vector<int> deg(v_num);
    int md = 0;
    for (int i = 0; i < v_num; ++i)
    {
        deg[i] = graph[i].deg;
        md = std::max(md, deg[i]);
    }

    std::vector<int> bin(md + 1);
    for (int i = 0; i <= md; ++i)
        bin[i] = 0;
    for (int i = 0; i < v_num; ++i)
        bin[deg[i]]++;

    int start = 0;
    for (int i = 0; i <= md; ++i)
    {
        int num = bin[i];
        bin[i] = start;
        start += num;
    }

    std::vector<int> vert(v_num), pos(v_num);
    for (int i = 0; i < v_num; ++i)
    {
        pos[i] = bin[deg[i]];
        vert[pos[i]] = i;
        bin[deg[i]]++;
    }
    for (int i = md; i > 0; --i)
        bin[i] = bin[i - 1];
    bin[0] = 0;

    std::vector<int> degen_order;
    degen_order.reserve(v_num);
    int degeneracy = 0;
    for (int i = 0; i < v_num; ++i)
    {
        int v = vert[i];
        degen_order.push_back(v);
        degeneracy = std::max(degeneracy, deg[v]);
        for (int j = 0; j < graph[v].deg; ++j)
        {
            int u = pool_edges[graph[v].start + j];
            if (deg[u] > deg[v])
            {
                int du = deg[u], pu = pos[u];
                int pw = bin[du], w = vert[pw];
                if (u != w)
                {
                    pos[u] = pw;
                    vert[pu] = w;
                    pos[w] = pu;
                    vert[pw] = u;
                }
                bin[du]++;
                deg[u]--;
            }
        }
    }

    printf("degeneracy=%d\n", degeneracy);

    return degen_order;
}

void LoiMaximalClique::save_answers(const char *file_path)
{
    FILE *fp = fopen(file_path, "w");
    if (fp == NULL)
    {
        std::cout << "fail to create " << file_path << std::endl;
        quit();
    }

    for (int i = 0; i < pool_mc_idx; ++i)
        if (pool_mc[i] == -1)
            fprintf(fp, "\n");
        else
            fprintf(fp, "%d ", pool_mc[i]);

    fclose(fp);
}
void LoiMaximalClique::start_report()
{
    std::thread(&LoiMaximalClique::report_mc_num, this).detach();
};
