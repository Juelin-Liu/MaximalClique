#include "loi_maximal_clique.hpp"
#define LOG(x) std::cout << x << std::endl
#define PADDING 16
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

void LoiMaximalClique::set_buffer_capacity(size_t deg)
{
    size_t aligned_buffer_vector_size = get_vector_size<AlignType>(max_deg);
    align_malloc((void **)&matrix, sizeof(AlignType), (deg + 1) * aligned_buffer_vector_size);
    align_malloc((void **)&P_vec_pool, sizeof(AlignType), (deg + 1) * aligned_buffer_vector_size);
    align_malloc((void **)&X_vec_pool, sizeof(AlignType), (deg + 1) * aligned_buffer_vector_size);
    align_malloc((void **)&X_init, sizeof(AlignType), aligned_buffer_vector_size);
    align_malloc((void **)&stack_pool, sizeof(AlignType), (deg + 1) * aligned_buffer_vector_size);
    align_malloc((void **)&simd_buffer, sizeof(AlignType), sizeof(AlignType));

    R = (int *)malloc(sizeof(int) * deg);
    X = (int *)malloc(sizeof(int) * deg);
    stack_set_pool = (int *)malloc(sizeof(int) * deg * deg);
    pool_sets = (int *)malloc(sizeof(int) * deg * deg);
    index_vec = (int *)malloc(sizeof(int) * (deg + 1));
    stack_set_size = (int *)malloc(sizeof(int) * (deg + 1));

    triangle_cnt = (int *)malloc(sizeof(int) * (deg + 1));
    pivot_inter_cnt = (int *)malloc(sizeof(int) * (deg + 1));
};

void LoiMaximalClique::free_buffer()
{
    if (matrix != NULL)
    {
        free(matrix);
    }
    if (P_vec_pool != NULL)
    {
        free(P_vec_pool);
    }
    if (X_vec_pool != NULL)
    {
        free(X_vec_pool);
    }
    if (stack_set_pool != NULL)
    {
        free(stack_set_pool);
    }
    if (index_vec != NULL)
    {
        free(index_vec);
    }
    if (R != NULL)
    {
        free(R);
    }
    if (X != NULL)
    {
        free(X);
    }
    if (triangle_cnt != NULL)
    {
        free(triangle_cnt);
    }
}
void LoiMaximalClique::init_R_X(const QVertex &a)
{
    aligned_root_vector_size = get_vector_size<AlignType>(a.deg);
    memset(get_xvec(0), 0, a.deg * aligned_root_vector_size * sizeof(Bitmap));
    memset(get_pvec(0), 0, a.deg * aligned_root_vector_size * sizeof(Bitmap));
    for (int b_idx = 0; b_idx < a.deg; b_idx++)
    {
        const int b_id = pool_edges[a.start + b_idx];
        if (visited[b_id])
        {
            mark_as_one(get_xvec(0), b_idx);
        }
        // else
        // {
        //     mark_as_one(get_pvec(0), b_idx);
        // }
    }
    fill_with_one(get_pvec(0), a.deg);
    bitwise_andn(get_pvec(0), get_xvec(0), get_pvec(0), aligned_root_vector_size);
}

void LoiMaximalClique::build_matrix(int u_id)
{
    const auto &a = graph[u_id];
    int triangles = 0;
    int p_index = 0, p_triangle_size = 0, max_inter_cnt = 0;
    root_triangle_cnt = 0;
    root_start = a.start;
    root_offset = a.offset;
    root_deg = a.deg;
    root_vector_size = get_vector_size<Bitmap>(a.deg);
    aligned_root_vector_size = get_vector_size<AlignType>(a.deg);
    memset(get_bitmap(0), 0, a.deg * aligned_root_vector_size * sizeof(Bitmap));
    memset(triangle_cnt, 0, a.deg * sizeof(int));
    memset(pivot_inter_cnt, 0, a.deg * sizeof(int));
    memset(stack_set_size, 0, a.deg * sizeof(int));
    memset(get_xvec(0), 0, a.deg * aligned_root_vector_size * sizeof(Bitmap));
    memset(get_pvec(0), 0, a.deg * aligned_root_vector_size * sizeof(Bitmap));
    for (int b_idx = 0; b_idx < a.deg; b_idx++)
    {
        const int b_id = pool_edges[a.start + b_idx];
        if (!visited[b_id])
        {
            mark_as_one(get_pvec(0), b_idx);
            const QVertex &b = graph[b_id];
            int t_num = mark_intersect_simd8x(pool_edges + a.start, a.deg,
                                              pool_edges + b.start, b.deg, get_bitmap(b_idx));
            root_triangle_cnt += t_num;
        }
        else
        {
            mark_as_one(get_xvec(0), b_idx);
        }
    }
    memcpy(X_init, get_xvec(0), aligned_root_vector_size);
}

int LoiMaximalClique::build_matrix(const QVertex &a)
{
    int triangles = 0;
    int p_index = 0, p_triangle_size = 0, max_inter_cnt = 0;
    root_triangle_cnt = 0;
    root_start = a.start;
    root_offset = a.offset;
    root_deg = a.deg;
    root_vector_size = get_vector_size<Bitmap>(a.deg);
    aligned_root_vector_size = get_vector_size<AlignType>(a.deg);
    memset(get_bitmap(0), 0, a.deg * aligned_root_vector_size * sizeof(Bitmap));
    memset(triangle_cnt, 0, a.deg * sizeof(int));
    memset(pivot_inter_cnt, 0, a.deg * sizeof(int));
    memset(stack_set_size, 0, a.deg * sizeof(int));
    for (int b_idx = 0; b_idx < a.deg; b_idx++)
    {
        const int b_id = pool_edges[a.start + b_idx];
        if (!visited[b_id])
        {
            const QVertex &b = graph[b_id];
            int t_num = mark_intersect_simd8x(pool_edges + a.start, a.deg,
                                              pool_edges + b.start, b.deg, get_bitmap(b_idx));
            //herustic 2 choose v with maximum triangle
            triangle_cnt[b_idx] = t_num;
            if (t_num > p_triangle_size)
            {
                p_triangle_size = t_num;
                p_index = b_idx;
            }
            root_triangle_cnt += t_num;

            // herustic 3 choose v that minimize P \ N(v)
            int b_intern_cnt = bitwise_and_count(get_bitmap(b_idx), get_pvec(0), root_vector_size);
            if (b_intern_cnt > max_inter_cnt)
            {
                p_index = b_idx;
                max_inter_cnt = b_intern_cnt;
            }
        }
    }
    // herustic 3 choose v that minimize P \ N(v)
    pivot_inter_cnt[0] = max_inter_cnt;
    return p_index;
}

long long LoiMaximalClique::maximal_clique_bk()
{
    set_buffer_capacity(max_deg);
    index_vec[0] = -1;
    mc_num = 0, u_cnt = 0, total_mc_size = 0;
    visited = new bool[v_num];
    for (int u_id = 0; u_id < v_num; u_id++)
    {
        R[0] = u_id;
        build_matrix(u_id);
        BronKerbosch(0, 0);
        visited[u_id] = true;
        u_cnt = u_id;
    }
    delete[] visited;
    free_buffer();
    // printf("max_pool_sets_idx=%d\n", max_pool_sets_idx);
    printf("maximum_clique_size=%d\n", maximum_clique_size);
    return mc_num;
}

long long LoiMaximalClique::maximal_clique_pivot()
{
    set_buffer_capacity(max_deg);
    index_vec[0] = -1;
    mc_num = 0, u_cnt = 0, total_mc_size = 0;
    visited = new bool[v_num];
    for (int u_id = 0; u_id < v_num; u_id++)
    {
        R[0] = u_id;
        build_matrix(u_id);
        Tomita(0, 0);
        visited[u_id] = true;
        u_cnt = u_id;
    }
    delete[] visited;
    free_buffer();
    // printf("max_pool_sets_idx=%d\n", max_pool_sets_idx);
    printf("maximum_clique_size=%d\n", maximum_clique_size);
    return mc_num;
}

long long LoiMaximalClique::maximal_clique_degen()
{
    // return maximal_clique_degen_onepunch();
    set_buffer_capacity(max_deg);
    index_vec[0] = -1;
    mc_num = 0, u_cnt = 0, total_mc_size = 0;
    visited = new bool[v_num];
    auto dorder = degeneracy_order();
    for (int u_id : dorder)
    {
        R[0] = u_id;
        build_matrix(u_id);
        Tomita(0, 0);
        visited[u_id] = true;
        u_cnt++;
    }
    delete[] visited;
    free_buffer();
    // printf("max_pool_sets_idx=%d\n", max_pool_sets_idx);
    printf("maximum_clique_size=%d\n", maximum_clique_size);
    return mc_num;
}

void LoiMaximalClique::BronKerbosch(int depth, int mem_idx)
{
    int num = expand_avx2_compress(get_pvec(depth), pool_sets + mem_idx, root_vector_size);
    for (int i = 0; i < num; i++)
    {
        int n_index = *(pool_sets + mem_idx + i);
        R[depth + 1] = pool_edges[root_start + n_index];
        int psize = bitwise_and_count(get_pvec(depth), get_bitmap(n_index), get_pvec(depth + 1), root_vector_size);
        int xsize = bitwise_and_count(get_xvec(depth), get_bitmap(n_index), get_xvec(depth + 1), root_vector_size);
        if (psize == 0)
        {
            if (xsize == 0)
            {
                mc_num++;
                total_mc_size += depth + 2;
                if (pool_mc_idx + depth + 2 < PACK_NODE_POOL_SIZE)
                {
                    memcpy(pool_mc + pool_mc_idx, R, (depth + 2) * sizeof(int));
                    pool_mc_idx += depth + 2;
                    pool_mc[pool_mc_idx++] = -1;
                }
            }
        }
        else
        {
            BronKerbosch(depth + 1, mem_idx + num);
        }
        mark_as_zero(get_pvec(depth), n_index);
        mark_as_one(get_xvec(depth), n_index);
    }
}

void LoiMaximalClique::Tomita(int depth, int mem_idx)
{
    // choose pivot vertex
    // herustic 3, choose v thath minimize P \ N(v)
    int num, pivot_index = -1;
    int *pnext = pool_sets + mem_idx;
    if (depth == 0 | pivot_inter_cnt[depth - 1] > 0)
    {
        int max_intersect_cnt = 0;
        // not considering x
        // int pnum = expand_avx2_compress(get_pvec(depth), get_stack(depth), root_vector_size);

        // considering x, X_init is ignored because all the indexes in the have been visited
        int pnum = maskzor_expand_avx2_compress(X_init, get_pvec(depth), get_xvec(depth), simd_buffer, pnext, root_vector_size);
        for (int i = 0; i < pnum; i++)
        {
            int p_idx = pnext[i];
            int p_cnt = bitwise_and_count(get_pvec(depth), get_bitmap(p_idx), root_vector_size);
            if (p_cnt >= max_intersect_cnt)
            {
                max_intersect_cnt = p_cnt;
                pivot_index = p_idx;
            }
        }
        pivot_inter_cnt[depth] = max_intersect_cnt;
        // overwrite next index
        num = maskz_expand_avx2_compress(get_bitmap(pivot_index), get_pvec(depth), simd_buffer, pnext, root_vector_size);
    }
    else
    {
        // overwrite next index
        num = expand_avx2_compress(get_pvec(depth), pnext, root_vector_size);
    }
    stack_set_size[depth] = num;
    cur_depth = depth + 1;
    for (int i = 0; i < num; i++)
    {
        int n_index = *(pool_sets + mem_idx + i);
        R[depth + 1] = pool_edges[root_start + n_index];
        int psize = bitwise_and_count(get_pvec(depth), get_bitmap(n_index), get_pvec(depth + 1), root_vector_size);
        int xsize = bitwise_and_count(get_xvec(depth), get_bitmap(n_index), get_xvec(depth + 1), root_vector_size);
        if (psize == 0)
        {
            if (xsize == 0)
            {
                mc_num++;
                total_mc_size += depth + 2;
                total_compressed_mc_size += root_vector_size;
                if (pool_mc_idx + depth + 2 < PACK_NODE_POOL_SIZE)
                {
                    memcpy(pool_mc + pool_mc_idx, R, (depth + 2) * sizeof(int));
                    pool_mc_idx += depth + 2;
                    pool_mc[pool_mc_idx++] = -1;
                }
            }
        }
        else
        {
            Tomita(depth + 1, mem_idx + num);
        }
        mark_as_zero(get_pvec(depth), n_index);
        mark_as_one(get_xvec(depth), n_index);
    }
}

// long long LoiMaximalClique::maximal_clique_bk()
// {
//     set_buffer_capacity(max_deg);
//     index_vec[0] = -1;
//     mc_num = 0, u_cnt = 0, p_set_idx = 0;
//     visited = new bool[v_num];
//     for (QVertex &u : graph)
//     {
//         R[0] = u_cnt;
//         init_R_X(u);
//         build_matrix(u);
//         int num = expand_avx2(get_nvec(0), get_stack(0), root_vector_size);
//         int *next = get_stack(0);
//         for (int i = 0; i < num; i++)
//         {
//             int v_index = next[i];
//             // int v_index = i;
//             dfs(v_index, 1);
//             mark_as_one(get_xvec(0), v_index);
//             mark_as_zero(get_pvec(0), v_index);
//         }
//         visited[u_cnt++] = true;
//     }
//     delete[] visited;
//     free_buffer();
//     // printf("max_pool_sets_idx=%d\n", max_pool_sets_idx);
//     printf("maximum_clique_size=%d\n", maximum_clique_size);
//     return mc_num;
// }

// TODO implement this
// long long LoiMaximalClique::maximal_clique_pivot()
// {
//     set_buffer_capacity(max_deg);
//     index_vec[0] = -1;
//     mc_num = 0, u_cnt = 0, p_set_idx = 0;
//     visited = new bool[v_num];
//     for (QVertex &u : graph)
//     {
//         R[0] = u_cnt;

//         if (visited[u_cnt])
//         {
//             u_cnt++;
//             continue;
//         }
//         init_R_X(u);
//         // choose the first vertex from X otherwise choose from P
//         int pivot_index = build_matrix(u);
//         // bitwise_andn(get_pvec(0), get_bitmap(pivot_index), get_nvec(0), aligned_root_vector_size);
//         // int num = expand_ctz(get_nvec(0), get_stack(0), root_vector_size);
//         int num = maskz_expand_avx2_compress(get_bitmap(pivot_index), get_pvec(0), simd_buffer, get_stack(0), root_vector_size);
//         int *next = get_stack(0);
//         for (int i = 0; i < num; i++)
//         {
//             int v_index = next[i];
//             // int v_index = i;
//             dfs_pivot(v_index, 1);
//             mark_as_one(get_xvec(0), v_index);
//             mark_as_zero(get_pvec(0), v_index);
//         }
//         visited[u_cnt++] = true;
//     }
//     delete[] visited;
//     free_buffer();
//     // printf("max_pool_sets_idx=%d\n", max_pool_sets_idx);
//     printf("maximum_clique_size=%d\n", maximum_clique_size);
//     return mc_num;
// }

// long long LoiMaximalClique::maximal_clique_degen()
// {
//     // return maximal_clique_degen_onepunch();
//     set_buffer_capacity(max_deg);
//     index_vec[0] = -1;
//     mc_num = 0, u_cnt = 0, p_set_idx = 0;
//     visited = new bool[v_num];
//     std::vector<int> dorder = degeneracy_order();
//     for (int u_id : dorder)
//     {
//         u_cnt++;
//         R[0] = u_id;
//         const QVertex &u = graph[u_id];
//         init_R_X(u);
//         // choose the v with max triangles
//         int pivot_index = build_matrix(u);
//         // bitwise_andn(get_pvec(0), get_bitmap(pivot_index), get_nvec(0), aligned_root_vector_size);
//         // int num = expand_ctz(get_nvec(0), get_stack(0), root_vector_size);
//         // int num = expand_avx2_compress(get_nvec(0), get_stack(0), root_vector_size);
//         int num = maskz_expand_avx2_compress(get_bitmap(pivot_index), get_pvec(0), simd_buffer, get_stack(0), root_vector_size);
//         int *next = get_stack(0);
//         for (int i = 0; i < num; i++)
//         {
//             int v_index = next[i];
//             // int v_index = i;
//             dfs_pivot(v_index, 1);
//             mark_as_one(get_xvec(0), v_index);
//             mark_as_zero(get_pvec(0), v_index);
//         }
//         visited[u_id] = true;
//     }
//     delete[] visited;
//     free_buffer();
//     // printf("max_pool_sets_idx=%d\n", max_pool_sets_idx);
//     printf("maximum_clique_size=%d\n", maximum_clique_size);
//     return mc_num;
// }

long long LoiMaximalClique::maximal_clique_degen_onepunch()
{
    set_buffer_capacity(max_deg);
    mc_num = 0, u_cnt = 0, p_set_idx = 0;
    visited = new bool[v_num];
    std::vector<int> dorder = degeneracy_order();
    for (int u_id : dorder)
    {
        u_cnt++;
        R[0] = u_id;
        const QVertex &u = graph[u_id];
        init_R_X(u);
        // choose the v with max triangles
        int pivot_index = build_matrix(u);

// bitwise_andn(get_pvec(0), get_bitmap(pivot_index), get_nvec(0), aligned_root_vector_size);
// int num = expand_ctz(get_nvec(0), get_stack(0), root_vector_size);
// int num = expand_avx2_compress(get_nvec(0), get_stack(0), root_vector_size);
#if __SIMD_LEVEL__ == 2
        int num = maskz_expand_avx2_compress(get_bitmap(pivot_index), get_pvec(0), simd_buffer, get_stack(0), root_vector_size);
#elif __SIMD_LEVEL__ == 0
        int num = maskz_expand_ctz(get_bitmap(pivot_index), get_pvec(0), get_stack(0), root_vector_size);
#endif
        stack_set_size[0] = num;
        dfs_pivot_serious_onepunch();
        int *next = get_stack(0);

        // for (int i = num - 1; i >= 0; i--)
        // {
        //     cur_depth = 0;
        //     int v_index = next[i];
        //     // int v_index = i;
        //     // dfs_pivot_normal_onepunch();
        //     mark_as_one(get_xvec(0), v_index);
        //     mark_as_zero(get_pvec(0), v_index);
        // }
        visited[u_id] = true;
    }
    delete[] visited;
    free_buffer();
    // printf("max_pool_sets_idx=%d\n", max_pool_sets_idx);
    printf("maximum_clique_size=%d\n", maximum_clique_size);
    return mc_num;
}

// void LoiMaximalClique::dfs(int v_index, int depth)
// {
//     // bookkeeping
//     R[depth] = pool_edges[root_start + v_index];
//     // compute the cliques formed with all visiting vertexes
//     bitwise_and(get_bitmap(v_index), get_pvec(depth - 1), get_pvec(depth), aligned_root_vector_size);
//     bitwise_and(get_bitmap(v_index), get_xvec(depth - 1), get_xvec(depth), aligned_root_vector_size);
//     if (all_zero(get_pvec(depth), root_vector_size))
//     {
//         // declare maximal clique is found
//         if (all_zero(get_xvec(depth), root_vector_size))
//         {
//             mc_num++;
//             if (pool_mc_idx + depth + 1 < PACK_NODE_POOL_SIZE)
//             {
//                 // std::sort(R, R + depth);
//                 memcpy(pool_mc + pool_mc_idx, R, (depth + 1) * sizeof(int));
//                 pool_mc_idx += depth + 1;
//                 pool_mc[pool_mc_idx++] = -1;
//             }
//             maximum_clique_size = std::max(maximum_clique_size, depth + 1);
//         }
//         return;
//     }
//     int num = expand_ctz(get_nvec(depth), get_stack(depth), root_vector_size);
//     int *next = get_stack(depth);
//     for (int i = 0; i < num; i++)
//     {
//         dfs(next[i], depth + 1);
//         mark_as_zero(get_pvec(depth), next[i]);
//         mark_as_one(get_xvec(depth), next[i]);
//     }
// }

// void LoiMaximalClique::dfs_pivot(int v_index, int depth)
// {
//     R[depth] = pool_edges[root_start + v_index];
//     // compute the cliques
//     bitwise_and(get_bitmap(v_index), get_pvec(depth - 1), get_pvec(depth), aligned_root_vector_size);
//     bitwise_and(get_bitmap(v_index), get_xvec(depth - 1), get_xvec(depth), aligned_root_vector_size);
//     if (all_zero(get_pvec(depth), root_vector_size))
//     {
//         // declare maximal clique is found
//         if (all_zero(get_xvec(depth), root_vector_size))
//         {
//             mc_num++;
//             total_mc_size += depth + 1;
//             if (pool_mc_idx + depth + 1 < PACK_NODE_POOL_SIZE)
//             {
//                 // std::sort(R, R + depth);
//                 memcpy(pool_mc + pool_mc_idx, R, (depth + 1) * sizeof(int));
//                 pool_mc_idx += depth + 1;
//                 pool_mc[pool_mc_idx++] = -1;
//             }
//             maximum_clique_size = std::max(maximum_clique_size, depth + 1);
//         }

//         // max_pool_sets_idx = std::max(max_pool_sets_idx, depth + 1);
//         return;
//     }
//     int num;
//     int *next = get_stack(depth);
//     // choose a pivot point
//     // herustic 1, choose v with max degree, (assume graph is reindexed, small id high deg)
//     // int pivot_index = find_first_index(get_pvec(depth), root_vector_size);

//     // herustic 2, choose v with most triangles
//     // int pivot_index = -1, max_cnt = 0;
//     // int pnum = expand_avx2_compress(get_pvec(depth), get_stack(depth), root_vector_size);
//     // int *pnext = get_stack(depth);
//     // for (int i = 0; i < pnum; i++)
//     // {
//     //     int p_idx = pnext[i];
//     //     if (triangle_cnt[p_idx] > max_cnt)
//     //     {
//     //         pivot_index = p_idx;
//     //         max_cnt = triangle_cnt[p_idx];
//     //     }
//     // }

//     // herustic 3, choose v thath minimize P \ N(v)
//     if (pivot_inter_cnt[depth - 1] > 0)
//     {
//         int pivot_index = -1, max_intersect_cnt = 0;
//         int *pnext = get_stack(depth);
//         // not considering x
//         // int pnum = expand_avx2_compress(get_pvec(depth), get_stack(depth), root_vector_size);

//         // considering x
//         // bitwise_or(get_pvec(depth), get_xvec(depth), get_nvec(depth), aligned_root_vector_size);
//         // bitwise_andn(get_nvec(depth), get_xvec(0), get_nvec(0), aligned_root_vector_size);
//         // int pnum = expand_avx2_compress(get_pvec(depth), get_stack(depth), root_vector_size);
//         int pnum = maskzor_expand_avx2_compress(get_xvec(0), get_pvec(depth), get_xvec(depth), simd_buffer, pnext, root_vector_size);
//         for (int i = 0; i < pnum; i++)
//         {
//             int p_idx = pnext[i];
//             // method 1 using avx2, minimum memory footprint
//             int p_cnt = bitwise_and_count(get_pvec(depth), get_bitmap(p_idx), root_vector_size);

//             // method 2 using avx2
//             // bitwise_andn(get_pvec(depth), get_bitmap(p_idx), get_nvec(depth), aligned_root_vector_size);
//             // int p_cnt = count_bitmap(get_nvec(depth), root_vector_size);
//             if (p_cnt >= max_intersect_cnt)
//             {
//                 max_intersect_cnt = p_cnt;
//                 pivot_index = p_idx;
//             }
//         }
//         pivot_inter_cnt[depth] = max_intersect_cnt;
//         num = maskz_expand_avx2_compress(get_bitmap(pivot_index), get_pvec(depth), simd_buffer, get_stack(depth), root_vector_size);
//     }
//     else
//     {
//         num = expand_avx2_compress(get_pvec(depth), next, root_vector_size);
//     }

//     // decode method 1, using avx2, minimum memory footprint
//     // decode method 2, using avx2
//     // bitwise_andn(get_pvec(depth), get_bitmap(pivot_index), get_nvec(depth), aligned_root_vector_size);
//     // int num = expand_avx2_compress(get_nvec(depth), get_stack(depth), root_vector_size);

//     // decode method 3, use ctz
//     // bitwise_andn(get_pvec(depth), get_bitmap(pivot_index), get_nvec(depth), aligned_root_vector_size);
//     // int num = expand_ctz(get_nvec(depth), get_stack(depth), root_vector_size);

//     for (int i = 0; i < num; i++)
//     {
//         dfs_pivot(next[i], depth + 1);
//         mark_as_zero(get_pvec(depth), next[i]);
//         mark_as_one(get_xvec(depth), next[i]);
//     }
// }

// void LoiMaximalClique::dfs_pivot_normal_onepunch()
// {
//     int depth = ++cur_depth;
//     int remain_stack_size = stack_set_size[depth - 1]--;
//     if (remain_stack_size <= 0)
//     {
//         --cur_depth;
//         return;
//     }
//     int v_index = get_stack(depth - 1)[remain_stack_size - 1];
//     R[depth] = pool_edges[root_start + v_index];
//     // compute the cliques
//     bitwise_and(get_bitmap(v_index), get_pvec(depth - 1), get_pvec(depth), aligned_root_vector_size);
//     bitwise_and(get_bitmap(v_index), get_xvec(depth - 1), get_xvec(depth), aligned_root_vector_size);
//     if (all_zero(get_pvec(depth), root_vector_size))
//     {
//         // declare maximal clique is found
//         if (all_zero(get_xvec(depth), root_vector_size))
//         {
//             mc_num++;
//             if (pool_mc_idx + depth + 1 < PACK_NODE_POOL_SIZE)
//             {
//                 // std::sort(R, R + depth);
//                 memcpy(pool_mc + pool_mc_idx, R, (depth + 1) * sizeof(int));
//                 pool_mc_idx += depth + 1;
//                 pool_mc[pool_mc_idx++] = -1;
//             }
//             maximum_clique_size = std::max(maximum_clique_size, depth + 1);
//         }

//         // max_pool_sets_idx = std::max(max_pool_sets_idx, depth + 1);
//         --cur_depth;
//         return;
//     }
//     int num;
//     int *next = get_stack(depth);
//     // choose a pivot point
//     // herustic 1, choose v with max degree, (assume graph is reindexed, small id high deg)
//     // int pivot_index = find_first_index(get_pvec(depth), root_vector_size);

//     // herustic 2, choose v with most triangles
//     // int pivot_index = -1, max_cnt = 0;
//     // int pnum = expand_avx2_compress(get_pvec(depth), get_stack(depth), root_vector_size);
//     // int *pnext = get_stack(depth);
//     // for (int i = 0; i < pnum; i++)
//     // {
//     //     int p_idx = pnext[i];
//     //     if (triangle_cnt[p_idx] > max_cnt)
//     //     {
//     //         pivot_index = p_idx;
//     //         max_cnt = triangle_cnt[p_idx];
//     //     }
//     // }

//     // herustic 3, choose v thath minimize P \ N(v)
//     if (pivot_inter_cnt[depth - 1] > 0)
//     {
//         int pivot_index = -1, max_intersect_cnt = 0;
//         int *pnext = get_stack(depth);
//         // not considering x
//         // int pnum = expand_avx2_compress(get_pvec(depth), get_stack(depth), root_vector_size);

//         // considering x
//         // bitwise_or(get_pvec(depth), get_xvec(depth), get_nvec(depth), aligned_root_vector_size);
//         // bitwise_andn(get_nvec(depth), get_xvec(0), get_nvec(0), aligned_root_vector_size);
//         // int pnum = expand_avx2_compress(get_pvec(depth), get_stack(depth), root_vector_size);
//         int pnum = maskzor_expand_avx2_compress(get_xvec(0), get_pvec(depth), get_xvec(depth), simd_buffer, pnext, root_vector_size);
//         for (int i = 0; i < pnum; i++)
//         {
//             int p_idx = pnext[i];
//             // method 1 using avx2, minimum memory footprint
//             int p_cnt = bitwise_and_count(get_pvec(depth), get_bitmap(p_idx), root_vector_size);

//             // method 2 using avx2
//             // bitwise_andn(get_pvec(depth), get_bitmap(p_idx), get_nvec(depth), aligned_root_vector_size);
//             // int p_cnt = count_bitmap(get_nvec(depth), root_vector_size);
//             if (p_cnt >= max_intersect_cnt)
//             {
//                 max_intersect_cnt = p_cnt;
//                 pivot_index = p_idx;
//             }
//         }
//         pivot_inter_cnt[depth] = max_intersect_cnt;
//         num = maskz_expand_avx2_compress(get_bitmap(pivot_index), get_pvec(depth), simd_buffer, get_stack(depth), root_vector_size);
//     }
//     else
//     {
//         num = expand_avx2_compress(get_pvec(depth), next, root_vector_size);
//     }

//     // decode method 1, using avx2, minimum memory footprint
//     // decode method 2, using avx2
//     // bitwise_andn(get_pvec(depth), get_bitmap(pivot_index), get_nvec(depth), aligned_root_vector_size);
//     // int num = expand_avx2_compress(get_nvec(depth), get_stack(depth), root_vector_size);

//     // decode method 3, use ctz
//     // bitwise_andn(get_pvec(depth), get_bitmap(pivot_index), get_nvec(depth), aligned_root_vector_size);
//     // int num = expand_ctz(get_nvec(depth), get_stack(depth), root_vector_size);
//     stack_set_size[depth] = num;
//     for (int i = num - 1; i >= 0; i--)
//     {
//         dfs_pivot_normal_onepunch();
//         mark_as_zero(get_pvec(depth), next[i]);
//         mark_as_one(get_xvec(depth), next[i]);
//     }
//     --cur_depth;
// }
void LoiMaximalClique::dfs_pivot_serious_onepunch()
{
    bool need_to_mask_last_index = false;
    cur_depth = 0;
    while (cur_depth >= 0)
    {
        int stack_idx = stack_set_size[cur_depth] - 1;
        if (stack_idx < 0)
        {
            --cur_depth;
            continue;
        }
        if (need_to_mask_last_index)
        {
            // current stack has been fully explored
            mark_as_zero(get_pvec(cur_depth), index_vec[cur_depth + 1]);
            mark_as_one(get_xvec(cur_depth), index_vec[cur_depth + 1]);
        }
        // visit a vertex in the stack
        int v_index = get_stack(cur_depth)[stack_idx];
        stack_set_size[cur_depth]--;

        // goes deeper
        ++cur_depth;
        index_vec[cur_depth] = v_index;
        R[cur_depth] = pool_edges[root_start + v_index];
        // compute the cliques
        // bitwise_and(get_bitmap(v_index), get_pvec(cur_depth - 1), get_pvec(cur_depth), aligned_root_vector_size);
        // bitwise_and(get_bitmap(v_index), get_xvec(cur_depth - 1), get_xvec(cur_depth), aligned_root_vector_size);
        double_bitwise_and(get_bitmap(v_index), get_pvec(cur_depth - 1), get_xvec(cur_depth - 1), get_pvec(cur_depth), get_xvec(cur_depth), aligned_root_vector_size);
        if (all_zero(get_pvec(cur_depth), root_vector_size))
        {
            // declare maximal clique is found
            if (all_zero(get_xvec(cur_depth), root_vector_size))
            {
                mc_num++;
                total_mc_size += cur_depth + 1;
                if (pool_mc_idx + cur_depth + 1 < PACK_NODE_POOL_SIZE)
                {
                    // std::sort(R, R + cur_depth);
                    memcpy(pool_mc + pool_mc_idx, R, (cur_depth + 1) * sizeof(int));
                    pool_mc_idx += cur_depth + 1;
                    pool_mc[pool_mc_idx++] = -1;
                }
                maximum_clique_size = std::max(maximum_clique_size, cur_depth + 1);
            }
            // done with the current vertex, mark as visited
            cur_depth--;
            need_to_mask_last_index = true;
            continue;
        }
        int num;
        int *next = get_stack(cur_depth);
        // choose a pivot point
        // herustic 1, choose v with max degree, (assume graph is reindexed, small id high deg)
        // int pivot_index = maskzor_find_first_index(get_xvec(0), get_pvec(cur_depth), get_xvec(cur_depth), root_vector_size);
        // assert(pivot_index != -1);
        // num = maskz_expand_avx2_compress(get_bitmap(pivot_index), get_pvec(cur_depth), simd_buffer, get_stack(cur_depth), root_vector_size);

        // herustic 2, choose v with most triangles
        // int pivot_index = -1, max_cnt = 0;
        // int pnum = expand_avx2_compress(get_pvec(depth), get_stack(depth), root_vector_size);
        // int *pnext = get_stack(depth);
        // for (int i = 0; i < pnum; i++)
        // {
        //     int p_idx = pnext[i];
        //     if (triangle_cnt[p_idx] > max_cnt)
        //     {
        //         pivot_index = p_idx;
        //         max_cnt = triangle_cnt[p_idx];
        //     }
        // }

        // herustic 3, choose v thath minimize P \ N(v)
        if (pivot_inter_cnt[cur_depth - 1] > 0)
        {
            int pivot_index = -1, max_intersect_cnt = 0;
            int *pnext = get_stack(cur_depth);
// not considering x
// int pnum = expand_avx2_compress(get_pvec(depth), get_stack(depth), root_vector_size);

// considering x
// bitwise_or(get_pvec(depth), get_xvec(depth), get_nvec(depth), aligned_root_vector_size);
// bitwise_andn(get_nvec(depth), get_xvec(0), get_nvec(0), aligned_root_vector_size);
// int pnum = expand_avx2_compress(get_pvec(depth), get_stack(depth), root_vector_size);
#if __SIMD_LEVEL__ == 2
            int pnum = maskzor_expand_avx2_compress(get_xvec(0), get_pvec(cur_depth), get_xvec(cur_depth), simd_buffer, pnext, root_vector_size);
#elif __SIMD_LEVEL__ == 0
            int pnum = maskzor_expand_ctz(get_xvec(0), get_pvec(cur_depth), get_xvec(cur_depth), pnext, root_vector_size);
#endif
            for (int i = 0; i < pnum; i++)
            {
                int p_idx = pnext[i];
                // method 1 using avx2, minimum memory footprint
                int p_cnt = bitwise_and_count(get_pvec(cur_depth), get_bitmap(p_idx), root_vector_size);

                // method 2 using avx2
                // bitwise_andn(get_pvec(depth), get_bitmap(p_idx), get_nvec(depth), aligned_root_vector_size);
                // int p_cnt = count_bitmap(get_nvec(depth), root_vector_size);
                if (p_cnt >= max_intersect_cnt)
                {
                    max_intersect_cnt = p_cnt;
                    pivot_index = p_idx;
                }
            }
            pivot_inter_cnt[cur_depth] = max_intersect_cnt;
#if __SIMD_LEVEL__ == 2
            num = maskz_expand_avx2_compress(get_bitmap(pivot_index), get_pvec(cur_depth), simd_buffer, pnext, root_vector_size);
#elif __SIMD_LEVEL__ == 0
            num = maskz_expand_ctz(get_bitmap(pivot_index), get_pvec(cur_depth), pnext, root_vector_size);
#endif
        }
        else
        {
#if __SIMD_LEVEL__ == 2

            num = expand_avx2_compress(get_pvec(cur_depth), next, root_vector_size);
#elif __SIMD_LEVEL__ == 0
            num = expand_ctz(get_pvec(cur_depth), next, root_vector_size);

#endif
        }
        // decode method 1, using avx2, minimum memory footprint
        // decode method 2, using avx2
        // bitwise_andn(get_pvec(depth), get_bitmap(pivot_index), get_nvec(depth), aligned_root_vector_size);
        // int num = expand_avx2_compress(get_nvec(depth), get_stack(depth), root_vector_size);

        // decode method 3, use ctz
        // bitwise_andn(get_pvec(depth), get_bitmap(pivot_index), get_nvec(depth), aligned_root_vector_size);
        // int num = expand_ctz(get_nvec(depth), get_stack(depth), root_vector_size);
        stack_set_size[cur_depth] = num;
        cur_depth++;
        need_to_mask_last_index = false;
    }
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
}

void LoiMaximalClique::report_mc_num()
{
    double counter = 0.0;
    long long last_mc_size = 0;
    long long last_mc_cnt = 0;
    long long last_compressed_mc_size = 0;
    std::cout << std::setw(15);
    std::cout << "executed |  mc number  | vertex num | total mc size | mc rate (K/s) | data (MB/s) | Xdata (MB/s) | current X ratio | Avg X ratio |" << std::endl;
    for (;;)
    {
        std::this_thread::sleep_for(std::chrono::seconds(REPORT_ELAPSE));
        counter += REPORT_ELAPSE;
        if (counter > MAX_REPORT_TIME)
        {
            break;
        }
        long long current_mc_size = total_mc_size;
        long long current_compressed_mc_size = total_compressed_mc_size;
        long long current_mc_cnt = mc_num;
        double data_rate = (current_mc_size - last_mc_size) * 4 / 1000000 / REPORT_ELAPSE;
        double compressed_data_rate = (current_compressed_mc_size - last_compressed_mc_size) / 1000000;
        double mc_rate = (current_mc_cnt - last_mc_cnt) / 1000;
        std::cout << std::setw(8) << counter << "s| ";
        std::cout << std::setw(11) << current_mc_cnt << " | " ;
        std::cout << std::setw(10) << u_cnt << " | " ;
        std::cout << std::setw(13) << current_mc_size << " | " ;
        std::cout << std::setw(13) << mc_rate << " | ";
        std::cout << std::setw(11) << data_rate << " | " ;
        std::cout << std::setw(12) << compressed_data_rate << " | " ;
        std::cout << std::setw(15) << data_rate / compressed_data_rate << " | ";
        std::cout << std::setw(11) << (double)(total_mc_size * 4) / (double)total_compressed_mc_size << " | "<< std::endl;
        last_mc_size = current_mc_size;
        last_compressed_mc_size = current_compressed_mc_size;
        last_mc_cnt = current_mc_cnt;
    }
}

std::string LoiMaximalClique::matrix_to_string()
{
    std::string result = "";
    for (int i = 0; i < root_deg; i++)
    {
        result += std::to_string(i) + " ";
        Bitmap *bitmap = &matrix[i * root_vector_size];
        result += bitmap_to_string(bitmap, root_vector_size);
        result += "\n";
    }
    result.pop_back();
    return result;
}
std::string LoiMaximalClique::to_string(int *vec, int size)
{
    std::string result = "";
    for (int i = 0; i < size; i++)
    {
        result += std::to_string(vec[i]) + " ";
    }
    return result;
}

std::string LoiMaximalClique::bitmap_to_string(Bitmap *bitmap_vec, int vector_size)
{
    std::string result = "";
    for (int j = 0; j < vector_size; j++)
    {
        Bitmap bitmap = bitmap_vec[j];
        for (int k = 0; k < 8; k++)
        {
            if (bitmap & (1llu << k))
            {
                result += "1";
            }
            else
            {
                result += "0";
            }
        }
        result += " ";
    }
    return result;
}