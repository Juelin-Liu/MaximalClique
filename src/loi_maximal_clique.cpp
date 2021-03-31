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

int LoiMaximalClique::build_matrix(const QVertex &a)
{
    int triangles = 0;
    root_start = a.start;
    root_offset = a.offset;
    root_deg = a.deg;
    root_vector_size = get_vector_size<Bitmap>(a.deg);
    memset(matrix, 0, a.deg * root_vector_size * sizeof(Bitmap));
    for (int b_idx = 0; b_idx < a.deg; b_idx++)
    {
        const int b_id = pool_edges[a.start + b_idx];
        const QVertex &b = graph[b_id];
        triangles += mark_intersect_simd8x(pool_edges + a.start, a.deg,
                                           pool_edges + b.start, b.deg, &matrix[b_idx * root_vector_size]);
    }
    return triangles;
}

int LoiMaximalClique::maximal_clique_bk()
{
    int buffer_vector_size = get_vector_size<Bitmap>(max_deg);
    matrix = (Bitmap *)calloc(max_deg * buffer_vector_size, sizeof(Bitmap));
    P_vec_pool = (Bitmap *)calloc((max_deg + 1) * buffer_vector_size, sizeof(Bitmap));
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
        memset(P_vec_pool, 0xff, root_vector_size * sizeof(Bitmap));
        for (int v_index = u.offset; v_index < u.deg; v_index++)
        {
            dfs_avx2(v_index, 1);
        }
    }

    free(matrix);
    free(P_vec_pool);
    free(id_vec);
    free(index_vec);
    free(index_pool);

    printf("max_pool_sets_idx=%d\n", max_pool_sets_idx);
    printf("maximum_clique_size=%d\n", maximum_clique_size);
    return mc_num;
}

// TODO implement this
int LoiMaximalClique::maximal_clique_pivot()
{
    int buffer_vector_size = get_vector_size<Bitmap>(max_deg);
    matrix = (Bitmap *)calloc(max_deg * buffer_vector_size, sizeof(Bitmap));
    P_vec_pool = (Bitmap *)calloc((max_deg + 1) * buffer_vector_size, sizeof(Bitmap));
    X_vec_pool = (Bitmap *)calloc((max_deg + 1) * buffer_vector_size, sizeof(Bitmap));
    index_pool = (int *)calloc(max_deg * buffer_vector_size, sizeof(int));
    // R_vec = (Bitmap *)calloc(buffer_vector_size, sizeof(Bitmap));
    // X_vec = (Bitmap *)calloc(buffer_vector_size, sizeof(Bitmap));
    id_vec = (int *)calloc(max_deg + 1, sizeof(int));
    index_vec = (int *)calloc(max_deg, sizeof(int));
    index_vec[0] = -1;
    mc_num = 0;
    u_id = 0;

    for (QVertex &u : graph)
    {
        id_vec[0] = u_id++;
        build_matrix(u);
        memset(get_xvec(0), 0, root_vector_size);
        memset(get_pvec(0), 0, root_vector_size);
        fill_with_one(get_pvec(0), root_deg);
        fill_with_one(get_xvec(0), root_offset);
        // choose the pivot as the first one
        int pivot_index = root_deg / 2; // (herustic)
        int next[root_deg + PADDING];
        Bitmap next_vec[root_vector_size];
        bitwise_andn(get_pvec(0), get_xvec(0), next_vec, root_vector_size);
        // // visit the candidates, ignore pivoting vertex's neighbor
        bitwise_andn(next_vec, get_bitmap(pivot_index), next_vec, root_vector_size);
        int num = expandToIndex(next_vec, next, root_vector_size);

        // LOG("id: " << id_vec[0] << " offset: " << root_offset << " deg: " << root_deg << " vector_size: " << root_vector_size);
        // LOG("matrix: \n"
        //     << matrix_to_string());
        // LOG("p_vec: " << bitmap_to_string(get_pvec(0), root_vector_size));
        // LOG("x_vec: " << bitmap_to_string(get_xvec(0), root_vector_size));
        // LOG("pool_edges: " << to_string(pool_edges + root_start, root_deg));
        // LOG("pivot index: " << pivot_index);
        // LOG("next: " << to_string(next, num));
        for (int i = 0; i < num; i++)
        {
            int v_index = next[i];
            // int v_index = i;
            dfs_avx2_pivot(v_index, 1);
            mark_as_one(get_xvec(0), v_index);
            mark_as_zero(get_pvec(0), v_index);
        }
    }

    free(matrix);
    free(P_vec_pool);
    free(id_vec);
    free(index_vec);
    free(index_pool);

    printf("max_pool_sets_idx=%d\n", max_pool_sets_idx);
    printf("maximum_clique_size=%d\n", maximum_clique_size);
    return mc_num;
}

// TODO implement this
int LoiMaximalClique::maximal_clique_degen()
{
    return mc_num;
}

void LoiMaximalClique::dfs_avx2(int v_index, int depth)
{
    // check if v is part of a maximal clique
    index_vec[depth] = v_index;
    id_vec[depth] = pool_edges[root_start + v_index];
    intersect(get_bitmap(v_index), get_pvec(depth - 1), get_pvec(depth), root_vector_size);
    if (all_zero(get_pvec(depth), root_vector_size))
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
        if (root_start == 0)
        {
            LOG("id: " << id_vec[0] << " offset: " << root_offset << " deg: " << root_deg << " vector_size: " << root_vector_size);
            LOG("pool_edges: " << to_string(pool_edges + root_start, root_deg));
            LOG("matrix: \n"
                << matrix_to_string());
            LOG("P_vec_pool: " << bitmap_to_string(get_pvec(depth), root_vector_size));
            LOG("index_vec: " << to_string(index_vec, depth + 1));
            LOG("id_vec: " << to_string(id_vec, depth + 1) << "\n");
        }
        return;
    }
    int num = expandToIndex(get_pvec(depth), get_index_by_depth(depth), v_index + 1, root_deg);
    int *next = get_index_by_depth(depth);
    for (int i = 0; i < num; i++)
    {
        dfs_avx2(next[i], depth + 1);
    }
}

void LoiMaximalClique::dfs_avx2_pivot(int v_index, int depth)
{

    // bookkeeping
    index_vec[depth] = v_index;
    id_vec[depth] = pool_edges[root_start + v_index];
    // compute the cliques formed with all visiting vertexes
    bitwise_and(get_bitmap(v_index), get_pvec(depth - 1), get_pvec(depth), root_vector_size);
    bitwise_and(get_bitmap(v_index), get_xvec(depth - 1), get_xvec(depth), root_vector_size);
    // LOG("depth: " << depth);
    // LOG("visiting indexes: " << to_string(index_vec, depth));
    // LOG("p_vec: " << bitmap_to_string(get_pvec(depth), root_vector_size));
    // LOG("x_vec: " << bitmap_to_string(get_xvec(depth), root_vector_size));

    // no more to explore
    // get candidates for next visit
    if (all_zero(get_pvec(depth), root_vector_size))
    {
        // declare maximal clique is found
        if (all_zero(get_xvec(depth), root_vector_size))
        {
            mc_num++;
            if (pool_mc_idx + depth + 1 < PACK_NODE_POOL_SIZE)
            {
                // std::sort(id_vec, id_vec + depth);
                memcpy(pool_mc + pool_mc_idx, id_vec, (depth + 1) * sizeof(int));
                pool_mc_idx += depth + 1;
                pool_mc[pool_mc_idx++] = -1;
            }
            maximum_clique_size = std::max(maximum_clique_size, depth + 1);
            // if (root_start == 0)
            // {
            //     LOG("id: " << id_vec[0] << " offset: " << root_offset << " deg: " << root_deg << " vector_size: " << root_vector_size);
            //     LOG("pool_edges: " << to_string(pool_edges + root_start, root_deg));
            //     LOG("matrix: \n"
            //         << matrix_to_string());
            //     LOG("P_vec_pool: " << bitmap_to_string(get_pvec(depth), root_vector_size));
            //     LOG("index_vec: " << to_string(index_vec, depth + 1));
            //     LOG("id_vec: " << to_string(id_vec, depth + 1) << "\n");
            // }
            // LOG("Found\n");
            // LOG("index_vec: " << to_string(index_vec, depth + 1));
        }

        // max_pool_sets_idx = std::max(max_pool_sets_idx, depth + 1);
        return;
    }
    // choose a pivot point
    Bitmap next_vec[root_vector_size];
    int pivot_index = find_first_index(get_pvec(depth), root_vector_size);
    assert(pivot_index != -1);
    // LOG("pivot: " << pivot_index);
    int next[root_deg + PADDING];
    bitwise_andn(get_pvec(depth), get_xvec(depth), next_vec, root_vector_size );
    bitwise_andn(next_vec, get_bitmap(pivot_index), next_vec, root_vector_size);
    int num = expandToIndex(next_vec, next, root_vector_size);
    // LOG("next: " << to_string(next, num) << "\n");
    for (int i = 0; i < num; i++)
    {
        dfs_avx2_pivot(next[i], depth + 1);
        mark_as_zero(get_pvec(depth), next[i]);
        mark_as_one(get_xvec(depth), next[i]);
    }
}

void LoiMaximalClique::dfs_clz(int v_index, int depth)
{
    // check if v is part of a maximal clique
    index_vec[depth] = v_index;
    id_vec[depth] = pool_edges[root_start + v_index];
    bool found = intersect_allzero(get_bitmap(v_index), get_pvec(depth - 1), get_pvec(depth), root_vector_size);
    if (found)
    {
        mc_num++;
        if (pool_mc_idx + depth + 1 < PACK_NODE_POOL_SIZE)
        {
            memcpy(pool_mc + pool_mc_idx, id_vec, (depth + 1) * sizeof(int));
            pool_mc_idx += depth + 1;
            pool_mc[pool_mc_idx++] = -1;
        }
        if (root_start == 0)
        {
            LOG("id: " << id_vec[0] << " offset: " << root_offset << " deg: " << root_deg << " vector_size: " << root_vector_size);
            LOG("pool_edges: " << to_string(pool_edges + root_start, root_deg));
            LOG("matrix: \n"
                << matrix_to_string());
            LOG("P_vec_pool: \n"
                << bitmap_to_string(get_pvec(depth), root_vector_size));
            LOG("index_vec: " << to_string(index_vec, depth + 1));
            LOG("id_vec: " << to_string(id_vec, depth + 1) << "\n");
        }
        return;
    }
    if (v_index + 1 == root_deg)
    {
        return;
    }
    Bitmap *bitvec = get_pvec(depth);
    int v_filter = v_index + 1;
    int v_start = v_filter / (sizeof(Bitmap) * 8);
    Bitmap bitset = bitvec[v_start];
    bitset = bitset & (0xffffffffffffffff >> (v_filter % 64));
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
            dfs_clz(new_index, depth + 1);
        }
    }
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