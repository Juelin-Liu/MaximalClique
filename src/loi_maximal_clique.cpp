#include "loi_maximal_clique.hpp"
#define LOG(x) std::cout << x << std::endl
#define PADDING 16
LoiMaximalClique::LoiMaximalClique()
{
    v_num = 0;
    e_num = 0;
    align_malloc((void **)&simd_buffer, sizeof(AlignType), sizeof(AlignType));
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
    align_malloc((void **)&next_vec_pool, sizeof(AlignType), (deg + 1) * aligned_buffer_vector_size);
    R = (int *)malloc(sizeof(int) * deg);
    X = (int *)malloc(sizeof(int) * deg);
    next_set_pool = (int *)malloc(sizeof(int) * deg * deg);
    index_vec = (int *)malloc(sizeof(int) * (deg + 1));
    triangle_cnt = (int *)malloc(sizeof(int) * (deg + 1));
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
    if (next_set_pool != NULL)
    {
        free(next_set_pool);
    }
    if (index_vec != NULL)
    {
        free(index_vec);
    }
    if (R != NULL)
    {
        free(R);
    }
    if(X != NULL){
        free(X);
    }
    if(triangle_cnt != NULL){
        free(triangle_cnt);
    }
}
void LoiMaximalClique::init_R_X(const QVertex &a){
    aligned_root_vector_size = get_vector_size<AlignType>(a.deg);
    memset(get_xvec(0), 0, aligned_root_vector_size* sizeof(Bitmap));
    memset(get_pvec(0), 0, aligned_root_vector_size* sizeof(Bitmap));
    for (int b_idx = 0; b_idx < a.deg; b_idx++)
    {
        const int b_id = pool_edges[a.start + b_idx];
        if (visited[b_id])
        {
            mark_as_one(get_xvec(0), b_idx);
        } else {
            mark_as_one(get_pvec(0), b_idx);
        }
    }

}
int LoiMaximalClique::build_matrix(const QVertex &a)
{
    int triangles = 0;
    int p_index = 0, p_triangle_size = 0;
    root_start = a.start;
    root_offset = a.offset;
    root_deg = a.deg;
    root_vector_size = get_vector_size<Bitmap>(a.deg);
    aligned_root_vector_size = get_vector_size<AlignType>(a.deg);
    memset(get_bitmap(0), 0, a.deg * aligned_root_vector_size * sizeof(Bitmap));
    memset(triangle_cnt, 0, a.deg * sizeof(int));
    for (int b_idx = 0; b_idx < a.deg; b_idx++)
    {
        const int b_id = pool_edges[a.start + b_idx];
        if (!visited[b_id])
        {
            const QVertex &b = graph[b_id];
            int t_num = mark_intersect_simd8x(pool_edges + a.start, a.deg,
                                               pool_edges + b.start, b.deg, get_bitmap(b_idx)); 
            triangle_cnt[b_idx] = t_num;
            if (t_num > p_triangle_size) {
                p_triangle_size = t_num;
                p_index = b_idx;
            }      
            // bitwise_andn(get_bitmap(b_idx), get_xvec(0), get_bitmap(b_idx), aligned_root_vector_size);
        }
    }
    return p_index;
}

int LoiMaximalClique::maximal_clique_bk2()
{
    int buffer_vector_size = get_vector_size<Bitmap>(max_deg);
    matrix = (Bitmap *)calloc(max_deg * buffer_vector_size, sizeof(Bitmap));
    P_vec_pool = (Bitmap *)calloc((max_deg + 1) * buffer_vector_size, sizeof(Bitmap));
    index_pool = (int *)calloc(max_deg * buffer_vector_size, sizeof(int));
    R = (int *)calloc(max_deg + 1, sizeof(int));
    index_vec = (int *)calloc(max_deg, sizeof(int));
    index_vec[0] = -1;
    mc_num = 0;
    u_cnt = 0;
    for (QVertex &u : graph)
    {
        R[0] = u_cnt++;
        build_matrix(u);
        memset(P_vec_pool, 0xff, root_vector_size * sizeof(Bitmap));
        for (int v_index = u.offset; v_index < u.deg; v_index++)
        {
            dfs(v_index, 1);
        }
    }

    free(matrix);
    free(P_vec_pool);
    free(R);
    free(index_vec);
    free(index_pool);

    printf("max_pool_sets_idx=%d\n", max_pool_sets_idx);
    printf("maximum_clique_size=%d\n", maximum_clique_size);
    return mc_num;
}

int LoiMaximalClique::maximal_clique_bk()
{
    set_buffer_capacity(max_deg);
    index_vec[0] = -1;
    mc_num = 0, u_cnt = 0, p_set_idx = 0;
    visited = new bool[v_num];
    for (QVertex &u : graph)
    {
        R[0] = u_cnt;
        build_matrix(u);
        // fill_with_one(get_xvec(0), root_offset);
        bitwise_andn(get_pvec(0), get_xvec(0), get_nvec(0), aligned_root_vector_size);
        // choose the pivot as the first one
        // int pivot_index = root_deg / 2; // (herustic)
        // int next[root_deg + PADDING];
        // // // visit the candidates, ignore pivoting vertex's neighbor
        // bitwise_andn(next_vec, get_bitmap(pivot_index), next_vec, aligned_root_vector_size);
        int num = expand_avx2(get_nvec(0), get_nset(0), root_vector_size);
        int *next = get_nset(0);
        // LOG("id: " << R[0] << " offset: " << root_offset << " deg: " << root_deg << " vector_size: " << root_vector_size);
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
            dfs(v_index, 1);
            mark_as_one(get_xvec(0), v_index);
            mark_as_zero(get_pvec(0), v_index);
        }
        visited[u_cnt++] = true;
    }
    delete[] visited;
    free_buffer();
    printf("max_pool_sets_idx=%d\n", max_pool_sets_idx);
    printf("maximum_clique_size=%d\n", maximum_clique_size);
    return mc_num;
}

// TODO implement this
int LoiMaximalClique::maximal_clique_pivot()
{
    set_buffer_capacity(max_deg);
    index_vec[0] = -1;
    mc_num = 0, u_cnt = 0, p_set_idx = 0;
    visited = new bool[v_num];
    for (QVertex &u : graph)
    {
        R[0] = u_cnt;
        
        if (visited[u_cnt]){
            u_cnt++;
            continue;
        }
        init_R_X(u);
        // choose the first vertex from X otherwise choose from P
        int pivot_index = build_matrix(u);
        bitwise_andn(get_pvec(0), get_bitmap(pivot_index), get_nvec(0), aligned_root_vector_size);
        int num = expand_ctz(get_nvec(0), get_nset(0), root_vector_size);
        int *next = get_nset(0);
        // LOG("id: " << R[0] << " offset: " << root_offset << " deg: " << root_deg << " vector_size: " << root_vector_size);
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
            dfs_pivot(v_index, 1);
            mark_as_one(get_xvec(0), v_index);
            mark_as_zero(get_pvec(0), v_index);
        }
        visited[u_cnt++] = true;
    }
    delete[] visited;
    free_buffer();
    printf("max_pool_sets_idx=%d\n", max_pool_sets_idx);
    printf("maximum_clique_size=%d\n", maximum_clique_size);
    return mc_num;
}

// TODO implement this
int LoiMaximalClique::maximal_clique_degen()
{
    set_buffer_capacity(max_deg);
    index_vec[0] = -1;
    mc_num = 0, u_cnt = 0, p_set_idx = 0;
    visited = new bool[v_num];
    std::vector<int> dorder = degeneracy_order();
    for (int v : dorder)
    {
        u_cnt++;
        R[0] = v;
        const QVertex &u = graph[v];
        init_R_X(u);
        // choose the first vertex from X otherwise choose from P
        int pivot_index = build_matrix(u);
        bitwise_andn(get_pvec(0), get_bitmap(pivot_index), get_nvec(0), aligned_root_vector_size);
        int num = expand_ctz(get_nvec(0), get_nset(0), root_vector_size);
        int *next = get_nset(0);
        // LOG("id: " << R[0] << " offset: " << root_offset << " deg: " << root_deg << " vector_size: " << root_vector_size);
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
            dfs_pivot(v_index, 1);
            mark_as_one(get_xvec(0), v_index);
            mark_as_zero(get_pvec(0), v_index);
        }
        visited[v] = true;
    }
    delete[] visited;
    free_buffer();
    printf("max_pool_sets_idx=%d\n", max_pool_sets_idx);
    printf("maximum_clique_size=%d\n", maximum_clique_size);
    return mc_num;
}

void LoiMaximalClique::dfs(int v_index, int depth)
{
    // bookkeeping
    index_vec[depth] = v_index;
    R[depth] = pool_edges[root_start + v_index];
    // compute the cliques formed with all visiting vertexes
    bitwise_and(get_bitmap(v_index), get_pvec(depth - 1), get_pvec(depth), aligned_root_vector_size);
    bitwise_and(get_bitmap(v_index), get_xvec(depth - 1), get_xvec(depth), aligned_root_vector_size);
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
                // std::sort(R, R + depth);
                memcpy(pool_mc + pool_mc_idx, R, (depth + 1) * sizeof(int));
                pool_mc_idx += depth + 1;
                pool_mc[pool_mc_idx++] = -1;
            }
            maximum_clique_size = std::max(maximum_clique_size, depth + 1);
        }
        return;
    }
    bitwise_andn(get_pvec(depth), get_xvec(depth), get_nvec(depth), aligned_root_vector_size);
    int num = expand_avx2(get_nvec(depth), get_nset(depth), root_vector_size);
    int *next = get_nset(depth);
    for (int i = 0; i < num; i++)
    {
        dfs(next[i], depth + 1);
        mark_as_zero(get_pvec(depth), next[i]);
        mark_as_one(get_xvec(depth), next[i]);
    }
}

void LoiMaximalClique::dfs_pivot(int v_index, int depth)
{

    // bookkeeping
    index_vec[depth] = v_index;
    R[depth] = pool_edges[root_start + v_index];
    // compute the cliques formed with all visiting vertexes
    bitwise_and(get_bitmap(v_index), get_pvec(depth - 1), get_pvec(depth), aligned_root_vector_size);
    bitwise_and(get_bitmap(v_index), get_xvec(depth - 1), get_xvec(depth), aligned_root_vector_size);
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
                // std::sort(R, R + depth);
                memcpy(pool_mc + pool_mc_idx, R, (depth + 1) * sizeof(int));
                pool_mc_idx += depth + 1;
                pool_mc[pool_mc_idx++] = -1;
            }
            maximum_clique_size = std::max(maximum_clique_size, depth + 1);
        }

        // max_pool_sets_idx = std::max(max_pool_sets_idx, depth + 1);
        return;
    }
    // choose a pivot point
    // herustic 1, choose v with max degree, (small id high deg)
    // int pivot_index = find_first_index(get_pvec(depth), root_vector_size);

    // herustic 2, choose v with most triangles 
    int pivot_index = -1, max_cnt = 0;
    int pnum = expand_avx2_compress(get_pvec(depth), get_nset(depth), root_vector_size);
    int* pnext = get_nset(depth);
    for (int i = 0; i < pnum; i++)
    {
        int p_idx = pnext[i];
        if (triangle_cnt[p_idx] >= max_cnt){
            pivot_index = p_idx;
            max_cnt = triangle_cnt[p_idx];
        }
    }

    // // herustic 3, choose v minimize P \ N(v) 
    // int pivot_index = -1, max_cnt = 0;
    // int pnum = expand_avx2_compress(get_pvec(depth), get_nset(depth), root_vector_size);
    // int* pnext = get_nset(depth);
    // for (int i = 0; i < pnum; i++)
    // {
    //     int p_idx = pnext[i];
    //     int p_cnt = bitwise_andn_count(get_pvec(depth), get_bitmap(p_idx), simd_buffer, root_vector_size);
    //     if(p_cnt >= max_cnt){
    //         max_cnt = p_cnt;
    //         pivot_index = p_idx;
    //     }
    // }
    assert(pivot_index != -1);
    bitwise_andn(get_pvec(depth), get_bitmap(pivot_index), get_nvec(depth), aligned_root_vector_size);
    int num = expand_avx2_compress(get_nvec(depth), get_nset(depth), root_vector_size);
    // int num = expand_ctz(get_nvec(depth), get_nset(depth), root_vector_size);
    int *next = get_nset(depth);
    for (int i = 0; i < num; i++)
    {
        dfs_pivot(next[i], depth + 1);
        mark_as_zero(get_pvec(depth), next[i]);
        mark_as_one(get_xvec(depth), next[i]);
    }
}

// void LoiMaximalClique::dfs_clz(int v_index, int depth)
// {
//     // check if v is part of a maximal clique
//     index_vec[depth] = v_index;
//     R[depth] = pool_edges[root_start + v_index];
//     bool found = intersect_allzero(get_bitmap(v_index), get_pvec(depth - 1), get_pvec(depth), root_vector_size);
//     if (found)
//     {
//         mc_num++;
//         if (pool_mc_idx + depth + 1 < PACK_NODE_POOL_SIZE)
//         {
//             memcpy(pool_mc + pool_mc_idx, R, (depth + 1) * sizeof(int));
//             pool_mc_idx += depth + 1;
//             pool_mc[pool_mc_idx++] = -1;
//         }
//         if (root_start == 0)
//         {
//             LOG("id: " << R[0] << " offset: " << root_offset << " deg: " << root_deg << " vector_size: " << root_vector_size);
//             LOG("pool_edges: " << to_string(pool_edges + root_start, root_deg));
//             LOG("matrix: \n"
//                 << matrix_to_string());
//             LOG("P_vec_pool: \n"
//                 << bitmap_to_string(get_pvec(depth), root_vector_size));
//             LOG("index_vec: " << to_string(index_vec, depth + 1));
//             LOG("R: " << to_string(R, depth + 1) << "\n");
//         }
//         return;
//     }
//     if (v_index + 1 == root_deg)
//     {
//         return;
//     }
//     Bitmap *bitvec = get_pvec(depth);
//     int v_filter = v_index + 1;
//     int v_start = v_filter / (sizeof(Bitmap) * 8);
//     Bitmap bitset = bitvec[v_start];
//     bitset = bitset & (0xffffffffffffffff >> (v_filter % 64));
//     while (bitset)
//     {
//         int r = __builtin_clzll(bitset);
//         int new_index = r + 64 * v_start;
//         bitset ^= 1LLU << (63 - r);
//         dfs_clz(new_index, depth + 1);
//     }
//     // COUNT LEADING ZERO
//     for (int i = v_start + 1; i < root_vector_size; i++)
//     {
//         bitset = bitvec[i];
//         while (bitset)
//         {
//             int r = __builtin_clzll(bitset);
//             int new_index = r + 64 * i;
//             bitset ^= 1LLU << (63 - r);
//             dfs_clz(new_index, depth + 1);
//         }
//     }
// }
std::vector<int> LoiMaximalClique::degeneracy_order()
{
    std::vector<int> deg(v_num);
    int md = 0;
    for (int i = 0; i < v_num; ++i) {
        deg[i] = graph[i].deg;
        md = std::max(md, deg[i]);
    }

    std::vector<int> bin(md + 1);
    for (int i = 0; i <= md; ++i) bin[i] = 0;
    for (int i = 0; i < v_num; ++i) bin[deg[i]]++;

    int start = 0;
    for (int i = 0; i <= md; ++i) {
        int num = bin[i];
        bin[i] = start;
        start += num;
    }

    std::vector<int> vert(v_num), pos(v_num);
    for (int i = 0; i < v_num; ++i) {
        pos[i] = bin[deg[i]];
        vert[pos[i]] = i;
        bin[deg[i]]++;
    }
    for (int i = md; i > 0; --i) bin[i] = bin[i - 1];
    bin[0] = 0;

    std::vector<int> degen_order;
    degen_order.reserve(v_num);
    int degeneracy = 0;
    for (int i = 0; i < v_num; ++i) {
        int v = vert[i];
        degen_order.push_back(v);
        degeneracy = std::max(degeneracy, deg[v]);
        for (int j = 0; j < graph[v].deg; ++j) {
            int u = pool_edges[graph[v].start + j];
            if (deg[u] > deg[v]) {
                int du = deg[u], pu = pos[u];
                int pw = bin[du], w = vert[pw];
                if (u != w) {
                    pos[u] = pw; vert[pu] = w;
                    pos[w] = pu; vert[pw] = u;
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
    int counter = 0;
    for (;;)
    {
        std::this_thread::sleep_for(std::chrono::seconds(REPORT_ELAPSE));
        counter++;
        if (counter > MAX_REPORT_TIME)
        {
            break;
        }
        std::cout << "executed: "<< counter << " s\t\tmc number: " << mc_num << "\t\tvertex num: " << u_cnt << std::endl;
    }
}