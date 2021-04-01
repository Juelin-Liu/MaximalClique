# CliqueX
Faster Maximal Clique Algorithm

### Project Structure
- `src` directory contains the source code. 
    - `org_maximal_clique` contains the naive implementation.
    - `bp_maximal_clique` contains the current state-of-the-art implementation.
    - `loi_maximal_clique` contains the this `CliqueX` implementation.`

Special thanks to Shuo Han, Lei Zou, and Jeffrey Xu Yu, who provided the [source code](https://github.com/pkumod/GraphSetIntersection) the current state-of-the-art implementation and many useful APIs.

- `data` contains the scripts for downloading data and reordering data.
- `paper` contains the research paper that builds the foundation of this project.

### How to Run
In `src` directory run `make`. It will generate two binaries: `mc`, which is for maximal clique finding, and `reorder` which is for reordering the graph. 

To find maximal clique in a graph use:
> ./mc [ path to graph file ] [ path to output file]

For instance find maximal clique in `data/reactome`:
> ./mc ../data/reactome

### Graphs Included
The graph included is `reactome`: 
| Graph Name | Vertex Number | Edge Number | Triangle Count |
| --- | --- | --- | --- |
| reactome | 6,327 | 147,547 | 4,187,734 |

### Introduction

This project aims to find faster algorithm to find maximal cliques in an undirected graph. Listing all maximal cliques is a NP hard problem. [Moon & Moser](https://link.springer.com/article/10.1007%2FBF02760024) have shown that every n-vertex graph has at most 3^(n/3) maximal cliques. 

They can be listed by the [Bron-Kerbosch algorithm](https://en.wikipedia.org/wiki/Bron–Kerbosch_algorithm). A recursive backtracking algorithm that is worst case optimal. 

Existing algorithms focus on speeding up maximal clique finding when its number is significant less than the worst case. Examples include [pivoting](https://en.wikipedia.org/wiki/Bron–Kerbosch_algorithm#With_pivoting), and [pivoting + ordering(degeneracy)](https://en.wikipedia.org/wiki/Bron–Kerbosch_algorithm#With_vertex_ordering) or [Tomita](https://snap.stanford.edu/class/cs224w-readings/tomita06cliques.pdf).

This project proposed a new DFS backtracking algorithm build on top of [Local Online Indexing (Loi)](./paper/LocalOnlineIndexing.pdf). `Loi` is a compact encoding of a vertex's triangles in binary format. It uses SIMD instructions to accelerate bitmap encoding and decoding. It uses bitwise AND operation to get set intersection instead of traditional vertex at a time approach. Simply put, it is designed to be cache friendly, SIMD friendly, and branchless.

### Bron-Kerbosch Recap - From WikiPedia
The core part of Bron-Kerbosch algorithm consists of three arrays: R, P, and X.

R contains all vertexes the current DFS tree path.
P contains all vertexes that is connected to all vertex in R.
X contains all vertexes that has been fully explored.

Initially, R is empty, P contains all vertexes in the graph, and X is empty.

The algorithm goes like this:
```
algorithm BronKerbosch1(R, P, X) is
    if P and X are both empty then
        report R as a maximal clique
    for each vertex v in P do
        BronKerbosch1(R ⋃ {v}, P ⋂ N(v), X ⋂ N(v))
        P := P \ {v}
        X := X ⋃ {v}
```
This algorithm will conduct a DFS on each vertex in the graph. It recursively finds vertex that can is connected to all vertexes in R.

When P is empty, the current DFS path has no more vertex to explore so it can return. 
But before that, it needs to check if R contains a maximal clique. The way to do that is check if X is empty. If X is not empty, this maximal clique must has been explored in the previous search. (Remember or X is builded, only the vertexes that have been fully explored are in X. Also they are connected to all vertex in R.)

### Pivoting Vertex

The pitfall of this algorithm is it needs to open a recursive call on every cliques. But many are not maximal cliques. A remedy is to choose a pivoting vertex.
> To save time and allow the algorithm to backtrack more quickly in branches of the search that contain no maximal cliques, Bron and Kerbosch introduced a variant of the algorithm involving a "pivot vertex" u, chosen from P (or more generally, as later investigators realized, from P ⋃ X). Any maximal clique must include either u or one of its non-neighbors, for otherwise the clique could be augmented by adding u to it. Therefore, only u and its non-neighbors need to be tested as the choices for the vertex v that is added to R in each recursive call to the algorithm. In pseudocode:

```
algorithm BronKerbosch2(R, P, X) is
    if P and X are both empty then
        report R as a maximal clique
    choose a pivot vertex u in P ⋃ X
    for each vertex v in P \ N(u) do
        BronKerbosch2(R ⋃ {v}, P ⋂ N(v), X ⋂ N(v))
        P := P \ {v}
        X := X ⋃ {v}
```

If the pivot is chosen to minimize the number of recursive calls made by the algorithm, the savings in running time compared to the non-pivoting version of the algorithm can be significant.

### Pivot Vertex - Choose The First One Wisely 

An alternative method for improving the basic form of the Bron–Kerbosch algorithm involves forgoing pivoting at the outermost level of recursion, and instead choosing the ordering of the recursive calls carefully in order to minimize the sizes of the sets P of candidate vertices within each recursive call.

The degeneracy of a graph G is the smallest number d such that every subgraph of G has a vertex with degree d or less. Every graph has a degeneracy ordering, an ordering of the vertices such that each vertex has d or fewer neighbors that come later in the ordering; a degeneracy ordering may be found in linear time by repeatedly selecting the vertex of minimum degree among the remaining vertices. If the order of the vertices v that the Bron–Kerbosch algorithm loops through is a degeneracy ordering, then the set P of candidate vertices in each call (the neighbors of v that are later in the ordering) will be guaranteed to have size at most d. The set X of excluded vertices will consist of all earlier neighbors of v, and may be much larger than d. In recursive calls to the algorithm below the topmost level of the recursion, the pivoting version can still be used.[6][7]

In pseudocode, the algorithm performs the following steps:

```
algorithm BronKerbosch3(G) is
    P = V(G)
    R = X = empty
    for each vertex v in a degeneracy ordering of G do
        BronKerbosch2({v}, P ⋂ N(v), X ⋂ N(v))
        P := P \ {v}
        X := X ⋃ {v}
```

### Pivot Vertex - Choose All Wisely
Glad you have read to this point because. Here is what this project really about.

The degeneracy ordering mentioned above only handles the first pivoting vertex. But what about the rest?

An intuitive way to think about pivoting vertex is it makes the DFS tree thinner by not visiting many possible branches. P becomes more and more selectively as the DFS depth increases, and probabliy we don't need to visit those branches at all. 

A greedy way of selecting the pivoting vertex tries to minimize the new DFS calls in the current depth. To do this, for all vertex in P we need to compute the intersection of this vertex's neighbors with P. Then we choose the vertex with minimum intersection size as the pivoting vertex.

In pseudocode:
```
algorithm BronKerboschGreedy(R, P, X) is
    if P and X are both empty then
        report R as a maximal clique
    
    for each vertex u in P do:
        U := P ⋂ N(u)
    choose u with Min(|U|)

    for each vertex v in P \ N(u) do
        BronKerboschGreedy(R ⋃ {v}, P ⋂ N(v), X ⋂ N(v))
        P := P \ {v}
        X := X ⋃ {v}
```
