# CliqueX
A new algorithm for clique finding in dense graphs.
### Project Structure
- `src` directory contains the source code. Many thanks to Shuo Han, Lei Zou, and Jeffrey Xu Yu, who provided the [baseline](https://github.com/pkumod/GraphSetIntersection) for comparison. 
- `data` contains the scripts for downloading data.
- `paper` contains the research paper that builds the foundation of this project.

### How to Run
In `src` directory run `make`. It will generate two binaries: `mc`, which is for maximal clique finding, and `reorder` which is for reordering the graph. 

To find maximal clique in a graph use:
> ./mc [ path to graph file ] [ path to output file]

For downloading the graph file. Switch to `data` directory, run:
> bash init.sh

The `init.sh` script will download four undirected datasets from [SNAP](http://snap.stanford.edu/data/) and preprocess the downloaded graphs. 

After the graph is downloaded, you can the following command to find all maximal clique in graph `amazon` and output to `amaozn_mc`:
> ./src/mc ./data/amazon ./mc_result/amazon_mc

### Graphs Downloaded
The four graphs are: 
| Graph Name | Vertex Number | Edge Number |
| --- | --- | --- |
| amazon | 334,863 | 925,872 |
| youtube | 1,134,890 | 2,987,624 |
| lj | 3,997,962 | 34,681,189 |
| orkut | 3,072,441 | 117,185,083 |


### Introduction
According to WikiPedia:
> In computer science, the clique problem is the computational problem of finding cliques (subsets of vertices, all adjacent to each other, also called complete subgraphs) in a graph.

This project aims to find faster algorithm to find maximal cliques in an undirected graph. Listing all maximal cliques is a NP hard problem. [Moon & Moser](https://link.springer.com/article/10.1007%2FBF02760024) have shown that every n-vertex graph has at most 3^(n/3) maximal cliques. They can be listed by the [Bron-Kerbosch algorithm](https://en.wikipedia.org/wiki/Bron–Kerbosch_algorithm). A recursive backtracking algorithm that is worst case optimal. 

Existing algorithms focus on speeding up maximal clique finding when its number is significant less than the worst case. Examples include [pivoting](https://en.wikipedia.org/wiki/Bron–Kerbosch_algorithm#With_pivoting), [vertex ordering(degeneracy)](https://en.wikipedia.org/wiki/Bron–Kerbosch_algorithm#With_vertex_ordering), and [Tomita](https://snap.stanford.edu/class/cs224w-readings/tomita06cliques.pdf).

This project proposed a new DFS backtracking algorithm build on top of [Local Online Indexing (Loi)](./paper/LocalOnlineIndexing.pdf). `Loi` is a compact encoding of a vertex's triangles in binary format. It uses SIMD instructions to accelerate bitmap encoding and decoding. It uses bitwise AND operation to get set intersection results rather than traditional vertex at a time approach. Simply put, it is designed to be cache friendly, SIMD friendly, and branchless.

