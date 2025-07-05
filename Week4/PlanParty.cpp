#include <iostream>
#include <vector>

struct Vertex {
    int weight;
    std::vector<int> adjacent_vertices;
};

using Tree = std::vector<Vertex>;
using MemoizationTable = std::vector<int>;

Tree ReadTree() noexcept
{
    std::ios::sync_with_stdio(false);

    int vertex_count;
    std::cin >> vertex_count;

    Tree tree(vertex_count);

    // Read weights for each vertex
    for (int i = 0; i < vertex_count; ++i) {
        std::cin >> tree[i].weight;
    }

    // Read edges and build adjacency lists
    for (int i = 1; i < vertex_count; ++i) {
        int from, to;
        std::cin >> from >> to;
        tree[from - 1].adjacent_vertices.push_back(to - 1);
        tree[to - 1].adjacent_vertices.push_back(from - 1);
    }

    return tree;
}

int ComputeMaxIndependentSetRecursive(const Tree& tree, int current_vertex, 
                                     MemoizationTable& memo, int parent_vertex) noexcept
{
    if (memo[current_vertex] == -1) {
        if (tree[current_vertex].adjacent_vertices.empty()) {
            memo[current_vertex] = tree[current_vertex].weight;
        }
        else {
            int include_current = tree[current_vertex].weight;
            for (int neighbor : tree[current_vertex].adjacent_vertices) {
                if (neighbor == parent_vertex)
                    continue;
                for (int grandchild : tree[neighbor].adjacent_vertices) {
                    if (grandchild != current_vertex) {
                        include_current += ComputeMaxIndependentSetRecursive(tree, grandchild, memo, neighbor);
                    }
                }
            }

            int exclude_current = 0;
            for (int neighbor : tree[current_vertex].adjacent_vertices) {
                if (neighbor != parent_vertex) {
                    exclude_current += ComputeMaxIndependentSetRecursive(tree, neighbor, memo, current_vertex);
                }
            }

            memo[current_vertex] = std::max(include_current, exclude_current);
        }
    }

    return memo[current_vertex];
}

int ComputeMaxIndependentSet(const Tree& tree) noexcept
{
    if (tree.empty()) {
        return 0;
    }

    MemoizationTable memo(tree.size(), -1);
    return ComputeMaxIndependentSetRecursive(tree, 0, memo, -1);
}

int main()
{
    const Tree tree = ReadTree();
    std::cout << ComputeMaxIndependentSet(tree) << std::endl;
    return 0;
}