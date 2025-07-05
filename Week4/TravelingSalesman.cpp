#include <iostream>
#include <vector>
#include <bitset>

constexpr int INFINITE = 1e9;

struct PathNode {
    int distance{ INFINITE };
    int predecessor{ 0 };
};

using AdjacencyMatrix = std::vector<std::vector<int>>;
using NodeSubset = std::vector<int>;
using PathCostMatrix = std::vector<std::vector<PathNode>>;
using TourPath = std::pair<int, std::vector<int>>;

AdjacencyMatrix read_input_data() noexcept
{
    std::ios::sync_with_stdio(false);

    int num_vertices, num_edges;
    std::cin >> num_vertices >> num_edges;
    AdjacencyMatrix graph(num_vertices, std::vector<int>(num_vertices, INFINITE));

    for (int i = 0, source, destination, weight; i < num_edges; ++i) {
        std::cin >> source >> destination >> weight;
        --source, --destination;
        graph[source][destination] = graph[destination][source] = weight;
    }

    return graph;
}

class TravelingSalesman {
public:
    TravelingSalesman(AdjacencyMatrix complete_graph)
        : num_nodes{ (int)complete_graph.size() }
        , path_cost{ PathCostMatrix(1 << (num_nodes), std::vector<PathNode>(num_nodes)) }
        , distance_matrix{ std::move(complete_graph) }
    {
        for (int node = 1; node < num_nodes; ++node) {
            path_cost[1 << node][node] = PathNode{ distance_matrix[0][node], 0 };
        }
    }

    TourPath find_optimal_tour() noexcept
    {
        for (int subset_size = 2; subset_size < 1 << (num_nodes - 1); ++subset_size) {
            if (subset_size & (subset_size - 1)) {
                NodeSubset current_subset = generate_subsets(subset_size);
                int bitmask = subset_size * 2;

                for (auto node_k : current_subset) {
                    int previous_bitmask = bitmask ^ (1 << node_k);
                    PathNode min_path;

                    for (auto node_m : current_subset) {
                        if (path_cost[previous_bitmask][node_m].distance + distance_matrix[node_m][node_k] < min_path.distance && node_k != node_m) {
                            min_path = PathNode{ path_cost[previous_bitmask][node_m].distance + distance_matrix[node_m][node_k], node_m };
                        }
                    }
                    path_cost[bitmask][node_k] = min_path;
                }
            }
        }

        return reconstruct_optimal_path();
    }

private:
    NodeSubset generate_subsets(const int size) noexcept
    {
        NodeSubset subsets;
        std::bitset<16> bit_representation = size;

        for (auto node = 0u; node < bit_representation.size(); ++node) {
            if (bit_representation[node]) {
                subsets.emplace_back(node + 1);
            }
        }

        return subsets;
    }

    TourPath reconstruct_optimal_path() noexcept
    {
        PathNode min_path;
        int full_bitmask = (1 << num_nodes) - 2;

        for (int node = 1; node < num_nodes; ++node) {
            if (min_path.distance > path_cost[full_bitmask][node].distance + distance_matrix[node][0]) {
                min_path = PathNode{ path_cost[full_bitmask][node].distance + distance_matrix[node][0], node };
            }
        }

        if (min_path.distance == INFINITE) {
            return { -1, {} };
        }

        TourPath optimal_tour;
        optimal_tour.second.reserve(num_nodes);
        optimal_tour.first = min_path.distance;
        optimal_tour.second.emplace_back(1);

        for (int i = 0; i < num_nodes - 1; ++i) {
            optimal_tour.second.emplace_back(min_path.predecessor + 1);
            min_path.predecessor = path_cost[full_bitmask][min_path.predecessor].predecessor;
            full_bitmask = full_bitmask ^ (1 << (optimal_tour.second.back() - 1));
        }

        return optimal_tour;
    }

    const int num_nodes;
    PathCostMatrix path_cost;
    AdjacencyMatrix distance_matrix;
};

void print_solution(const TourPath& solution) noexcept
{
    printf("%d\n", solution.first);
    if (solution.first == -1)
        return;
    for (auto node : solution.second)
        printf("%d ", node);
    printf("\n");
}

int main()
{
    TravelingSalesman tsp_solver(read_input_data());
    print_solution(tsp_solver.find_optimal_tour());
    return 0;
}