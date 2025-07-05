#include <iostream>
#include <vector>
#include <queue>
#include <limits>

using std::vector;
using std::queue;
using std::size_t;

class FlowGraph {
public:
    struct Edge {
        int from, to, capacity, flow;
    };

private:
    vector<Edge> edges;
    vector<vector<size_t>> adjacencyList;

public:
    explicit FlowGraph(size_t vertexCount) : adjacencyList(vertexCount) {
        edges.reserve(vertexCount / 2);
    }

    void add_edge(int from, int to, int capacity) {
        Edge forward = { from, to, capacity, 0 };
        Edge backward = { to, from, 0, 0 };

        adjacencyList[from].push_back(edges.size());
        edges.push_back(forward);

        adjacencyList[to].push_back(edges.size());
        edges.push_back(backward);
    }

    size_t size() const {
        return adjacencyList.size();
    }

    const vector<size_t>& get_edge_ids(int from) const {
        return adjacencyList[from];
    }

    const Edge& get_edge(size_t id) const {
        return edges[id];
    }

    void add_flow(size_t id, int flowAmount) {
        edges[id].flow += flowAmount;
        edges[id ^ 1].flow -= flowAmount; // ^1 toggles between forward/backward
    }
};

// Reads stock chart data: stockCount x timePoints
vector<vector<int>> read_stock_data(size_t stockCount, size_t timePoints) {
    vector<vector<int>> stockData(stockCount, vector<int>(timePoints));
    for (size_t i = 0; i < stockCount; ++i)
        for (size_t j = 0; j < timePoints; ++j)
            std::cin >> stockData[i][j];
    return stockData;
}

// Builds a bipartite DAG from stock chart data
FlowGraph build_flow_graph(size_t stockCount, size_t timePoints) {
    vector<vector<int>> stockData = read_stock_data(stockCount, timePoints);
    FlowGraph graph(stockCount * 2 + 2);

    int source = 0;
    int sink = stockCount * 2 + 1;

    // Connect source to left side of bipartite graph
    for (int i = 0; i < stockCount; ++i)
        graph.add_edge(source, i + 1, 1);

    // Connect left to right side of bipartite graph based on chart dominance
    for (int i = 0; i < stockCount; ++i) {
        for (int j = 0; j < stockCount; ++j) {
            if (i == j) continue;

            bool strictlyLess = true;
            for (int k = 0; k < timePoints; ++k) {
                if (stockData[i][k] >= stockData[j][k]) {
                    strictlyLess = false;
                    break;
                }
            }

            if (strictlyLess) {
                graph.add_edge(i + 1, stockCount + j + 1, 1);
            }
        }
    }

    // Connect right side of bipartite graph to sink
    for (int i = stockCount + 1; i <= stockCount * 2; ++i)
        graph.add_edge(i, sink, 1);

    return graph;
}

// Performs BFS to find augmenting path in residual graph
void bfs_augmenting_path(const FlowGraph& graph, int source, int sink, vector<int>& parentEdgeId) {
    queue<int> q;
    q.push(source);
    std::fill(parentEdgeId.begin(), parentEdgeId.end(), -1);

    while (!q.empty()) {
        int current = q.front();
        q.pop();

        for (auto edgeId : graph.get_edge_ids(current)) {
            const FlowGraph::Edge& edge = graph.get_edge(edgeId);

            if (parentEdgeId[edge.to] == -1 && edge.capacity > edge.flow && edge.to != source) {
                parentEdgeId[edge.to] = edgeId;
                q.push(edge.to);
            }
        }
    }
}

// Runs Edmonds-Karp to compute maximum flow
void compute_max_flow(FlowGraph& graph, int source, int sink) {
    vector<int> parentEdgeId(graph.size());

    while (true) {
        bfs_augmenting_path(graph, source, sink, parentEdgeId);

        if (parentEdgeId[sink] == -1)
            break;

        int pathFlow = std::numeric_limits<int>::max();

        // Trace minimum residual capacity along path
        for (int edgeId = parentEdgeId[sink]; edgeId != -1; edgeId = parentEdgeId[graph.get_edge(edgeId).from]) {
            const FlowGraph::Edge& edge = graph.get_edge(edgeId);
            pathFlow = std::min(pathFlow, edge.capacity - edge.flow);
        }

        // Augment flow along the path
        for (int edgeId = parentEdgeId[sink]; edgeId != -1; edgeId = parentEdgeId[graph.get_edge(edgeId).from]) {
            graph.add_flow(edgeId, pathFlow);
        }
    }
}

// Computes the minimum number of overlaid charts (minimum path cover in DAG)
int compute_minimum_overlaid_charts(const FlowGraph& graph, int stockCount) {
    int matched = 0;

    for (int i = 1; i <= stockCount; ++i) {
        for (auto edgeId : graph.get_edge_ids(i)) {
            const FlowGraph::Edge& edge = graph.get_edge(edgeId);
            if (edge.flow > 0) {
                ++matched;
                break;
            }
        }
    }

    return stockCount - matched;
}

// Main entry point
int main() {
    std::ios_base::sync_with_stdio(false);

    size_t stockCount, timePoints;
    std::cin >> stockCount >> timePoints;

    FlowGraph graph = build_flow_graph(stockCount, timePoints);

    int source = 0;
    int sink = graph.size() - 1;

    compute_max_flow(graph, source, sink);

    int result = compute_minimum_overlaid_charts(graph, stockCount);
    std::cout << result << std::endl;

    return 0;
}
