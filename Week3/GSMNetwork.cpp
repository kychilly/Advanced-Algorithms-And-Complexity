#include <iostream>
#include <vector>
#include <cstdio>

using std::cin;
using std::vector;
using std::ios;

struct GraphColoringToSAT {
    static constexpr auto colorCount = 3;

    struct Edge {
        int from;
        int to;
    };

    GraphColoringToSAT(int vertexCountInput, int edgeCountInput)
        : vertexCount(vertexCountInput)
        , graphEdges(edgeCountInput)
    {
    }

    void print_SAT_formula() const noexcept
    {
        const int clauseCount = colorCount * graphEdges.size() + vertexCount;
        const int variableCount = vertexCount * colorCount;

        printf("%d %d\n", clauseCount, variableCount);

        for (int j = 0, currentVariable = 1; j < vertexCount; ++j, currentVariable += colorCount) {
            printf("%d %d %d 0\n", currentVariable, currentVariable + 1, currentVariable + 2);
        }

        for (const Edge& currentEdge : graphEdges) {
            printf("%d %d 0\n", -((currentEdge.from - 1) * colorCount + 1), -((currentEdge.to - 1) * colorCount + 1));
            printf("%d %d 0\n", -((currentEdge.from - 1) * colorCount + 2), -((currentEdge.to - 1) * colorCount + 2));
            printf("%d %d 0\n", -((currentEdge.from - 1) * colorCount + 3), -((currentEdge.to - 1) * colorCount + 3));
        }
    }

    int vertexCount;
    vector<Edge> graphEdges;
};

int main()
{
    ios::sync_with_stdio(false);

    int vertexCountInput, edgeCountInput;
    cin >> vertexCountInput >> edgeCountInput;

    GraphColoringToSAT transformer(vertexCountInput, edgeCountInput);

    for (int i = 0; i < edgeCountInput; ++i) {
        cin >> transformer.graphEdges[i].from >> transformer.graphEdges[i].to;
    }

    transformer.print_SAT_formula();

    return 0;
}