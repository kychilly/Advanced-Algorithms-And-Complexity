#include <iostream>
#include <vector>
#include <cstdio>

using std::cin;
using std::vector;
using std::ios;

struct HamiltonianPathToSAT {

    HamiltonianPathToSAT(int vertexCountInput, int edgeCountInput)
        : clauseCount(0)
        , vertexCount(vertexCountInput)
        , adjacencyMatrix(vertexCountInput, vector<bool>(vertexCountInput, false))
        , variableMapping(vertexCountInput, vector<int>(vertexCountInput))
    {
        for (int i = 0, currentVariable = 1; i < vertexCount; ++i) {
            for (int j = 0; j < vertexCount; ++j) {
                variableMapping[i][j] = currentVariable++;
            }
        }
    }

    void print_SAT_formula(const int maxClauses) noexcept
    {
        clausesOutput.reserve(maxClauses * 3);

        eachVertexInPath();
        eachVertexInPathOnlyOnce();
        eachPathPositionOccupied();
        verticesOccupyDifferentPositions();
        nonadjacentVerticesNonadjacentInPath();

        printf("%d %d \n%s", clauseCount, vertexCount * vertexCount, clausesOutput.c_str());
    }

    void eachVertexInPath() noexcept
    {
        for (int i = 0; i < vertexCount; ++i, ++clauseCount) {
            for (int j = 0; j < vertexCount; ++j) {
                clausesOutput += std::to_string(variableMapping[i][j]) + " ";
            }
            clausesOutput += "0\n";
        }
    }

    void eachVertexInPathOnlyOnce() noexcept
    {
        for (int i = 0; i < vertexCount; ++i) {
            for (int j = 0; j < vertexCount; ++j) {
                for (int k = i + 1; k < vertexCount; ++k, ++clauseCount) {
                    clausesOutput += std::to_string(-variableMapping[i][j]) + " " + std::to_string(-variableMapping[k][j]) + " 0\n";
                }
            }
        }
    }

    void eachPathPositionOccupied() noexcept
    {
        for (int i = 0; i < vertexCount; ++i, ++clauseCount) {
            for (int j = 0; j < vertexCount; ++j) {
                clausesOutput += std::to_string(variableMapping[j][i]) + " ";
            }
            clausesOutput += "0\n";
        }
    }

    void verticesOccupyDifferentPositions() noexcept
    {
        for (int i = 0; i < vertexCount; ++i) {
            for (int j = 0; j < vertexCount; ++j) {
                for (int k = j + 1; k < vertexCount; ++k, ++clauseCount) {
                    clausesOutput += std::to_string(-variableMapping[i][j]) + " " + std::to_string(-variableMapping[i][k]) + " 0\n";
                }
            }
        }
    }

    void nonadjacentVerticesNonadjacentInPath() noexcept
    {
        for (int i = 0; i < vertexCount; ++i) {
            for (int j = 0; j < vertexCount; ++j) {
                if (!adjacencyMatrix[i][j] && j != i) {
                    for (int k = 0; k < vertexCount - 1; ++k, ++clauseCount) {
                        clausesOutput += std::to_string(-variableMapping[i][k]) + " " + std::to_string(-variableMapping[j][k + 1]) + " 0\n";
                    }
                }
            }
        }
    }

    unsigned int clauseCount;
    const unsigned int vertexCount;
    vector<vector<bool>> adjacencyMatrix;
    vector<vector<int>> variableMapping;
    std::string clausesOutput;
};

int main()
{
    ios::sync_with_stdio(false);

    int vertexCountInput, edgeCountInput;
    cin >> vertexCountInput >> edgeCountInput;

    HamiltonianPathToSAT transformer(vertexCountInput, edgeCountInput);

    for (int edgeIndex = 0; edgeIndex < edgeCountInput; ++edgeIndex) {
        int vertexFrom, vertexTo;
        cin >> vertexFrom >> vertexTo;
        transformer.adjacencyMatrix[vertexFrom - 1][vertexTo - 1] = true;
        transformer.adjacencyMatrix[vertexTo - 1][vertexFrom - 1] = true;
    }

    transformer.print_SAT_formula(120000);

    return 0;
}