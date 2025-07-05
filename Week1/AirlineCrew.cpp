#include <iostream>
#include <vector>
#include <cstring>

using namespace std;

const int MAXN = 100;

int leftSize, rightSize;
int matchLeft[MAXN], matchRight[MAXN];
int visited[MAXN];

bool find_augmenting_path(int u, const vector<vector<bool>>& adjacencyMatrix) {
    for (int v = 0; v < rightSize; ++v) {
        if (adjacencyMatrix[u][v] && !visited[v]) {
            visited[v] = 1;

            if (matchRight[v] == -1 || find_augmenting_path(matchRight[v], adjacencyMatrix)) {
                matchLeft[u] = v;
                matchRight[v] = u;
                return true;
            }
        }
    }
    return false;
}

int compute_maximum_matching(const vector<vector<bool>>& adjacencyMatrix) {
    memset(matchLeft, -1, sizeof(matchLeft));
    memset(matchRight, -1, sizeof(matchRight));

    int matchCount = 0;

    for (int u = 0; u < leftSize; ++u) {
        if (matchLeft[u] == -1) {
            memset(visited, 0, sizeof(visited));
            if (find_augmenting_path(u, adjacencyMatrix))
                ++matchCount;
        }
    }

    return matchCount;
}

vector<vector<bool>> read_input_matrix() {
    int nLeft, nRight;
    cin >> nLeft >> nRight;
    vector<vector<bool>> matrix(nLeft, vector<bool>(nRight));

    for (int i = 0; i < nLeft; ++i)
        for (int j = 0; j < nRight; ++j) {
            int value;
            cin >> value;
            matrix[i][j] = (value == 1);
        }

    return matrix;
}

int main() {
    vector<vector<bool>> adjacencyMatrix = read_input_matrix();

    leftSize = adjacencyMatrix.size();
    rightSize = adjacencyMatrix[0].size();

    int matchCount = compute_maximum_matching(adjacencyMatrix);

    for (int i = 0; i < leftSize; ++i) {
        if (matchLeft[i] != -1)
            cout << matchLeft[i] + 1 << ' ';
        else
            cout << -1 << ' ';
    }

    return 0;
}
