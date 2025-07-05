#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

struct ILPToSATTransformer {
    vector<vector<int>> constraintMatrix;
    vector<int> bounds;

    ILPToSATTransformer(int numConstraints, int numVariables) : 
        constraintMatrix(numConstraints, vector<int>(numVariables)), bounds(numConstraints) {}

    void transformToSAT() {
        vector<string> clauses;
        int numVariables = constraintMatrix[0].size();

        for (int constraintIdx = 0; constraintIdx < constraintMatrix.size(); ++constraintIdx) {
            const auto &currentConstraint = constraintMatrix[constraintIdx];
            vector<int> nonZeroIndices;
            for (int varIdx = 0; varIdx < numVariables; ++varIdx) {
                if (currentConstraint[varIdx] != 0) {
                    nonZeroIndices.push_back(varIdx);
                }
            }
            int nonZeroCount = nonZeroIndices.size();

            for (int mask = 0; mask < (1 << nonZeroCount); ++mask) {
                int total = 0;
                for (int bit = 0; bit < nonZeroCount; ++bit) {
                    if (mask & (1 << bit)) {
                        total += currentConstraint[nonZeroIndices[bit]];
                    }
                }
                if (total > bounds[constraintIdx]) {
                    string clause;
                    for (int bit = 0; bit < nonZeroCount; ++bit) {
                        int varIdx = nonZeroIndices[bit];
                        if (mask & (1 << bit)) {
                            clause += to_string(-(varIdx + 1)) + " ";
                        } else {
                            clause += to_string(varIdx + 1) + " ";
                        }
                    }
                    clause += "0";
                    clauses.push_back(clause);
                }
            }
        }

        if (clauses.empty()) {
            clauses.push_back("1 -1 0");
        }

        cout << clauses.size() << " " << numVariables << endl;
        for (const auto &clause : clauses) {
            cout << clause << endl;
        }
    }
};

int main() {
    ios::sync_with_stdio(false);

    int numConstraints, numVariables;
    cin >> numConstraints >> numVariables;
    ILPToSATTransformer transformer(numConstraints, numVariables);
    for (int i = 0; i < numConstraints; i++) {
        for (int j = 0; j < numVariables; j++) {
            cin >> transformer.constraintMatrix[i][j];
        }
    }
    for (int i = 0; i < numConstraints; i++) {
        cin >> transformer.bounds[i];
    }

    transformer.transformToSAT();

    return 0;
}