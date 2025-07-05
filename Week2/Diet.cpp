#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include <algorithm>
#include <limits>
#include <numeric>

using namespace std;

const long double EPSILON = 1e-3;
const long double INF_LIMIT = 1e9L;

using Vector = vector<long double>;
using Matrix = vector<Vector>;

struct LinearSystem {
    Matrix coefficients;
    Vector rhs;
};

struct PivotPosition {
    short col = 0;
    short row = 0;
};

PivotPosition find_pivot(const Matrix& coeffs, vector<bool>& used_rows, vector<bool>& used_cols) {
    PivotPosition pivot;

    while (used_rows[pivot.row]) ++pivot.row;
    while (used_cols[pivot.col]) ++pivot.col;

    long double max_val = 0.0;
    for (int i = pivot.row; i < coeffs.size(); ++i) {
        if (fabs(max_val) - fabs(coeffs[i][pivot.col]) < EPSILON) {
            max_val = coeffs[i][pivot.col];
            pivot.row = i;
        }
    }

    return pivot;
}

void swap_rows(Matrix& coeffs, Vector& rhs, vector<bool>& used_rows, PivotPosition& pivot) {
    swap(coeffs[pivot.col], coeffs[pivot.row]);
    swap(rhs[pivot.col], rhs[pivot.row]);
    swap(used_rows[pivot.col], used_rows[pivot.row]);
    pivot.row = pivot.col;
}

int perform_back_substitution(Matrix& coeffs, Vector& rhs) {
    for (int i = coeffs.size() - 1; i > 0; --i) {
        long double value = rhs[i];
        for (int j = 0; j < i; ++j) {
            rhs[j] -= coeffs[j][i] * value;
        }
    }
    return 0;
}

void normalize_pivot_row(Matrix& coeffs, Vector& rhs, const PivotPosition& pivot) {
    long double divisor = coeffs[pivot.row][pivot.col];
    int n = coeffs.size();

    for (int j = pivot.col; j < n; ++j)
        coeffs[pivot.row][j] /= divisor;

    rhs[pivot.row] /= divisor;
}

void eliminate_column(Matrix& coeffs, Vector& rhs, const PivotPosition& pivot) {
    int size = coeffs.size();
    normalize_pivot_row(coeffs, rhs, pivot);

    for (int i = pivot.row + 1; i < size; ++i) {
        long double multiplier = coeffs[i][pivot.col];
        for (int j = pivot.col; j < size; ++j)
            coeffs[i][j] -= coeffs[pivot.row][j] * multiplier;
        rhs[i] -= rhs[pivot.row] * multiplier;
    }
}

void mark_used(const PivotPosition& pivot, vector<bool>& used_rows, vector<bool>& used_cols) {
    used_rows[pivot.row] = true;
    used_cols[pivot.col] = true;
}

pair<int, Vector> solve_system(LinearSystem system) {
    Matrix& coeffs = system.coefficients;
    Vector& rhs = system.rhs;

    if (coeffs.empty()) return {};

    int n = coeffs.size();
    vector<bool> used_rows(n, false), used_cols(n, false);

    for (int step = n; step > 0; --step) {
        PivotPosition pivot = find_pivot(coeffs, used_rows, used_cols);
        if (coeffs[pivot.row][pivot.col] == 0.0) break;

        swap_rows(coeffs, rhs, used_rows, pivot);
        eliminate_column(coeffs, rhs, pivot);
        mark_used(pivot, used_rows, used_cols);
    }

    int result = perform_back_substitution(coeffs, rhs);
    return { result, move(rhs) };
}

vector<vector<int>> generate_subsets(const vector<int>& base_set, int size_required) {
    int n = 1 << base_set.size();
    vector<vector<int>> subsets;

    for (int mask = 0; mask < n; ++mask) {
        vector<int> subset;
        for (int j = 0; j < base_set.size(); ++j) {
            if (mask & (1 << j))
                subset.push_back(base_set[j]);
            if (subset.size() > size_required)
                break;
        }
        if (subset.size() == size_required)
            subsets.push_back(move(subset));
    }
    return subsets;
}

vector<Vector> solve_all_subsystems(int variable_count, const Matrix& A, const Vector& b) {
    vector<Vector> all_solutions;
    vector<int> row_indices(A.size());
    iota(row_indices.begin(), row_indices.end(), 0);

    auto combinations = generate_subsets(row_indices, variable_count);
    if (A.size() == 1) combinations.emplace_back(vector<int>{0});

    for (const auto& subset : combinations) {
        Matrix sub_matrix;
        Vector sub_rhs;
        for (int idx : subset) {
            sub_matrix.push_back(A[idx]);
            sub_rhs.push_back(b[idx]);
        }
        auto result = solve_system({ sub_matrix, sub_rhs });
        int status = result.first;
        Vector solution = std::move(result.second);

        if (status == 0 && !solution.empty())
            all_solutions.push_back(move(solution));
    }
    return all_solutions;
}

void augment_constraints(int& constraint_count, int var_count, Matrix& A, Vector& b) {
    while (constraint_count < var_count) {
        A.emplace_back(var_count, 0);
        b.push_back(0);
        ++constraint_count;
    }

    A.emplace_back(Vector(var_count, 1.0));
    b.push_back(INF_LIMIT);
    ++constraint_count;

    for (int i = 0; i < var_count; ++i) {
        Vector constraint(var_count, 0.0);
        constraint[i] = -1.0;
        A.push_back(move(constraint));
        b.push_back(0.0);
        ++constraint_count;
    }
}

pair<int, Vector> solve_diet(int n, int m, Matrix A, Vector b, Vector c) {
    augment_constraints(n, m, A, b);

    auto candidate_solutions = solve_all_subsystems(m, A, b);
    if (candidate_solutions.empty()) return { -1, {} };

    long double max_utility = -numeric_limits<long double>::max();
    int best_idx = -1;

    for (int i = 0; i < candidate_solutions.size(); ++i) {
        auto& solution = candidate_solutions[i];
        bool valid = true;

        for (int j = 0; j < n && valid; ++j) {
            long double lhs = inner_product(A[j].begin(), A[j].end(), solution.begin(), 0.0L);
            if (lhs - b[j] > EPSILON) valid = false;
        }

        if (valid) {
            long double value = inner_product(c.begin(), c.end(), solution.begin(), 0.0L);
            if (value > max_utility) {
                max_utility = value;
                best_idx = i;
            }
        }
    }

    if (best_idx == -1) return { -1, {} };

    auto& best_solution = candidate_solutions[best_idx];
    if (accumulate(best_solution.begin(), best_solution.end(), 0.0L) + EPSILON >= INF_LIMIT)
        return { 1, {} };

    return { 0, move(best_solution) };
}

int main() {
    ios::sync_with_stdio(false);

    int n, m;
    cin >> n >> m;

    Matrix A(n, Vector(m));
    for (auto& row : A)
        for (auto& val : row)
            cin >> val;

    Vector b(n), c(m);
    for (auto& bi : b) cin >> bi;
    for (auto& ci : c) cin >> ci;

    auto result_pair = solve_diet(n, m, move(A), move(b), move(c));
    int result_code = result_pair.first;
    Vector solution = std::move(result_pair.second);


    switch (result_code) {
        case -1:
            cout << "No solution\n";
            break;
        case 1:
            cout << "Infinity\n";
            break;
        case 0:
            cout << "Bounded solution\n";
            for (int i = 0; i < m; ++i)
                printf("%.18Lf%c", solution[i], " \n"[i + 1 == m]);
            break;
    }

    return 0;
}
