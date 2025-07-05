#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>

constexpr int OUTPUT_PRECISION = 10;

using Vector = std::vector<double>;
using Matrix = std::vector<Vector>;

struct LinearSystem {
    Matrix matrix;
    Vector rhs;
};

struct PivotPosition {
    short col = 0;
    short row = 0;
};

LinearSystem read_system() {
    int dimension;
    std::cin >> dimension;

    Matrix matrix(dimension, Vector(dimension, 0.0));
    Vector rhs(dimension, 0.0);

    for (int row = 0; row < dimension; ++row) {
        for (int col = 0; col < dimension; ++col)
            std::cin >> matrix[row][col];
        std::cin >> rhs[row];
    }

    return { std::move(matrix), std::move(rhs) };
}

PivotPosition choose_pivot(const Matrix& matrix, const std::vector<bool>& used_rows, const std::vector<bool>& used_cols) {
    PivotPosition pivot;

    while (used_rows[pivot.row]) ++pivot.row;
    while (used_cols[pivot.col]) ++pivot.col;

    double max_abs_val = 0.0;
    for (int i = pivot.row; i < matrix.size(); ++i) {
        if (std::fabs(matrix[i][pivot.col]) > std::fabs(max_abs_val)) {
            max_abs_val = matrix[i][pivot.col];
            pivot.row = i;
        }
    }

    return pivot;
}

void swap_rows(Matrix& matrix, Vector& rhs, std::vector<bool>& used_rows, PivotPosition& pivot) {
    std::swap(matrix[pivot.col], matrix[pivot.row]);
    std::swap(rhs[pivot.col], rhs[pivot.row]);
    std::swap(used_rows[pivot.col], used_rows[pivot.row]);
    pivot.row = pivot.col;
}

void scale_pivot_row(Matrix& matrix, Vector& rhs, const PivotPosition& pivot) {
    double pivot_val = matrix[pivot.row][pivot.col];
    int size = matrix.size();

    for (int j = pivot.col; j < size; ++j)
        matrix[pivot.row][j] /= pivot_val;

    rhs[pivot.row] /= pivot_val;
}

void eliminate_below(Matrix& matrix, Vector& rhs, const PivotPosition& pivot) {
    int size = matrix.size();
    for (int i = pivot.row + 1; i < size; ++i) {
        double factor = matrix[i][pivot.col];
        for (int j = pivot.col; j < size; ++j)
            matrix[i][j] -= matrix[pivot.row][j] * factor;
        rhs[i] -= rhs[pivot.row] * factor;
    }
}

void mark_pivot(const PivotPosition& pivot, std::vector<bool>& used_rows, std::vector<bool>& used_cols) {
    used_rows[pivot.row] = true;
    used_cols[pivot.col] = true;
}

void back_substitute(Matrix& matrix, Vector& rhs) {
    for (int i = matrix.size() - 1; i > 0; --i) {
        double value = rhs[i];
        for (int j = 0; j < i; ++j) {
            rhs[j] -= matrix[j][i] * value;
            matrix[j][i] = 0.0;
        }
    }
}

Vector solve_system(LinearSystem system) {
    Matrix& matrix = system.matrix;
    Vector& rhs = system.rhs;

    int n = matrix.size();
    std::vector<bool> used_rows(n, false);
    std::vector<bool> used_cols(n, false);

    for (int step = 0; step < n; ++step) {
        PivotPosition pivot = choose_pivot(matrix, used_rows, used_cols);
        swap_rows(matrix, rhs, used_rows, pivot);
        scale_pivot_row(matrix, rhs, pivot);
        eliminate_below(matrix, rhs, pivot);
        mark_pivot(pivot, used_rows, used_cols);
    }

    back_substitute(matrix, rhs);
    return rhs;
}

void print_solution(const Vector& solution) {
    std::cout << std::fixed << std::setprecision(OUTPUT_PRECISION);
    for (double val : solution)
        std::cout << val << ' ';
    std::cout << '\n';
}

int main() {
    std::ios_base::sync_with_stdio(false);

    LinearSystem system = read_system();

    if (!system.matrix.empty()) {
        Vector solution = solve_system(std::move(system));
        print_solution(solution);
    }

    return 0;
}
