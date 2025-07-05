#include <algorithm>
#include <iostream>
#include <vector>
#include <cstdio>
#include <functional>
#include <limits>
#include <functional>
#include <valarray>
#include <iomanip>

const long double EPS = std::numeric_limits<long double>::epsilon();

enum class SolutionState { optimal,
                          infeasible,
                          unbounded };
enum class MethodPhase { phase1,
                         phase2 };

using namespace std;

using Matrix = vector<vector<long double>>;
using Column = vector<long double>;
using Row = vector<long double>;

struct PivotPosition {
    short row;
    short column;
    bool is_optimal() { return row == -1 || column == -1; }
};

struct SimplexSolver {

    void debug_print() const
    {
        for (std::size_t i = 0; i < constraint_matrix.size(); ++i) {
            for (std::size_t j = 0; j < constraint_matrix[i].size(); ++j) {
                std::cout << std::fixed << std::setw(5) << std::setprecision(2) << constraint_matrix[i][j] << ' ';
            }
            std::cout << fixed << std::setw(5) << " | " << std::setprecision(2) << rhs_values[i] << std::endl;
        }

        for (const auto& v : objective_coeffs) {
            std::cout << std::setw(7) << v << ' ';
        }

        std::cout << " | " << rhs_values[rhs_values.size() - 2] << std::endl;

        if (current_phase == MethodPhase::phase1) {
            for (const auto& v : artificial_obj_coeffs) {
                cout << std::setw(5) << v << ' ';
            }

            std::cout << " | " << rhs_values.back() << endl;
        }
    }

    void handle_artificial_vars()
    {
        for (std::size_t i = 0, j = num_original_vars + num_constraints; i < rhs_values.size() - 1; ++i, ++j) {
            if (rhs_values[i] < 0.0) {
                basic_vars[i] = -2;
                constraint_matrix[i][j] = -1;

                rhs_values.back() += rhs_values[i];
                rhs_values[i] = -rhs_values[i];

                std::transform(constraint_matrix[i].begin(), constraint_matrix[i].end(),
                              constraint_matrix[i].begin(), std::negate<long double>());

                for (int k = 0; k < num_original_vars + num_constraints; ++k) {
                    artificial_obj_coeffs[k] += constraint_matrix[i][k];
                }
            }
        }

        std::transform(artificial_obj_coeffs.begin(), artificial_obj_coeffs.end(),
                       artificial_obj_coeffs.begin(), std::negate<long double>());
    }

    void handle_slack_vars()
    {
        for (std::size_t i = 0; i < constraint_matrix.size(); ++i) {
            constraint_matrix[i][i + num_original_vars] = 1;
        }
    }

    void perform_gauss_jordan_eliminations(Row& objective_row)
    {
        while (true) {
            PivotPosition pivot = find_pivot(objective_row);

            if (pivot.is_optimal() || current_state == SolutionState::unbounded) {
                break;
            }

            basic_vars[pivot.column] = pivot.row;

            scale_pivot(pivot);

            process_pivot(pivot, objective_row);
        }
    }

    void remove_artificial_vars()
    {
        objective_coeffs.resize(objective_coeffs.size() - num_constraints);
        rhs_values.pop_back();

        for (auto& row : constraint_matrix) {
            row.resize(row.size() - num_constraints);
        }
    }

    void run_phase2()
    {
        current_phase = MethodPhase::phase2;
        remove_artificial_vars();
        perform_gauss_jordan_eliminations(objective_coeffs);
    }

    void run_phase1()
    {
        current_phase = MethodPhase::phase1;
        perform_gauss_jordan_eliminations(artificial_obj_coeffs);
        current_state = is_zero(rhs_values.back()) ? SolutionState::optimal : SolutionState::infeasible;
    }

    void initialize_tableau()
    {
        basic_vars = vector<int>(constraint_matrix[0].size(), -1);
        std::transform(objective_coeffs.begin(), objective_coeffs.end(),
                       objective_coeffs.begin(), std::negate<double>());
        handle_slack_vars();
        if (has_negative_rhs) {
            artificial_obj_coeffs = Row(objective_coeffs.size());
            handle_artificial_vars();
        }
    }

    bool is_equal(double a, double b, double epsilon = 0.001)
    {
        return std::abs(a - b) < epsilon;
    }

    bool has_negative_constraints() const
    {
        auto it = std::find_if(rhs_values.cbegin(), rhs_values.cend(), [](auto j) { return j < 0.0; });
        return it == rhs_values.cend() ? false : true;
    }

    pair<int, vector<long double>> get_solution()
    {
        vector<long double> solution(num_original_vars);

        for (int i = 0; i < num_original_vars; ++i) {
            long double sum = 0.0;
            int k = 0;
            for (std::size_t j = 0; j < constraint_matrix.size(); ++j) {
                if (basic_vars[j] >= 0.0)
                    sum += fabs(constraint_matrix[j][i]);
                if (is_equal(constraint_matrix[j][i], 1.0)) {
                    k = j;
                }
            }

            solution[i] = (sum - 1.0 > EPS) ? 0.0 : rhs_values[k];
        }

        return { 0, std::move(solution) };
    }

    pair<int, vector<long double>> solve()
    {
        has_negative_rhs = has_negative_constraints();

        initialize_tableau();

        if (has_negative_rhs) {
            run_phase1();
            if (current_state == SolutionState::infeasible) {
                return { -1, {} };
            }
        }

        run_phase2();

        if (current_state == SolutionState::unbounded) {
            return { 1, {} };
        }

        return get_solution();
    }

    bool is_zero(long double a, long double epsilon = 0.001)
    {
        return std::abs(a - 0.0) < epsilon;
    }

    PivotPosition find_pivot(Row& objective_row)
    {
        short i = 0, j = distance(objective_row.begin(), min_element(objective_row.begin(), objective_row.end()));
        long double ratio = numeric_limits<long double>::max() - 1;

        if (objective_row[j] >= 0.0) {
            return { -1, -1 };
        }

        for (std::size_t k = 0; k < constraint_matrix.size(); ++k) {
            long double r = rhs_values[k] / constraint_matrix[k][j];
            if ((constraint_matrix[k][j] > 0.0 || constraint_matrix[k][j] < 0.0) && r - ratio < EPS && r > 0.0) {
                ratio = r;
                i = k;
            }
        }

        if (ratio == numeric_limits<long double>::max() - 1) {
            current_state = SolutionState::unbounded;
        }

        return { i, j };
    }

    void process_pivot(PivotPosition pivot, Row& objective_row)
    {
        long double multiplier{ 0.0 };

        for (int i = 0, s = constraint_matrix.size(); i < s; ++i) {
            if (pivot.row != i && !is_zero(constraint_matrix[i][pivot.column], EPS)) {

                multiplier = constraint_matrix[i][pivot.column];

                for (std::size_t j = 0; j < constraint_matrix[i].size(); ++j) {
                    constraint_matrix[i][j] -= constraint_matrix[pivot.row][j] * multiplier;
                }

                rhs_values[i] -= rhs_values[pivot.row] * multiplier;
            }
        }

        if (current_phase == MethodPhase::phase1) {
            rhs_values[rhs_values.size() - 2] -= rhs_values[pivot.row] * objective_coeffs[pivot.column];
            rhs_values[rhs_values.size() - 1] -= rhs_values[pivot.row] * objective_row[pivot.column];
            auto obj_mult = objective_row[pivot.column];
            auto coeff_mult = objective_coeffs[pivot.column];
            for (std::size_t i = 0; i < objective_row.size(); ++i) {
                objective_row[i] -= constraint_matrix[pivot.row][i] * obj_mult;
                objective_coeffs[i] -= constraint_matrix[pivot.row][i] * coeff_mult;
            }
        }
        else {
            rhs_values[rhs_values.size() - 1] -= rhs_values[pivot.row] * objective_coeffs[pivot.column];
            auto coeff_mult = objective_coeffs[pivot.column];
            for (std::size_t i = 0; i < objective_row.size(); ++i) {
                objective_coeffs[i] -= constraint_matrix[pivot.row][i] * coeff_mult;
            }
        }
    }

    void scale_pivot(PivotPosition pivot)
    {
        auto divisor = constraint_matrix[pivot.row][pivot.column];
        rhs_values[pivot.row] /= divisor;
        for (auto& val : constraint_matrix[pivot.row]) {
            val /= divisor;
        }
    }

    int num_constraints, num_original_vars;
    Matrix constraint_matrix;
    vector<long double> rhs_values, objective_coeffs, artificial_obj_coeffs;
    vector<int> basic_vars;
    SolutionState current_state;
    MethodPhase current_phase;
    bool has_negative_rhs;
};

int main()
{
    std::ios_base::sync_with_stdio(false);

    int num_constraints, num_vars;
    cin >> num_constraints >> num_vars;

    Matrix A(num_constraints, vector<long double>(num_constraints + num_vars + num_constraints, 0.0));
    for (int i = 0; i < num_constraints; i++) {
        for (int j = 0; j < num_vars; j++) {
            cin >> A[i][j];
        }
    }
    vector<long double> b(num_constraints + 2);
    for (int i = 0; i < num_constraints; i++) {
        cin >> b[i];
    }

    vector<long double> c(num_constraints + num_vars + num_constraints);
    for (int i = 0; i < num_vars; i++) {
        cin >> c[i];
    }

    SimplexSolver solver{ num_constraints, num_vars, std::move(A), std::move(b), std::move(c) };

    pair<int, vector<long double>> ans = solver.solve();

    switch (ans.first) {
    case -1:
        printf("No solution\n");
        break;
    case 0:
        printf("Bounded solution\n");
        for (int i = 0; i < num_vars; i++) {
            printf("%.18Lf%c", ans.second[i], " \n"[i + 1 == num_vars]);
        }
        break;
    case 1:
        printf("Infinity\n");
        break;
    }
    return 0;
}