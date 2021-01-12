#pragma once

#include <map>
#include <vector>
#include <iostream>
#include <algorithm>
#include <chrono>

typedef std::vector<double> Vector;
typedef std::vector<std::vector<double>> Matrix;

enum SolveResult {
    Infeasible  = -1,
    Solvable    = 1,
    Unbounded   = 0
};

class Simplex
{
private:
    constexpr static double eps = 1e-9;
    int sgn(double x) {
        return x < -eps ? -1 : (x > eps ? 1 : 0);
    }
    void pivot(int n, int m, Vector &c, Matrix &A, Vector &b, double &ans, int pivot_i, int pivot_j);
public:
    Simplex();
    ~Simplex();
    SolveResult solveStandardForm(int n, int m, Vector &c, Matrix &A, Vector &b, double &ans, Vector &x, double &ans_dual, Vector &x_dual, bool &dual_infeasible, double &t_sp, double &t_dsp);
    SolveResult solveNonStandardForm(int n, int m, Vector &c, Matrix &A, Vector &b, std::vector<int> &d, std::vector<int> &e, double &ans, Vector &x, double &ans_dual, Vector &x_dual, bool &dual_infeasible, double &t_sp, double &t_dsp);
    SolveResult dualSimplexMethod(int n, int m, Vector c, Matrix A, Vector b, double &ans, Vector &x);
    SolveResult simplexMethod(int n, int m, Vector c, Matrix A, Vector b, double &ans, Vector &x);
};

Simplex::Simplex()
{
}

Simplex::~Simplex()
{
}

/*
 * Solve linear programming in standard form:
 * Input: c_n, A_{m*n}, b_m
 * Output: ans, x_n 
 * Maximize c^T dot x
 * s.t. A dot x <= b, x_i >= 0
 * it will add slack variables and solve LP using simplex method and dual-simplex method
 */
SolveResult Simplex::solveStandardForm(int n, int m, Vector &c, Matrix &A, Vector &b, double &ans, Vector &x, double &ans_dual, Vector &x_dual, bool &dual_infeasible, double &t_sp, double &t_dsp) {
    if (n != c.size() || m != A.size() || m != b.size()) {
        std::cerr << "Wrong size of params" << std::endl;
        return Infeasible;
    }

    // Convert all constraints to equations with slack variables, and then write the problem as a tableau with some negative right sides, with or without a Z-COLUMN.
    int n_0 = n;
    for (auto &Ai : A) {
        if (n != Ai.size()) {
            std::cerr << "Wrong size of params" << std::endl;
            return Infeasible;
        }
        for (int j = n; j < n + m; j++) {
            Ai.push_back(0);
        }
        Ai[n_0++] = 1;
    }

    for (auto &ci : c) {
        ci = -ci;
    }
    for (int j = n; j < n + m; j++) {
        c.push_back(0);
    }

    auto start0 = std::chrono::system_clock::now();
    SolveResult dual_rst = dualSimplexMethod(n_0, m, c, A, b, ans_dual, x_dual);
    auto end0 = std::chrono::system_clock::now();
    auto t0 = std::chrono::duration_cast<std::chrono::milliseconds>(end0 - start0);
    t_dsp += t0.count();

    if (dual_rst == Infeasible) {
        dual_infeasible = true;
    }

    auto start1 = std::chrono::system_clock::now();
    SolveResult simplex_rst = simplexMethod(n_0, m, c, A, b, ans, x);
    auto end1 = std::chrono::system_clock::now();
    auto t1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1);
    t_sp = t1.count();

    return simplex_rst;
}


/*
 * Pivoting operation:
 * uses row operations to change one matrix entry (the PIVOT) to "1" 
 * and then to change all other entries in the pivot's column into ZERO's.
 * 
 */

void Simplex::pivot(int n, int m, Vector &c, Matrix &A, Vector &b, double &ans, int pivot_i, int pivot_j) {
    auto &A_p = A[pivot_i];
    
    // change the pivot to "1"
    double factor = A_p[pivot_j];
    for (int j = 0; j < n; j++) {
        A_p[j] /= factor;
    }
    b[pivot_i] /= factor;
    
    // change all other entries in the pivot's column into ZERO's
    for (int i = 0; i < m; i++) {
        if (i == pivot_i || sgn(A[i][pivot_j]) == 0) {
            continue;
        }
        factor = A[i][pivot_j];
        for (int j = 0; j < n; j++) {
            A[i][j] -= factor * A_p[j];
        }
        b[i] -= factor * b[pivot_i];
    }

    factor = c[pivot_j];
    for (int j = 0; j < n; j++) {
        c[j] -= factor * A_p[j];
    }
    ans -= factor * b[pivot_i];
}

/*
 * Solve linear programming in non-standard form:
 * Input: c_n, A_{m*n}, d_m, b_m, e_n
 * Output: ans, x_n 
 * Minimize c^T dot x
 * s.t. A dot x <= or == or >= b, x_i <= or == or >= 0
 * it will transform LP into standard form and solve it
 */

SolveResult Simplex::solveNonStandardForm(int n, int m, Vector &c, Matrix &A, Vector &b, std::vector<int> &d, std::vector<int> &e, double &ans, Vector &x, double &ans_dual, Vector &x_dual, bool &dual_infeasible, double &t_sp, double &t_dsp) {
    // inverse c
    for (int i = 0; i < n; i++) {
        c[i] = -c[i];
    }

    // turn "==" constraint to two "<=" constraints
    // turn ">=" constraint to one "<=" constrint
    int m_0 = m;
    for (int i = 0; i < m; i++) {
        if (d[i] == 1) {
            for (int j = 0; j < n; j++) {
                A[i][j] = -A[i][j];
            }
            b[i] = -b[i];
        } else if (d[i] == 0) {
            m_0++;
            Vector *vec = new Vector();
            for (int j = 0; j < n; j++) {
                vec->push_back(-A[i][j]);
            }
            A.push_back(*vec);
            b.push_back(-b[i]);
        }
    }

    // turn x_i <= 0 into -x_i >= 0
    // turn x_i with no bound into x_i1 >= 0 and x_i2 >= 0, x_i = x_i1 - x_i2
    int n_0 = n;
    std::map<int, int> extended_map;
    for (int j = 0; j < n; j++) {
        if (e[j] == -1) {
            extended_map[j] = j;
            for (int i = 0; i < m_0; i++) {
                A[i][j] = -A[i][j];
            }
            c[j] = -c[j];
        } else if (e[j] == 0) {
            extended_map[j] = n_0;
            n_0++;
            for (int i = 0; i < m_0; i++) {
                A[i].push_back(-A[i][j]);
            }
            c.push_back(-c[j]);
            x.push_back(0);
            x_dual.push_back(0);
        }
    }

    SolveResult rst = solveStandardForm(n_0, m_0, c, A, b, ans, x, ans_dual, x_dual, dual_infeasible, t_sp, t_dsp);
    
    ans = -ans;
    ans_dual = -ans_dual;

    for (auto &item : extended_map) {
        if (item.first == item.second) {
            x[item.first] *= -1;
            x_dual[item.first] *= -1;
        } else {
            x[item.first] -= x[item.second];
            x_dual[item.first] -= x_dual[item.second];
        }
    }

    return rst;
}

/*
 * Solve LP using simplex method:
 * Input: c_n, A_{m*n}, b_m
 * Output: ans, x_n 
 * Maximize c^T dot x
 * s.t. A dot x == b, x >= 0
 * using two-phase method
 */ 

SolveResult Simplex::simplexMethod(int n, int m, Vector c, Matrix A, Vector b, double &ans, Vector &x) {
    
    // Phase I
    // int k = 0;
    int pivot_i = std::min_element(b.begin(), b.end()) - b.begin(), pivot_j = -1;
    while (sgn(b[pivot_i]) < 0) {
        // std::cout << "phase 1: " << ++k << std::endl;
        auto &pivot_row = A[pivot_i];
        pivot_j = std::min_element(pivot_row.begin(), pivot_row.end()) - pivot_row.begin();
        if (sgn(pivot_row[pivot_j]) >= 0) {
            return Infeasible;
        }
        pivot(n, m, c, A, b, ans, pivot_i, pivot_j);
        pivot_i = std::min_element(b.begin(), b.end()) - b.begin();
    }

    // Phase II
    pivot_j = std::min_element(c.begin(), c.end()) - c.begin();
    while (sgn(c[pivot_j]) < 0) {
        pivot_i = -1;
        double min_ratio;
        for (int i = 0; i < m; i++) {
            if (sgn(A[i][pivot_j]) <= 0) {
                continue;
            }
            double ratio = b[i] / A[i][pivot_j];
            if (pivot_i == -1 || ratio < min_ratio) {
                pivot_i = i;
                min_ratio = ratio;
            }
        }
        if (pivot_i == -1) {
            return Unbounded;
        }
        pivot(n, m, c, A, b, ans, pivot_i, pivot_j);
        pivot_j = std::min_element(c.begin(), c.end()) - c.begin();
    }

    for (int i = 0; i < n - m; i++) {
        if (sgn(c[i]) != 0) {
            x[i] = 0;
            continue;
        }
        for (int j = 0; j < m; j++) {
            if (sgn(A[j][i] - 1) == 0) {
                x[i] = b[j];
                break;
            }
        }
    }
    return Solvable;
}

/*
 * Solve LP using dual-simplex method:
 * Input: c_n, A_{m*n}, b_m
 * Output: ans, x_n 
 * Maximize c^T dot x
 * s.t. A dot x == b, x_i >= 0
 */ 

SolveResult Simplex::dualSimplexMethod(int n, int m, Vector c, Matrix A, Vector b, double &ans, Vector &x) {
    if (sgn(*std::min_element(c.begin(), c.end())) < 0) {
        return Infeasible;
    }
    
    int pivot_i = std::min_element(b.begin(), b.end()) - b.begin(), pivot_j = -1;
    while (sgn(b[pivot_i]) < 0) {
        auto &pivot_row = A[pivot_i];
        pivot_j = -1;
        for (int j = 1; j < n; j++) {
            if (sgn(pivot_row[j]) >= 0 || sgn(c[j]) <= 0) {
                continue;
            }
            if (pivot_j == -1 || c[j] / pivot_row[j] > c[pivot_j] / pivot_row[pivot_j]) {
                pivot_j = j;
            }
        }
        if (pivot_j == -1) {
            return Infeasible;
        }
        pivot(n, m, c, A, b, ans, pivot_i, pivot_j);
        pivot_i = std::min_element(b.begin(), b.end()) - b.begin();
    }

    for (int i = 0; i < n - m; i++) {
        if (sgn(c[i]) != 0) {
            x[i] = 0;
            continue;
        }
        for (int j = 0; j < m; j++) {
            if (sgn(A[j][i] - 1) == 0) {
                x[i] = b[j];
                break;
            }
        }
    }

    return Solvable;
}