#pragma once

#include <map>
#include <vector>
#include <iostream>
#include <algorithm>

typedef std::vector<double> Vector;
typedef std::vector<std::vector<double>> Matrix;

enum SolveResult {
    Infeasible  = -1,
    Solvable    = 0,
    Unbounded   = 1
};

class Simplex
{
private:
    constexpr static double eps = 1e-9;
    int sgn(double x) {
        return x < -eps ? -1 : (x > eps ? 1 : 0);
    }
    void pivot(int n, int m, Vector &c, Matrix &A, Vector &d, double &ans, int pivot_i, int pivot_j);
public:
    Simplex();
    ~Simplex();
    SolveResult solveStandardForm(int n, int m, Vector &c, Matrix &A, Vector &d, double &ans, Vector &x);
    SolveResult solveNonStandardForm(int n, int m, Vector &c, Matrix &A, Vector &d, std::vector<int> &b, std::vector<int> &e, double &ans, Vector &x);
};

Simplex::Simplex()
{
}

Simplex::~Simplex()
{
}

/*
 * Solve linear programming in standard form:
 * Input: c_n, A_{m*n}, d_m
 * Output: ans, x_n 
 * Maximize c^T dot x
 * s.t. A dot x \leq d
 * it will add slack variables
 */
SolveResult Simplex::solveStandardForm(int n, int m, Vector &c, Matrix &A, Vector &d, double &ans, Vector &x) {
    if (n != c.size() || m != A.size() || m != d.size()) {
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

    // Phase I
    int pivot_i = std::min_element(d.begin(), d.end()) - d.begin(), pivot_j = -1;
    while (sgn(d[pivot_i]) < 0) {
        auto &pivot_row = A[pivot_i];
        pivot_j = std::min_element(pivot_row.begin(), pivot_row.end()) - pivot_row.begin();
        if (sgn(pivot_row[pivot_j]) >= 0) {
            return Infeasible;
        }
        pivot(n_0, m, c, A, d, ans, pivot_i, pivot_j);
        pivot_i = std::min_element(d.begin(), d.end()) - d.begin();
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
            double ratio = d[i] / A[i][pivot_j];
            if (pivot_i == -1 || ratio < min_ratio) {
                pivot_i = i;
                min_ratio = ratio;
            }
        }
        if (pivot_i == -1) {
            return Unbounded;
        }
        pivot(n_0, m, c, A, d, ans, pivot_i, pivot_j);
        pivot_j = std::min_element(c.begin(), c.end()) - c.begin();
    }

    for (int i = 0; i < n; i++) {
        if (sgn(c[i]) != 0) {
            x[i] = 0;
            continue;
        }
        for (int j = 0; j < m; j++) {
            if (sgn(A[j][i] - 1) == 0) {
                x[i] = d[j];
                break;
            }
        }
    }

    // std::cout << "A:" << std::endl;
    // for (int i = 0; i < m; i++) {
    //     for (int j = 0; j < n_0; j++) {
    //         std::cout << A[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "c:" << std::endl;
    // for (int j = 0; j < n_0; j++) {
    //     std::cout << c[j] << " ";
    // }
    // std::cout << std::endl;

    // std::cout << "d:" << std::endl;
    // for (int i = 0; i < m; i++) {
    //     std::cout << d[i] << " ";
    // }
    // std::cout << std::endl;

    // std::cout << "ans: " << ans << std::endl;
    
    return Solvable;
}


/*
 *  Pivoting operation:
 * 
 * 
 * 
 */

void Simplex::pivot(int n, int m, Vector &c, Matrix &A, Vector &d, double &ans, int pivot_i, int pivot_j) {
    auto &A_p = A[pivot_i];
    
    double factor = A_p[pivot_j];
    for (int j = 0; j < n; j++) {
        A_p[j] /= factor;
    }
    d[pivot_i] /= factor;
    
    for (int i = 0; i < m; i++) {
        if (i == pivot_i || sgn(A[i][pivot_j]) == 0) {
            continue;
        }
        factor = A[i][pivot_j];
        for (int j = 0; j < n; j++) {
            A[i][j] -= factor * A_p[j];
        }
        d[i] -= factor * d[pivot_i];
    }

    factor = c[pivot_j];
    for (int j = 0; j < n; j++) {
        c[j] -= factor * A_p[j];
    }
    ans -= factor * d[pivot_i];
}

/*
 * Solve linear programming in non-standard form:
 * Input: c_n, A_{m*n}, d_m, b_m, e_n
 * Output: ans, x_n 
 * Minimize c^T dot x
 * s.t. A dot x \leq or \eq or \geq d
 * it will add slack variables
 */

SolveResult Simplex::solveNonStandardForm(int n, int m, Vector &c, Matrix &A, Vector &d, std::vector<int> &b, std::vector<int> &e, double &ans, Vector &x) {
    for (int i = 0; i < n; i++) {
        c[i] = -c[i];
    }

    int m_0 = m;
    for (int i = 0; i < m; i++) {
        if (b[i] == -1) {
            for (int j = 0; j < n; j++) {
                A[i][j] = -A[i][j];
            }
            d[i] = -d[i];
        } else if (b[i] == 0) {
            m_0++;
            Vector *vec = new Vector();
            for (int j = 0; j < n; j++) {
                vec->push_back(-A[i][j]);
            }
            A.push_back(*vec);
            d.push_back(-d[i]);
        }
    }

    int n_0 = n;
    std::map<int, int> extended_map;
    for (int j = 0; j < n; j++) {
        if (e[j] == -1) {
            for (int i = 0; i < m; i++) {
                A[i][j] = -A[i][j];
            }
            c[j] = -c[j];
        } else if (e[j] == 0) {
            extended_map[j] = n_0;
            n_0++;
            for (int i = 0; i < m; i++) {
                A[i].push_back(-A[i][j]);
            }
            c.push_back(-c[j]);
            x.push_back(0);
        }
    }

    SolveResult rst = solveStandardForm(n_0, m_0, c, A, d, ans, x);
    for (auto &item : extended_map) {
        x[item.first] -= x[item.second];
    }
    return rst;
}