#include "simplex.hpp"

int main() {
    Simplex simplex;
    int n = 2, m = 3;
    double ans = 0;
    Vector x = {0, 0}, c = {-2, -3}, d = {10, -12, -12};
    Matrix A = {{1, 1}, {-1, -2}, {-2, -1}};
    simplex.solveStandardForm(n, m, c, A, d, ans, x);
    std::cout << ans << std::endl;
    for (int i = 0; i < n; i++) {
        std::cout << x[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}