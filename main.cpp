#include "simplex.hpp"
#include <fstream>
#include <sstream>
#include <string>

int main()
{
    Simplex simplex;
    int n = 0, m = 0;
    double ans = 0, tmp;
    Vector x, c, b;
    std::vector<int> d, e;
    Matrix A;

    std::ifstream infile("F:\\in.txt");
    infile >> n >> m;

    for (int i = 0; i < n; i++)
    {
        infile >> tmp;
        c.push_back(tmp);
    }

    for (int j = 0; j < m; j++)
    {
        A.push_back(*(new Vector()));
        for (int i = 0; i < n; i++)
        {
            infile >> tmp;
            A[j].push_back(tmp);
        }
        infile >> tmp;
        b.push_back(tmp);
        infile >> tmp;
        d.push_back(tmp);
    }

    for (int i = 0; i < n; i++)
    {
        infile >> tmp;
        e.push_back(tmp);
    }

    x.resize(n);
    SolveResult k = simplex.solveNonStandardForm(n, m, c, A, b, d, e, ans, x);

    std::cout << k << std::endl;
    if (k == Solvable)
    {
        std::cout << ans << std::endl;

        for (int i = 0; i < n; i++)
        {
            std::cout << x[i] << " ";
        }
    }
    std::cout << std::endl;

    return 0;
}