#include "simplex.hpp"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

int main()
{
    //Initialization
    Simplex simplex;
    int n = 0, m = 0;
    double ans = 0, ans_dual, tmp;
    Vector x, x_dual, c, b;
    std::vector<int> d, e;
    Matrix A;

    //Get the input from stdin
    std::cin >> n >> m;

    for (int i = 0; i < n; i++)
    {
        std::cin >> tmp;
        c.push_back(tmp);
    }

    for (int j = 0; j < m; j++)
    {
        A.push_back(*(new Vector()));
        for (int i = 0; i < n; i++)
        {
            std::cin >> tmp;
            A[j].push_back(tmp);
        }
        std::cin >> tmp;
        b.push_back(tmp);
        std::cin >> tmp;
        d.push_back(tmp);
    }

    for (int i = 0; i < n; i++)
    {
        std::cin >> tmp;
        e.push_back(tmp);
    }
    
    x.resize(n);
    x_dual.resize(n);
    bool dual_infeasible = false;
    //Call solveNonStandardForm()
    //k == -1: Infeasible
    //k == 0: Unbounded
    //k == 1: Solvable
    SolveResult k = simplex.solveNonStandardForm(n, m, c, A, b, d, e, ans, x, ans_dual, x_dual, dual_infeasible);

    //Output the result of Simplex and Dual simplex
    std::cout << k << std::endl;
    if (k == Solvable)
    {
        std::cout << "Simplex: " << std::endl;
        std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(6) << ans << std::endl;

        for (int i = 0; i < n; i++)
        {
            std::cout << x[i] << " ";
        }

        if (!dual_infeasible)
        {
            std::cout << std::endl
                      << "Dual simplex: " << std::endl;
            std::cout << ans_dual << std::endl;

            for (int i = 0; i < n; i++)
            {
                std::cout << x_dual[i] << " ";
            }
        }
        else
            std::cout << std::endl
                      << "Dual simplex infeasible" << std::endl;
    }
    std::cout << std::endl;

    return 0;
}