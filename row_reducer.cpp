#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

std::vector<double> gaussian(std::vector<std::vector<double>> &A)
{
    int numRows = A.size();
    int numCols = A[0].size();
    const double EPS = 0;

    for (int j = 0; j < numCols - 1; ++j)
    {
        for (int i = j + 1; i < numRows; ++i)
        {
            if (std::abs(A[j][j]) < EPS)
            {
                throw std::runtime_error("Zero pivot encountered");
            }
            double M = (A[i][j] != 0) ? -(A[i][j] / A[j][j]) : 0;
            for (int k = j; k < numCols; ++k)
            {
                A[i][k] += M * A[j][k];
            }
        }
    }

    std::vector<double> sol(numCols - 1, 0);
    for (int i = numRows - 1; i >= 0; --i)
    {
        double total = 0;
        for (int j = numCols - 2; j >= 0; --j)
        {
            total += sol[j] * A[i][j];
        }
        double b = A[i][numCols - 1];
        if (std::abs(A[i][i]) < EPS)
        {
            throw std::runtime_error("Zero pivot encountered");
        }
        sol[i] = (b - total) / A[i][i];
    }

    return sol;
}

int main()
{
    try
    {
        std::vector<std::vector<double>> A = {{1, 2, -1, 2}, {0, 3, 1, 4}, {2, -1, 1, 2}};
        std::vector<double> resultA = gaussian(A);
        std::cout << "Solution for A: ";
        for (double val : resultA)
        {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
    catch (const std::runtime_error &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}