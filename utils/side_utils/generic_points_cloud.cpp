#include "utils.h"

#include <fstream>
#include <vector>

#include <math.h>
#include <math_constants.h>

int main(int argc, char *argv[])
{
    int N = 10;
    int M = 10;
    int R = 1;

    double points[N][M][3];

    double r = 0.5;
    double fi = 0;
    double th = 0;

    std::ofstream ofile("size.msz");

    ofile << R * M * N << '\n';
    for (size_t k = 0; k < R; k++)
    {
        r -= (k * 0.01);

        for (size_t i = 0; i < N; i++)
        {
            fi += (2 * M_PI * i / N);
            for (size_t j = 0; j < M; j++)
            {
                th += (2 * M_PI * j / M);
                points[i][j][0] = r * sin(th) * cos(fi);
                points[i][j][1] = r * sin(th) * sin(fi);
                points[i][j][2] = r * cos(th);

                ofile << points[i][j][0] << ' '
                      << points[i][j][1] << ' '
                      << points[i][j][2] << ' ' << 0.05 << '\n';
            }
        }
    }
    ofile << 0;

    ofile.close();

    return 0;
}