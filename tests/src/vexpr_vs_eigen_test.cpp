#include <iostream>
#include <iomanip>
#include <cassert>
#include <chrono>

#include "Eigen/Dense"
#include "vexpr/vexpr.h"

// Функция для измерения времени выполнения
template <typename Func>
double measure_time(Func f, int iterations = 1000)
{
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i)
    {
        f();
    }
    auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(end - start).count() / iterations;
}

// Тест корректности
template <size_t N>
void test_correctness()
{
    // Инициализация
    Eigen::VectorXd eigen_v1(N), eigen_v2(N);
    vexpr::VectorX<double, N> my_v1, my_v2;

    for (size_t i = 0; i < N; ++i)
    {
        eigen_v1[i] = my_v1[i] = i + 1;
        eigen_v2[i] = my_v2[i] = N - i;
    }

    // Тестирование операций
    auto eigen_add = eigen_v1 + eigen_v2;
    auto my_add = my_v1 + my_v2;

    for (size_t i = 0; i < N; ++i)
    {
        assert(std::abs(eigen_add[i] - my_add[i]) < 1e-9);
    }

    auto eigen_sub = eigen_v1 - eigen_v2;
    auto my_sub = my_v1 - my_v2;

    for (size_t i = 0; i < N; ++i)
    {
        assert(std::abs(eigen_sub[i] - my_sub[i]) < 1e-9);
    }

    double eigen_dot = eigen_v1.dot(eigen_v2);
    double my_dot = my_v1.dot(my_v2);
    assert(std::abs(eigen_dot - my_dot) < 1e-9);

    double scalar = 2.5;
    auto eigen_mul = eigen_v1 * scalar;
    auto my_mul = my_v1 * scalar;

    for (size_t i = 0; i < N; ++i)
    {
        assert(std::abs(eigen_mul[i] - my_mul[i]) < 1e-9);
    }

    double eigen_norm = eigen_v1.norm();
    double my_norm = my_v1.norm();
    assert(std::abs(eigen_norm - my_norm) < 1e-9);

    std::cout << "All correctness tests passed for N = " << N << std::endl;
}

// Тест производительности
template <size_t N>
void test_performance()
{
    // Инициализация
    Eigen::VectorXd eigen_v1 = Eigen::VectorXd::Random(N);
    Eigen::VectorXd eigen_v2 = Eigen::VectorXd::Random(N);

    vexpr::VectorX<double, N> my_v1, my_v2;
    for (size_t i = 0; i < N; ++i)
    {
        my_v1[i] = eigen_v1[i];
        my_v2[i] = eigen_v2[i];
    }

    // Измерение времени
    double scalar = 2.5;

    double eigen_add_time = measure_time([&]()
                                         { volatile auto result = eigen_v1 + eigen_v2; });

    double my_add_time = measure_time([&]()
                                      { volatile auto result = my_v1 + my_v2; });

    double eigen_dot_time = measure_time([&]()
                                         { volatile auto result = eigen_v1.dot(eigen_v2); });

    double my_dot_time = measure_time([&]()
                                      { volatile auto result = my_v1.dot(my_v2); });

    double eigen_mul_time = measure_time([&]()
                                         { volatile auto result = eigen_v1 * scalar; });

    double my_mul_time = measure_time([&]()
                                      { volatile auto result = my_v1 * scalar; });

    double eigen_div_time = measure_time([&]()
                                         { volatile auto result = eigen_v1 / scalar; });

    double my_div_time = measure_time([&]()
                                      { volatile auto result = my_v1 / scalar; });

    double eigen_norm_time = measure_time([&]()
                                          { volatile auto result = eigen_v1.norm(); });

    double my_norm_time = measure_time([&]()
                                       { volatile auto result = my_v1.norm(); });

    // Настройки формата
    std::cout << "\nPerformance results for N = " << N << ":\n";
    std::cout << std::string(65, '-') << '\n';
    std::cout << std::fixed << std::setprecision(9);
    const int w_op = 15;
    const int w_time = 18;
    const int w_speedup = 15;

    std::cout << std::left
              << std::setw(w_op) << "Operation"
              << std::right
              << std::setw(w_time) << "Eigen3 (s)"
              << std::setw(w_time) << "MyVector (s)"
              << std::setw(w_speedup) << "Speedup"
              << "\n"
              << std::string(65, '-') << '\n';

    auto print_row = [&](const std::string &op, double e, double m)
    {
        std::cout << std::left << std::setw(w_op) << op;
        std::cout << std::right << std::fixed
                  << std::setw(w_time) << e
                  << std::setw(w_time) << m
                  << std::setw(w_speedup) << std::setprecision(2) << (e / m)
                  << "\n";
        std::cout << std::setprecision(9); // Восстанавливаем точность
    };

    print_row("Addition", eigen_add_time, my_add_time);
    print_row("Dot Product", eigen_dot_time, my_dot_time);
    print_row("Scalar Mult", eigen_mul_time, my_mul_time);
    print_row("Div", eigen_div_time, my_div_time);
    print_row("norm", eigen_norm_time, my_norm_time);
}

int main()
{
    // Тестируем для разных размерностей
    test_correctness<3>();
    test_correctness<10>();
    test_correctness<100>();

    test_performance<3>();
    test_performance<10>();
    test_performance<100>();
    test_performance<1000>();

    return 0;
}
