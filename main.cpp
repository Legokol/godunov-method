#include <Eigen/Dense>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>

using Eigen::Vector3d;

const double kappa = 1.4;

struct Point {
    double x;
    double t;
};

struct State {
    Vector3d x;
    double p;
};

struct Solution {
    Point point;
    State state;
};

// Преобразование между переменными U и W
Vector3d getW(const State &state) {
    return Vector3d(state.x(0), state.x(1) / state.x(0), state.p);
}

State getState(const Vector3d &w) {
    return {Vector3d(w(0), w(0) * w(1), w(2) / (kappa - 1)), w(2)};
}

// Метод Ньютона для решения уравнений
double newtonMethod(double x0, const std::function<double(double)> &f, const std::function<double(double)> &df, const double relativeError = 1e-6) {
    double x = x0 - f(x0) / df(x0);
    double error = abs((x - x0) / 0.5 / (x + x0));
    while (error > relativeError) {
        x0 = x;
        x = x0 - f(x0) / df(x0);
        error = abs((x - x0) / 0.5 / (x + x0));
    }
    return x;
}

// Вспомогательные функции
double a(double rho) {
    return 2 / (kappa + 1) / rho;
}

double b(double p) {
    return p * (kappa - 1) / (kappa + 1);
}


// Функции для определения давления на разрыве
double funcP(double p, const Vector3d &w) {
    if (p >= w(2)) {
        return (p - w(2)) * sqrt(a(w(0)) / (p + b(w(2))));
    }
    double c = sqrt(kappa * w(2) / w(0));
    return 2 * c / (kappa - 1) * (pow(p / w(2), (kappa - 1) / 2 / kappa) - 1);
}

double dfuncP(double p, const Vector3d &w) {
    if (p > w(2)) {
        return sqrt(a(w(0)) / (p + b(w(2)))) * (1 - (p - w(2)) / 2 / (p + b(w(2))));
    }
    double c = sqrt(kappa * w(2) / w(0));
    return 2 * c / (kappa - 1) * (pow(p / w(2), -(kappa - 1) / 2 / kappa) - 1);
}


// Функция расчёта потока
Vector3d calcW(const Vector3d &wLeft, const Vector3d &wRight) {
    if ((wLeft - wRight).cwiseAbs().maxCoeff() < 1e-6) {
            return wLeft;
        }
    double p0 = 0.5 * (wLeft(2) + wRight(2));

    auto f = [&, wLeft, wRight](double p) { return funcP(p, wLeft) + funcP(p, wRight) + wRight(1) - wLeft(1); };
    auto df = [&, wLeft, wRight](double p) { return dfuncP(p, wLeft) + dfuncP(p, wRight); };

    if (p0 - f(p0) / df(p0) < 0) {
        double cLeft = sqrt(kappa * wLeft(2) / wLeft(0));
        double cRight = sqrt(kappa * wRight(2) / wRight(0));
        p0 = pow((cLeft + cRight - 0.5 * (kappa - 1) * (wRight(1) - wLeft(1))) / (cLeft / pow(wLeft(2), (kappa - 1) / 2 / kappa) + cRight / pow(wRight(2), (kappa - 1) / 2 / kappa)), 2 * kappa / (kappa - 1));
    }

    double p = newtonMethod(p0, f, df);
    double u = 0.5 * (wLeft(1) + wRight(1)) + 0.5 * (funcP(p, wRight) - funcP(p, wLeft));

    if (u >= 0) {
        double cLeft = sqrt(kappa * wLeft(2) / wLeft(0));
        if (p > wLeft(2)) {
            double shockVelocity = wLeft(1) - cLeft * sqrt((kappa + 1) / 2 / kappa * p / wLeft(2) + (kappa - 1) / 2 / kappa);
            if (shockVelocity > 0) {
                return wLeft;
            } else {
                double rho = wLeft(0) * (p / wLeft(2) + (kappa - 1) / (kappa + 1)) / ((kappa - 1) / (kappa + 1) * p / wLeft(2) + 1);
                return Vector3d(rho, u, p);
            }
        } else {
            double headVelocity = wLeft(1) - cLeft;
            double tailVelocity = u - cLeft * pow(p / wLeft(2), (kappa - 1) / 2 / kappa);
            if (headVelocity > 0) {
                return wLeft;
            } else if (tailVelocity < 0) {
                double rho = wLeft(0) * pow(p / wLeft(2), 1 / kappa);
                return Vector3d(rho, u, p);
            } else {
                double rho = wLeft(0) * pow(2 / (kappa + 1) + (kappa - 1) / (kappa + 1) / cLeft * wLeft(1), 2 / (kappa - 1));
                double p = wLeft(2) * pow(2 / (kappa + 1) + (kappa - 1) / (kappa + 1) / cLeft * wLeft(1), 2 * kappa / (kappa - 1));
                double u = 2 / (kappa + 1) * (cLeft + (kappa - 1) / 2 * wLeft(1));
                return Vector3d(rho, u, p);
            }
        }
    } else {
        double cRight = sqrt(kappa * wRight(2) / wRight(0));
        if (p > wRight(2)) {
            double shockVelocity = wRight(1) + cRight * sqrt((kappa + 1) / 2 / kappa * p / wRight(2) + (kappa - 1) / 2 / kappa);
            if (shockVelocity < 0) {
                return wRight;
            } else {
                double rho = wRight(0) * (p / wRight(2) + (kappa - 1) / (kappa + 1)) / ((kappa - 1) / (kappa + 1) * p / wRight(2) + 1);
                return Vector3d(rho, u, p);
            }
        } else {
            double headVelocity = wRight(1) + cRight;
            double tailVelocity = u + cRight * pow(p / wRight(2), (kappa - 1) / 2 / kappa);
            if (headVelocity < 0) {
                return wRight;
            } else if (tailVelocity > 0) {
                double rho = wRight(0) * pow(p / wRight(2), 1 / kappa);
                return Vector3d(rho, u, p);
            } else {
                double rho = wRight(0) * pow(2 / (kappa + 1) - (kappa - 1) / (kappa + 1) / cRight * wRight(1), 2 / (kappa - 1));
                double p = wRight(2) * pow(2 / (kappa + 1) - (kappa - 1) / (kappa + 1) / cRight * wRight(1), 2 * kappa / (kappa - 1));
                double u = 2 / (kappa + 1) * (-cRight + (kappa - 1) / 2 * wRight(1));
                return Vector3d(rho, u, p);
            }
        }
    }
    return Vector3d();
}

Vector3d fluxVector(const State &state) {
    return Vector3d(state.x(1), state.x(1) * state.x(1) / state.x(0) + state.p, state.x(1) / state.x(0) * (state.x(2) + state.p));
}


int main() {
    // Параметры задачи
    double L = 10;
    double T = 2e-2;

    // Начальные условия
    double rhoLeft = 13;
    double vLeft = 0;
    double pLeft = 10e5;

    double rhoRight = 1.3;
    double vRight = 0;
    double pRight = 1e5;

    // Шаги по пространству и времени
    double h = 0.2;
    double tau = 1e-5;

    // Тестирование
    // Vector3d wLeft{1, 0, 1};
    // Vector3d wRight{0.125, 0, 0.1};

    std::vector<std::vector<Solution>> solution(1);
    solution[0].resize(2 * L / h);
    for (int j = 0; j < solution[0].size(); ++j) {
        solution[0][j].point = {-L + (0.5 + j) * h, 0};
        if (solution[0][j].point.x < 0) {
            solution[0][j].state = {{rhoLeft, rhoLeft * vLeft, pLeft / (kappa - 1)}, pLeft};
        } else {
            solution[0][j].state = {{rhoRight, rhoRight * vRight, pRight / (kappa - 1)}, pRight};
        }
    }

    double t = tau;
    std::vector<Solution> timeLayer(solution[0].size());
    while (t <= T) {
        Vector3d wCurrent = getW(solution.back()[0].state);
        Vector3d fluxLeft = fluxVector(getState(wCurrent));
        for (int j = 0; j < timeLayer.size() - 1; ++j) {
            // std::cout << j << std::endl;
            Vector3d wRight = getW(solution.back()[j + 1].state);
            Vector3d fluxRight = fluxVector(getState(calcW(wCurrent, wRight)));

            timeLayer[j].point = {-L + (0.5 + j) * h, t};
            timeLayer[j].state.x = solution.back()[j].state.x + tau / h * (fluxLeft - fluxRight);
            timeLayer[j].state.p = timeLayer[j].state.x(2) * (kappa - 1);

            fluxLeft = fluxRight;
            wCurrent = wRight;
        }

        Vector3d fluxRight = fluxVector(getState(wCurrent));
        timeLayer.back().point = {L - 0.5 * h, t};
        timeLayer.back().state.x = solution.back().back().state.x + tau / h * (fluxLeft - fluxRight);
        timeLayer.back().state.p = timeLayer.back().state.x(2) * (kappa - 1);

        solution.push_back(timeLayer);
        t += tau;
        // for (int i = 0; i < 100; i++)
        // {
        //     if (abs(i * 2e-3 - t) < 1e-6) {
        //         std::cout << t << std::endl;
        //     }
        // }
        // std::cout << t << std::endl;
    }

    // Запись результатов расчёта
    // std::vector<std::ofstream> writer(4);
    // writer[0].open("results/p.txt");
    // writer[1].open("results/rho.txt");
    // writer[2].open("results/u.txt");
    // writer[3].open("results/e.txt");
    // for (int i = 0; i < solution.size(); ++i) {
    //     writer[0] << solution[i][0].point.t << ' ';
    //     writer[1] << solution[i][0].point.t << ' ';
    //     writer[2] << solution[i][0].point.t << ' ';
    //     writer[3] << solution[i][0].point.t << ' ';
    //     for (int j = 0; j < solution[i].size(); ++j) {
    //         writer[0] << solution[i][j].state.p << ' ';
    //         writer[1] << solution[i][j].state.x(0) << ' ';
    //         writer[2] << solution[i][j].state.x(1) / solution[i][j].state.x(0) << ' ';
    //         writer[3] << solution[i][j].state.x(2) / solution[i][j].state.x(0) << ' ';
    //     }
    //     writer[0] << '\n';
    //     writer[1] << '\n';
    //     writer[2] << '\n';
    //     writer[3] << '\n';
    // }
    // writer[0].close();
    // writer[1].close();
    // writer[2].close();
    // writer[3].close();
    return 0;
}