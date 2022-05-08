#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

constexpr unsigned int iterations{50};
constexpr double TOLX{1e-10};
constexpr double TOLF{1e-10};

void jacobiMethod(vector<vector<double>> &matrixA, vector<double> &vectorB, vector<double> vectorX);

void gaussSeidelMethod(vector<vector<double>> &matrixA, vector<double> &vectorB, vector<double> vectorX);

void SORMethod(vector<vector<double>> &matrixA, vector<double> &vectorB, vector<double> vectorX, double &&omega);

int main()
{
    vector<vector<double>> matrixA{{100, -1,  2,   -3},
                                   {1,   200, -4,  5},
                                   {-2,  4,   300, -6},
                                   {3,   -5,  6,   400}};

    vector<double> vectorB{116, -226, 912, -1174};
    vector<double> vectorX{2, 2, 2, 2};

    jacobiMethod(matrixA, vectorB, vectorX);
    gaussSeidelMethod(matrixA, vectorB, vectorX);
    SORMethod(matrixA, vectorB, vectorX, 0.5);
}

void jacobiMethod(vector<vector<double>> &matrixA, vector<double> &vectorB, vector<double> vectorX)
{
    double errorEstimator, residue, errorEstimatorTmp, residueTmp, tmp;
    vector<double> vectorX2{0, 0, 0, 0};

    cout << "Jacobi's method" << endl;
    cout << setw(4) << "No.";
    cout << "|" << setw(12) << "VectorX[0]";
    cout << "|" << setw(12) << "VectorX[1]";
    cout << "|" << setw(12) << "VectorX[2]";
    cout << "|" << setw(12) << "VectorX[3]";
    cout << "|" << setw(20) << "Error estimator";
    cout << "|" << setw(12) << "Residue" << endl;

    for (int i = 0; i < iterations; ++i)
    {
        errorEstimator = 0;
        residue = 0;
        for (int j = 0; j < 4; ++j)
        {
            tmp = 0;
            for (int k = 0; k < 4; ++k)
                if (j != k)
                    tmp += matrixA[j][k] * vectorX[k];

            vectorX2[j] = (vectorB[j] - tmp) / matrixA[j][j];

            tmp = 0;
            for (int k = 0; k < 4; ++k)
                tmp += matrixA[j][k] * vectorX[k];


            errorEstimatorTmp = fabs(vectorX2[j] - vectorX[j]);
            if (errorEstimatorTmp > errorEstimator)
                errorEstimator = errorEstimatorTmp;

            residueTmp = fabs(tmp - vectorB[j]);
            if (residueTmp > residue)
                residue = residueTmp;

            vectorX[j] = vectorX2[j];
        }

        cout << setw(4) << i + 1;
        cout << "|" << setw(12) << vectorX[0];
        cout << "|" << setw(12) << vectorX[1];
        cout << "|" << setw(12) << vectorX[2];
        cout << "|" << setw(12) << vectorX[3];
        cout << "|" << setw(20) << errorEstimator;
        cout << "|" << setw(12) << residue << endl;

        if (errorEstimator < TOLX && residue < TOLF)
        {
            cout << "The calculations are accurate" << endl;
            break;
        }
    }
    cout << endl;
}

void gaussSeidelMethod(vector<vector<double>> &matrixA, vector<double> &vectorB, vector<double> vectorX)
{
    double errorEstimator, residue, errorEstimatorTmp, residueTmp, tmp;
    vector<double> vectorX2{0, 0, 0, 0};

    cout << "Gauss-Seidel method" << endl;
    cout << setw(4) << "No.";
    cout << "|" << setw(12) << "VectorX[0]";
    cout << "|" << setw(12) << "VectorX[1]";
    cout << "|" << setw(12) << "VectorX[2]";
    cout << "|" << setw(12) << "VectorX[3]";
    cout << "|" << setw(20) << "Error estimator";
    cout << "|" << setw(12) << "Residue" << endl;

    for (int i = 0; i < iterations; ++i)
    {
        errorEstimator = 0;
        residue = 0;
        for (int j = 0; j < 4; ++j)
        {
            tmp = 0;
            for (int k = 0; k < 4; ++k)
                if (j != k)
                    tmp += matrixA[j][k] * vectorX[k];

            vectorX2[j] = vectorX[j];
            vectorX[j] = (vectorB[j] - tmp) / matrixA[j][j];

            errorEstimatorTmp = fabs(vectorX2[j] - vectorX[j]);
            if (errorEstimatorTmp > errorEstimator)
                errorEstimator = errorEstimatorTmp;
        }
        for (int j = 0; j < 4; ++j)
        {
            tmp = 0;
            for (int k = 0; k < 4; ++k)
                tmp += matrixA[j][k] * vectorX[k];

            residueTmp = fabs(tmp - vectorB[j]);
            if (residueTmp > residue)
                residue = residueTmp;
        }

        cout << setw(4) << i + 1;
        cout << "|" << setw(12) << vectorX[0];
        cout << "|" << setw(12) << vectorX[1];
        cout << "|" << setw(12) << vectorX[2];
        cout << "|" << setw(12) << vectorX[3];
        cout << "|" << setw(20) << errorEstimator;
        cout << "|" << setw(12) << residue << endl;

        if (errorEstimator < TOLX && residue < TOLF)
        {
            cout << "The calculations are accurate" << endl;
            break;
        }
    }
    cout << endl;
}

void SORMethod(vector<vector<double>> &matrixA, vector<double> &vectorB, vector<double> vectorX, double &&omega)
{
    double errorEstimator, residue, errorEstimatorTmp, residueTmp, tmp;
    vector<double> vectorX2{0, 0, 0, 0};

    cout << "SOR method" << endl;
    cout << setw(4) << "No.";
    cout << "|" << setw(12) << "VectorX[0]";
    cout << "|" << setw(12) << "VectorX[1]";
    cout << "|" << setw(12) << "VectorX[2]";
    cout << "|" << setw(12) << "VectorX[3]";
    cout << "|" << setw(20) << "Error estimator";
    cout << "|" << setw(12) << "Residue" << endl;

    for (int i = 0; i < iterations; ++i)
    {
        errorEstimator = 0;
        residue = 0;
        for (int j = 0; j < 4; ++j)
        {
            tmp = 0;
            for (int k = 0; k < 4; ++k)
                if (j != k)
                    tmp += matrixA[j][k] * vectorX[k];

            vectorX2[j] = vectorX[j];
            vectorX[j] = (1 - omega) * vectorX[j] + (omega * (vectorB[j] - tmp) / matrixA[j][j]);

            tmp = 0;
            for (int k = 0; k < 4; ++k)
                tmp += matrixA[j][k] * vectorX[k];


            errorEstimatorTmp = fabs(vectorX2[j] - vectorX[j]);
            if (errorEstimatorTmp > errorEstimator)
                errorEstimator = errorEstimatorTmp;

            residueTmp = fabs(tmp - vectorB[j]);
            if (residueTmp > residue)
                residue = residueTmp;
        }

        cout << setw(4) << i + 1;
        cout << "|" << setw(12) << vectorX[0];
        cout << "|" << setw(12) << vectorX[1];
        cout << "|" << setw(12) << vectorX[2];
        cout << "|" << setw(12) << vectorX[3];
        cout << "|" << setw(20) << errorEstimator;
        cout << "|" << setw(12) << residue << endl;

        if (errorEstimator < TOLX && residue < TOLF)
        {
            cout << "The calculations are accurate" << endl;
            break;
        }
    }
    cout << endl;
}