#include <regex>
#include <iostream>
#include <fstream>
#include "eigen/Eigen/Sparse"
#include "eigen/Eigen/Dense"
#include "MobMatrix.hpp"
#include "eigen/Eigen/SparseLU"

typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SpMat;
typedef Eigen::Triplet<double> Triplet;

// method for calculating the Moore-Penrose pseudo-Inverse as recommended by Eigen developers
template<typename _Matrix_Type_>
    _Matrix_Type_ pseudoInverse(const _Matrix_Type_ &a, double epsilon = 0.00001)
    {

        Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeFullU | Eigen::ComputeFullV);
            // For a non-square matrix
            // Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeThinU | Eigen::ComputeThinV);
        double tolerance = epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);
        return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
    }

int main(){

    const std::string state = "tx";
    const std::string citPat = "citiesMult/"+ state +"/Citypatch.txt";
    const std::string mobNet = "citiesMult/"+ state +"/mobnetwork.txt";
    const std::string popAr = "citiesMult/"+ state +"/Poparea.txt";

    MobMatrix T{citPat, mobNet, popAr};

    SpMat Matrix(T.N, T.N);
    SpMat Laplaciano(T.N, T.N);
    std::vector<Triplet> tripletes1, tripletes2;


    for(int i = 0; i < T.N; i++){ //Construccion de la matriz
        int laplDiag = 0;

        for(int j = 0; j < T.vecinos[i]; j++){
            double value = T.Mpesos[i][j];
            if(i != T.Mvecinos[i][j]){
                tripletes1.emplace_back(i, T.Mvecinos[i][j], value/2);//mitad en i,j
                tripletes1.emplace_back(T.Mvecinos[i][j], i, value/2);//mitad en j,i (simetrica)
            }
        }
    }

    Matrix.setFromTriplets(tripletes1.begin(), tripletes1.end());
    Matrix.makeCompressed();

    for(int i = 0; i < Matrix.outerSize(); i++){
        int vecinos = 0;
        for(SpMat::InnerIterator it(Matrix, i); it; ++it){
            tripletes2.emplace_back(it.row(), it.col(), -1);
            vecinos++;
        }
        tripletes2.emplace_back(i, i, vecinos);
    }

    
    Laplaciano.setFromTriplets(tripletes2.begin(), tripletes2.end());
    Laplaciano.makeCompressed();

    auto LaplacianoDense = Laplaciano.toDense();
    auto x = pseudoInverse(LaplacianoDense);

    
    SpMat identidad(T.N,T.N); identidad.setIdentity();

    SpMat Resistances = Matrix;
    for(int i = 0; i < Matrix.outerSize(); i++){
        SpMat::InnerIterator itRes(Resistances,i);
        for(SpMat::InnerIterator it(Matrix, i); it; ++it){
            itRes.valueRef() = (identidad.col(it.row()) - identidad.col(it.col())).transpose() * x * (identidad.col(it.row()) - identidad.col(it.col()));
            ++itRes;
        }
    }
    
    std::ofstream file("citiesMult/"+ state +"/results.txt");
    file << Resistances;
    file.close();

    return 0;
}