#include <regex>
#include <iostream>
#include <fstream>
#include "eigen/Eigen/Sparse"
#include "eigen/Eigen/Dense"
#include "MobMatrix.hpp"
#include "eigen/Eigen/SparseLU"
#include <random>
#include <time.h>

typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SpMat;
typedef Eigen::Triplet<double> Triplet;

// method for calculating the Moore-Penrose pseudo-Inverse as recommended by Eigen developers
template<typename _Matrix_Type_>
_Matrix_Type_ pseudoInverse(const _Matrix_Type_ &a, double epsilon = 0.000001)
{

    Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeFullU | Eigen::ComputeFullV);
        // For a non-square matrix
        // Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeThinU | Eigen::ComputeThinV);
    double tolerance = epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);
    return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
}

double contarPesos(const SpMat& Matrix);
SpMat::InnerIterator chooseLink(const SpMat& Matrix, double p);

#define LEER_MATRIZ
//#define ESCRIBIR_MATRIZ

int main(){

    const std::string state = "tx";
    const std::string citPat = "citiesMult/"+ state +"/Citypatch.txt";
    const std::string mobNet = "citiesMult/"+ state +"/mobnetwork.txt";
    const std::string popAr = "citiesMult/"+ state +"/Poparea.txt";

    MobMatrix T{citPat, mobNet, popAr};

    SpMat Matrix(T.N, T.N);
    SpMat Laplaciano(T.N, T.N);
    SpMat Resistances(T.N, T.N);
    std::vector<Triplet> tripletes1, tripletes2;

#ifdef ESCRIBIR_MATRIZ
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

    Resistances = Matrix;
    for(int i = 0; i < Matrix.outerSize(); i++){
        SpMat::InnerIterator itRes(Resistances,i);
        for(SpMat::InnerIterator it(Matrix, i); it; ++it){
            itRes.valueRef() = (identidad.col(it.row()) - identidad.col(it.col())).transpose() * x * (identidad.col(it.row()) - identidad.col(it.col()));
            ++itRes;
        }
        Resistances.coeffRef(i,i) = 1;
    }

    //Reescribir la matriz Matrix, esta vez identica a MobMatrix
    Matrix.setZero(); Matrix.data().squeeze();
    tripletes1.clear();
    for(int i = 0; i < T.N; i++){ //Construccion de la matriz
        for(int j = 0; j < T.vecinos[i]; j++){
            tripletes1.emplace_back(i, T.Mvecinos[i][j], T.Mpesos[i][j]/T.population[i]); //Pesos relativos
        }
    }
    
    Matrix.setFromTriplets(tripletes1.begin(), tripletes1.end());
    Matrix.makeCompressed();

    Resistances = Resistances.cwiseProduct(Matrix);
    //Guardar la matriz de resistencias
    std::ofstream file("citiesMult/"+ state +"/results.txt");
    for(int i = 0; i < T.N; i++){
        for(int j = 0; j < T.N; j++){
            file << Resistances.coeff(i,j) << "\t";
        }
        file << "\n";
    }
    file.flush(); file.close();
#endif
    //Leer la matriz de resistencias
#ifdef LEER_MATRIZ
    Matrix.setZero(); Matrix.data().squeeze();
    tripletes1.clear();
    for(int i = 0; i < T.N; i++){ //Construccion de la matriz
        for(int j = 0; j < T.vecinos[i]; j++){
            if(T.population[i] != 0)
                tripletes1.emplace_back(i, T.Mvecinos[i][j], T.Mpesos[i][j]/T.population[i]); //Pesos relativos
            else tripletes1.emplace_back(i, T.Mvecinos[i][j], 0);
        }
    }
    Matrix.setFromTriplets(tripletes1.begin(), tripletes1.end());
    Matrix.makeCompressed();

    std::ifstream fileread("citiesMult/"+ state +"/results.txt");

    tripletes1.clear();
    double a;
    for(int i = 0; i < T.N; i++){
        for(int j = 0; j < T.N; j++){
            fileread >> a;
            tripletes1.emplace_back(i,j,a);
        }
    }
    Resistances.setFromTriplets(tripletes1.begin(),tripletes1.end());
    fileread.close();
    //std::cout << Resistances << std::endl;

#endif
    //Generador numeros aleatorios
    static std::default_random_engine generator;
    static std::uniform_real_distribution<double> distribution(0.0, contarPesos(Resistances)); //Hasta que llegue a todos los pesos de la red
    generator.seed(static_cast<unsigned int>(time(NULL)));

    tripletes1.clear();
    const int LINKS = T.Links;
    double p;
    for(int s = 0; s < LINKS; s++){
        p = distribution(generator);
        auto it = chooseLink(Resistances, p);
        tripletes1.emplace_back(it.row(), it.col(), Matrix.coeff(it.row(),it.col())); //Meter el peso
    }
    Matrix.setZero(); Matrix.data().squeeze(); //Rehacer la matriz de cero
    Matrix.setFromTriplets(tripletes1.begin(), tripletes1.end());

    //Renormalizar los pesos a la poblacion total

    for(int i = 0; i < Matrix.outerSize(); i++){
        double norma = 0;
        for(SpMat::InnerIterator it(Matrix, i); it; ++it){
            norma += it.value();
        }
        for(SpMat::InnerIterator it(Matrix, i); it; ++it){
            it.valueRef() *= T.population[it.row()] / norma; //Todos los pesos ahora suman la poblacion
        }
    }

    std::ofstream file2{"citiesMult/"+state+"/newmobnetwork.txt"};
    for(int i = 0; i < Matrix.outerSize(); i++){
        for(SpMat::InnerIterator it(Matrix, i); it; ++it){
            file2 << it.row() << " " << it.col() << " " << static_cast<int>(it.value()) << std::endl;
        }
    }
    file2.close();

    return 0;
}

double contarPesos(const SpMat& Matrix){
    double peso = 0;
    for(int i = 0; i < Matrix.outerSize(); i++){
        for(SpMat::InnerIterator it(Matrix, i); it; ++it){
            peso += it.value();
        }
    }
    return peso;
}

SpMat::InnerIterator chooseLink(const SpMat& Matrix, double p){
    double temp = 0;
    for(int i = 0; i < Matrix.outerSize(); i++){
        for(SpMat::InnerIterator it(Matrix, i); it; ++it){
            temp += it.value();
            if(temp > p) return it;
        }
    }
    exit(1);
}