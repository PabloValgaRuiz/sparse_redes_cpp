#include <regex>
#include "eigen/Eigen/Sparse"
#include "MobMatrix.hpp"
#include "eigen/Eigen/SparseLU"

typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;
typedef Eigen::Triplet<double> Triplet;


int main(){

    const std::string state = "tx";
    const std::string citPat = "citiesMult/"+ state +"/Citypatch.txt";
    const std::string mobNet = "citiesMult/"+ state +"/mobnetwork.txt";
    const std::string popAr = "citiesMult/"+ state +"/Poparea.txt";

    MobMatrix T{citPat, mobNet, popAr};
    std::vector<size_t> vecinos;
    for(int i = 0; i < T.N; i++)
        vecinos.push_back(T.vecinos[i]);

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
            else
                tripletes1.emplace_back(i, T.Mvecinos[i][j], value); //diagonal
        }
    }
    for(int i = 0; i < T.N; i++){ //Construccion del laplaciano
        int laplDiag = 0;
        for(int j = 0; j < T.vecinos[i]; j++){
            if(i == T.Mvecinos[i][j]){
                laplDiag = 1; //Ver si hay un autoloop
                tripletes2.emplace_back(i,i, T.vecinos[i] - 1);
            }
            else if((i > T.Mvecinos[i][j]) && (laplDiag == 0)){
                tripletes2.emplace_back(i,i, T.vecinos[i]);
                laplDiag = 1;
            }
            else
                tripletes2.emplace_back(i, T.Mvecinos[i][j], 1);
        }
    }

    Matrix.setFromTriplets(tripletes1.begin(), tripletes1.end());
    Matrix.makeCompressed();
    Laplaciano.setFromTriplets(tripletes2.begin(), tripletes2.end());
    Laplaciano.makeCompressed();

    return 0;
}