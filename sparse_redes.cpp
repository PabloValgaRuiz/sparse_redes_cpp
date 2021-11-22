#include <regex>
#include "eigen/Eigen/SparseCore"
#include "MobMatrix.hpp"


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
    std::vector<Triplet> tripletes;

    Matrix.reserve(vecinos);

    for(int i = 0; i < T.N; i++){
        for(int j = 0; j < T.vecinos[i]; j++){
            Matrix.insert(i, T.Mvecinos[i][j]) = T.Mpesos[i][j];
        }
    }
    Matrix.makeCompressed();


    return 0;
}