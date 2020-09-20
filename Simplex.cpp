#include<iostream>
#include<vector>

/*
Reference: ISBN	978-4-320-01616-3, p37
*/

class Simplex{
    const int NOTFOUND = -1;
    const double INIT = 1e9;
    int n, m;
    std::vector<std::vector<double>> dict;
    std::vector<int> basis, nonbasis;
    std::vector<double> parameter;
    int searchNonbasis();
    int searchBasis(int nonbasisID);
    void pivot(int b, int nb);
public:
    Simplex(){};
    Simplex(int n, int m, const std::vector<std::vector<double>> &A);
    bool run();
    double getParameter(int i);
};

Simplex::Simplex(int N, int M, const std::vector<std::vector<double>> &A){
    n = N;
    m = M;
    dict = A;
    parameter.assign(n + m, 0.0);
    for(int i = 0; i < n; i++){
        nonbasis.push_back(i);
    }
    for(int i = 0; i < m; i++){
        basis.push_back(n + i);
    }
}

int Simplex::searchNonbasis(){
    for(int i = 0; i < nonbasis.size(); i++){
        if(dict[0][i + 1] < 0){
            return i + 1;
        }
    }
    return NOTFOUND;
}

int Simplex::searchBasis(int nonbasisID){
    int r = -1;
    double ratio = INIT;
    for(int i = 1; i < 1 + basis.size(); i++){
        if(dict[i][nonbasisID] < 0){
            double tmp_ratio = -dict[i][0] / dict[i][nonbasisID];
            if(tmp_ratio < ratio){
                ratio = tmp_ratio;
                r = i;
            }
        }
    }
    return r;
}

void Simplex::pivot(int bID, int nbID){
    double ratio = 0.0;
    for(int i = 0; i < 1 + basis.size(); i++){
        if(bID != i){
            ratio = dict[i][nbID] / dict[bID][nbID];
            for(int j = 0; j < 1 + nonbasis.size(); j++){
                dict[i][j] -= dict[bID][j] * ratio;
            }
        }
    }
    // for bID row
    ratio = -1.0 / dict[bID][nbID];
    for(int j = 0; j < 1 + nonbasis.size(); j++){
        if(j == nbID){
            dict[bID][j] = -1.0 * ratio;
        }else{
            dict[bID][j] *= ratio;
        }
    }
    // swap bID and nbID
    std::swap(basis[bID - 1], nonbasis[nbID - 1]);
}

bool Simplex::run(){
    int basisID = NOTFOUND;
    int nonBasisID = searchNonbasis();
    // run simplex algorithm
    while(nonBasisID != NOTFOUND){
        basisID = searchBasis(nonBasisID);
        if(basisID == NOTFOUND){
            return false;
        }
        pivot(basisID, nonBasisID);
        nonBasisID = searchNonbasis();
    }

    // calculate parameter
    for(auto nb: nonbasis){
        parameter[nb] = 0.0;
    }
    for(int i = 0; i < basis.size(); i++){
        int b = basis[i];
        parameter[b] = dict[i + 1][0];
        for(int j = 1; j < 1 + nonbasis.size(); j++){
            parameter[b] += dict[i][j] * parameter[nonbasis[j - 1]];
        }
    }
    return true;
}

double Simplex::getParameter(int i){
    return parameter[i];
}

int main(){
    int n = 3;
    int m = 3;
    double array[4][4] = {
        {0, -2, -1, -1},
        {4, -2, -2, 1},
        {4, -2, 0, -4},
        {1, 4, -3, 1}
    };
    std::vector<std::vector<double>> A(n + 1, std::vector<double>(m + 1));
    for(int i = 0; i < n + 1; i++){
        for(int j = 0; j < m + 1; j++){
            A[i][j] = array[i][j];
        }
    }
    Simplex simplex(3, 3, A);
    simplex.run();
    for(int i = 0; i < n + m; i++){
        std::cout << simplex.getParameter(i) << " ";
    }
    return 0;
}
