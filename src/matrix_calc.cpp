#include <vector>
#include <iostream>
#include <cmath>
#include <omp.h>

#include "matrix_calc.h"

/**
 *
 * CAUTION about matrix format
 * matrix[i][j]:= value at (i)th row and (j)th column
 * 
*/

double multi_vec_vec(const std::vector<double>& a, const std::vector<double>& b, int n){
    int i;
    double res=0;
    for(i=0;i<n;i++){
        res+=a[i]*b[i];
    }
    return res;
}

void multi_mat_vec(std::vector<double>& res,
    const std::vector< std::vector<double> >& a,
    const std::vector<double>& b, int n){

    int i, j;
    res = std::vector<double>(n,0);
    #ifdef _OPENMP
    #pragma omp parallel for private(j)
    #endif
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            res[i] += a[i][j] * b[j];
        }
    }
}

void multi_mat_mat(
    std::vector< std::vector<double> >& res, 
    const std::vector< std::vector<double> >& mat_left,
    const std::vector< std::vector<double> >& mat_right){

    int i,j,k;

    if(res.size() != mat_left.size() || res[0].size() != mat_right[0].size() || mat_left[0].size() != mat_right.size()){
        std::cerr<<"matrix size error"<<std::endl;
        std::cout<<res.size()<<"="<<mat_left.size()<<std::endl;
        std::cout<<res[0].size()<<"="<<mat_right[0].size()<<std::endl;
        std::cout<<mat_left[0].size()<<"="<<mat_right.size()<<std::endl;
        return;
    }
    
    for(i=0;i<(int)(res.size());i++){
        for(j=0;j<(int)(res[0].size());j++){
            res[i][j]=0;
        }
    }
    for(i=0;i<(int)(mat_left.size());i++){
        for(j=0;j<(int)(mat_right[0].size());j++){
            for(k=0;k<(int)(mat_left[0].size());k++){
                res[i][j] += mat_left[i][k] * mat_right[k][j];
            }
        }
    }
}

void transpose_mat(
    std::vector< std::vector<double> >& trans,
    const std::vector< std::vector<double> >& mat_in){

    int i,j;
    for(i=0;i<(int)(mat_in.size());i++){
        for(j=0;j<(int)(mat_in[0].size());j++){
            trans[j][i] = mat_in[i][j];
        }
    }
}

bool inv_mat_2d(
    std::vector< std::vector<double> >& inv,
    const std::vector< std::vector<double> >& mat_in
    ){

    double det;

    det=mat_in[0][0]*mat_in[1][1]-mat_in[0][1]*mat_in[1][0];

    if(det==0){
        std::cout<<"no inverse matrix"<<std::endl;
    }else{
        //std::cout<<"ok inv"<<std::endl;
        inv[0][0] = mat_in[1][1] / det;
        inv[0][1] = mat_in[0][1] / det * (-1);
        inv[1][0] = mat_in[1][0] / det * (-1);
        inv[1][1] = mat_in[0][0] / det;
    }
    return true;
}

bool inv_mat_3d(
    std::vector< std::vector<double> >& inv,
    const std::vector< std::vector<double> >& mat_in
    ){
    double det = calc_det_3d(mat_in);
    if(det==0){
        return false;
    }else{
        double inv_det = (det==0)? 0 : 1.0 / det;

        inv[0][0] = inv_det * (mat_in[1][1] * mat_in[2][2] - mat_in[1][2] * mat_in[2][1]);
        inv[0][1] = inv_det * (mat_in[0][2] * mat_in[2][1] - mat_in[0][1] * mat_in[2][2]);
        inv[0][2] = inv_det * (mat_in[0][1] * mat_in[1][2] - mat_in[0][2] * mat_in[1][1]);
 
        inv[1][0] = inv_det * (mat_in[1][2] * mat_in[2][0] - mat_in[1][0] * mat_in[2][2]);
        inv[1][1] = inv_det * (mat_in[0][0] * mat_in[2][2] - mat_in[0][2] * mat_in[2][0]);
        inv[1][2] = inv_det * (mat_in[0][2] * mat_in[1][0] - mat_in[0][0] * mat_in[1][2]);
 
        inv[2][0] = inv_det * (mat_in[1][0] * mat_in[2][1] - mat_in[1][1] * mat_in[2][0]);
        inv[2][1] = inv_det * (mat_in[0][1] * mat_in[2][0] - mat_in[0][0] * mat_in[2][1]);
        inv[2][2] = inv_det * (mat_in[0][0] * mat_in[1][1] - mat_in[0][1] * mat_in[1][0]);
    }
    return true;
}

bool inv_mat_4d(
    std::vector< std::vector<double> >& inv,
    const std::vector< std::vector<double> >& mat_in
    ){
    int i,j;
    double det = calc_det_4d(mat_in);
    if(det==0){
        return false;
    }
    inv[0][0] = (
        mat_in[1][1] * mat_in[2][2] * mat_in[3][3] +
        mat_in[1][2] * mat_in[2][3] * mat_in[3][1] +
        mat_in[1][3] * mat_in[2][1] * mat_in[3][2] -
        mat_in[1][1] * mat_in[2][3] * mat_in[3][2] -
        mat_in[1][2] * mat_in[2][1] * mat_in[3][3] -
        mat_in[1][3] * mat_in[2][2] * mat_in[3][1]
        );

    inv[0][1] = (
        mat_in[0][1] * mat_in[2][3] * mat_in[3][2] +
        mat_in[0][2] * mat_in[2][1] * mat_in[3][3] +
        mat_in[0][3] * mat_in[2][2] * mat_in[3][1] -
        mat_in[0][1] * mat_in[2][2] * mat_in[3][3] -
        mat_in[0][2] * mat_in[2][3] * mat_in[3][1] -
        mat_in[0][3] * mat_in[2][1] * mat_in[3][2]
        );

    inv[0][2] = (
        mat_in[0][1] * mat_in[1][2] * mat_in[3][3] +
        mat_in[0][2] * mat_in[1][3] * mat_in[3][1] +
        mat_in[0][3] * mat_in[1][1] * mat_in[3][2] -
        mat_in[0][1] * mat_in[1][3] * mat_in[3][2] -
        mat_in[0][2] * mat_in[1][1] * mat_in[3][3] -
        mat_in[0][3] * mat_in[1][2] * mat_in[3][1]
        );

    inv[0][3] = (
        mat_in[0][1] * mat_in[1][3] * mat_in[2][2] +
        mat_in[0][2] * mat_in[1][1] * mat_in[2][3] +
        mat_in[0][3] * mat_in[1][2] * mat_in[2][1] -
        mat_in[0][1] * mat_in[1][2] * mat_in[2][3] -
        mat_in[0][2] * mat_in[1][3] * mat_in[2][1] -
        mat_in[0][3] * mat_in[1][1] * mat_in[2][2]
        );

    inv[1][0] = (
        mat_in[1][0] * mat_in[2][3] * mat_in[3][2] +
        mat_in[1][2] * mat_in[2][0] * mat_in[3][3] +
        mat_in[1][3] * mat_in[2][2] * mat_in[3][0] -
        mat_in[1][0] * mat_in[2][2] * mat_in[2][3] -
        mat_in[1][2] * mat_in[2][3] * mat_in[3][0] -
        mat_in[1][3] * mat_in[2][0] * mat_in[3][2]
        );

    inv[1][1] = (
        mat_in[0][0] * mat_in[2][2] * mat_in[3][3] +
        mat_in[0][2] * mat_in[2][3] * mat_in[3][0] +
        mat_in[0][3] * mat_in[2][0] * mat_in[3][2] -
        mat_in[0][0] * mat_in[2][3] * mat_in[3][2] -
        mat_in[0][2] * mat_in[2][0] * mat_in[3][3] -
        mat_in[0][3] * mat_in[2][2] * mat_in[3][0]
        );

    inv[1][2] = (
        mat_in[0][0] * mat_in[1][3] * mat_in[3][2] +
        mat_in[0][2] * mat_in[1][0] * mat_in[3][3] +
        mat_in[0][3] * mat_in[1][2] * mat_in[3][0] -
        mat_in[0][0] * mat_in[1][2] * mat_in[3][3] -
        mat_in[0][2] * mat_in[1][3] * mat_in[3][0] -
        mat_in[0][3] * mat_in[1][0] * mat_in[3][2]
        );

    inv[1][3] = (
        mat_in[0][0] * mat_in[1][2] * mat_in[2][3] +
        mat_in[0][2] * mat_in[1][3] * mat_in[2][0] +
        mat_in[0][3] * mat_in[1][0] * mat_in[2][2] -
        mat_in[0][0] * mat_in[1][3] * mat_in[2][2] -
        mat_in[0][2] * mat_in[1][0] * mat_in[2][3] -
        mat_in[0][3] * mat_in[1][2] * mat_in[2][0]
        );

    inv[2][0] = (
        mat_in[1][0] * mat_in[2][1] * mat_in[3][3] +
        mat_in[1][1] * mat_in[2][3] * mat_in[3][0] +
        mat_in[1][3] * mat_in[2][0] * mat_in[3][1] -
        mat_in[1][0] * mat_in[2][3] * mat_in[3][1] -
        mat_in[1][1] * mat_in[2][0] * mat_in[3][3] -
        mat_in[1][3] * mat_in[2][1] * mat_in[3][0]
        );

    inv[2][1] = (
        mat_in[0][0] * mat_in[2][3] * mat_in[3][1] +
        mat_in[0][1] * mat_in[2][0] * mat_in[3][3] +
        mat_in[0][3] * mat_in[2][1] * mat_in[3][0] -
        mat_in[0][0] * mat_in[2][1] * mat_in[3][3] -
        mat_in[0][1] * mat_in[2][3] * mat_in[3][0] -
        mat_in[0][3] * mat_in[2][0] * mat_in[3][1]
        );

    inv[2][2] = (
        mat_in[0][0] * mat_in[1][1] * mat_in[3][3] +
        mat_in[0][1] * mat_in[1][3] * mat_in[3][0] +
        mat_in[0][3] * mat_in[1][0] * mat_in[3][1] -
        mat_in[0][0] * mat_in[1][3] * mat_in[3][1] -
        mat_in[0][1] * mat_in[1][0] * mat_in[3][3] -
        mat_in[0][3] * mat_in[1][1] * mat_in[3][0]
        );

    inv[2][3] = (
        mat_in[0][0] * mat_in[1][3] * mat_in[2][1] +
        mat_in[0][1] * mat_in[1][0] * mat_in[2][3] +
        mat_in[0][3] * mat_in[1][1] * mat_in[2][0] -
        mat_in[0][0] * mat_in[1][1] * mat_in[2][3] -
        mat_in[0][1] * mat_in[1][3] * mat_in[2][0] -
        mat_in[0][3] * mat_in[1][0] * mat_in[2][1]
        );

    inv[3][0] = (
        mat_in[1][0] * mat_in[2][2] * mat_in[3][1] +
        mat_in[1][1] * mat_in[2][0] * mat_in[3][2] +
        mat_in[1][2] * mat_in[2][1] * mat_in[3][0] -
        mat_in[1][0] * mat_in[2][1] * mat_in[3][2] -
        mat_in[1][1] * mat_in[2][2] * mat_in[3][0] -
        mat_in[1][2] * mat_in[2][0] * mat_in[3][1]
        );

    inv[3][1] = (
        mat_in[0][0] * mat_in[2][1] * mat_in[3][2] +
        mat_in[0][1] * mat_in[2][2] * mat_in[3][0] +
        mat_in[0][2] * mat_in[2][0] * mat_in[3][1] -
        mat_in[0][0] * mat_in[2][2] * mat_in[3][1] -
        mat_in[0][1] * mat_in[2][0] * mat_in[3][2] -
        mat_in[0][2] * mat_in[2][1] * mat_in[3][0]
        );

    inv[3][2] = (
        mat_in[0][0] * mat_in[1][2] * mat_in[3][1] +
        mat_in[0][1] * mat_in[1][0] * mat_in[3][2] +
        mat_in[0][2] * mat_in[1][1] * mat_in[3][0] -
        mat_in[0][0] * mat_in[1][1] * mat_in[3][2] -
        mat_in[0][1] * mat_in[1][2] * mat_in[3][0] -
        mat_in[0][2] * mat_in[1][0] * mat_in[3][1]
        );

    inv[3][3] = (
        mat_in[0][0] * mat_in[1][1] * mat_in[2][2] +
        mat_in[0][1] * mat_in[1][2] * mat_in[2][0] +
        mat_in[0][2] * mat_in[1][0] * mat_in[2][1] -
        mat_in[0][0] * mat_in[1][2] * mat_in[2][1] -
        mat_in[0][1] * mat_in[1][0] * mat_in[2][2] -
        mat_in[0][2] * mat_in[1][1] * mat_in[2][0]
        );

    double inv_det = 1.0 / det;
    for(i=0;i<4;i++){
        for(j=0;j<4;j++){
            inv[i][j] *= inv_det;
        }
    }
    return true;
}

double calc_trace(const std::vector< std::vector<double> >& mat_in, int n){
    int i;
    double trace=0;
    for(i=0;i<n;i++){
        trace+=mat_in[i][i];
    }
    return trace;
}

double calc_det_2d(const std::vector< std::vector<double> >& mat_in){
    return mat_in[0][0] * mat_in[1][1] - mat_in[0][1] * mat_in[1][0];
}

double calc_det_3d(const std::vector< std::vector<double> >& mat_in){
    return mat_in[0][0] * mat_in[1][1] * mat_in[2][2]    //   a11a22a33
         + mat_in[1][0] * mat_in[2][1] * mat_in[0][2]    // + a21a32a13
         + mat_in[2][0] * mat_in[0][1] * mat_in[1][2]    // + a31a12a23
         - mat_in[0][0] * mat_in[2][1] * mat_in[1][2]    // - a11a32a23
         - mat_in[2][0] * mat_in[1][1] * mat_in[0][2]    // - a31a22a13
         - mat_in[1][0] * mat_in[0][1] * mat_in[2][2];   // - a21a12a33
}

/*
double calc_det_4d(const std::vector< std::vector<double> >& mat_in){
    double val1, val2, val3, val4;
    std::vector<std::vector<double>> mat = std::vector< std::vector<double> >(3,std::vector<double>(3,0));

    mat = {
        {mat_in[1][1], mat_in[1][2], mat_in[1][3]},
        {mat_in[2][1], mat_in[2][2], mat_in[2][3]},
        {mat_in[3][1], mat_in[3][2], mat_in[3][3]}
    };
    val1 = mat_in[0][0] * calc_det_3d(mat);

    mat = {
        {mat_in[0][1], mat_in[0][2], mat_in[0][3]},
        {mat_in[2][1], mat_in[2][2], mat_in[2][3]},
        {mat_in[3][1], mat_in[3][2], mat_in[3][3]}
    };
    val2 = mat_in[1][0] * calc_det_3d(mat);

    mat = {
        {mat_in[0][1], mat_in[0][2], mat_in[0][3]},
        {mat_in[1][1], mat_in[1][2], mat_in[1][3]},
        {mat_in[3][1], mat_in[3][2], mat_in[3][3]}
    };
    val3 = mat_in[2][0] * calc_det_3d(mat);

    mat = {
        {mat_in[0][1], mat_in[0][2], mat_in[0][3]},
        {mat_in[1][1], mat_in[1][2], mat_in[1][3]},
        {mat_in[2][1], mat_in[2][2], mat_in[2][3]}
    };
    val4 = mat_in[3][0] * calc_det_3d(mat);

    return val1 - val2 + val3 - val4;
}
//*/
// This is faster than the function above (commented out) when n >= 4
// PROBLEM: sometimes fail to calculate and output nan value
//*
double calc_det_4d(const std::vector< std::vector<double> >& mat_in){
    std::vector<std::vector<double>> mat = mat_in;
    int x, y, i, n = 4;
    double det = 1, r;
    // convert to upper triangle matrix, and multiply by diagonal element
    for(y=0;y<n-1;y++){
        if(std::abs(mat[y][y]) < 1e-6){
            // if diagonal element is zero switch with nonzero row
            for(i=y+1;i<n;i++){
                if(mat[i][y] != 0){
                    break;
                }
            }
            if(i<n){
                for(x=0;x<n;x++){
                    std::swap(mat[i][x], mat[y][x]);
                }
                // reverse the sign because column was switched
                det = - det;
            }
        }
        for(i=y+1;i<n;i++){
            r = mat[i][y] / mat[y][y];
            for(x=y;x<n;x++){
                mat[i][x] -= r * mat[y][x];
            }
        }
        det *= mat[y][y];
    }
    det *= mat[y][y];
    return det;
}
//*/