#include<cmath>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>
#include<iostream>

int main(){

    // Basic Example of cpp
    std::cout << "Example of cpp \n";
    float a = 1.0, b = 2.0;
    std::cout << a << std::endl;
    std::cout << a/b << std::endl;
    std::cout << std::sqrt(b) << std::endl;
    std::cout << std::acos(-1) << std::endl;
    std::cout << std::sin(30.0/180.0*acos(-1)) << std::endl;

    // Example of vector
    std::cout << "Example of vector \n";
    // vector definition
    Eigen::Vector3f v(1.0f,2.0f,3.0f);
    Eigen::Vector3f w(1.0f,0.0f,0.0f);
    // vector output
    std::cout << "Example of output \n";
    std::cout << v << std::endl;
    // vector add
    std::cout << "Example of add \n";
    std::cout << v + w << std::endl;
    // vector scalar multiply
    std::cout << "Example of scalar multiply \n";
    std::cout << v * 3.0f << std::endl;
    std::cout << 2.0f * v << std::endl;

    // Example of matrix
    std::cout << "Example of matrix \n";
    // matrix definition
    Eigen::Matrix3f i,j;
    i << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
    j << 2.0, 3.0, 1.0, 4.0, 6.0, 5.0, 9.0, 7.0, 8.0;
    // matrix output
    std::cout << "Example of output \n";
    std::cout << i << std::endl;
    // matrix add i + j
    // matrix scalar multiply i * 2.0
    // matrix multiply i * j
    // matrix multiply vector i * v
    std::cout<< "P(2,1) \n";
    Eigen::Matrix3f rotate;
    Eigen::Vector3f P(2.0f,1.0f,1.0f);
    float sin = std::sin((45.0/180.0 * std::acos(-1)));
    float cos = std::cos((45.0/180.0 * std::acos(-1)));
    rotate << cos, -sin, 1, 
              sin, cos,  2, 
              0,   0,    1;
    std::cout<< rotate * P << std::endl;

    std::cout<<"Po"<<std::endl;
    Eigen::Vector3f origin(2,1,1);
    Eigen::Vector3f point(0, 0, 0);
    Eigen::Matrix3f trans;
    trans << std::cos(45),std::cos(-45),1,std::cos(-45),std::cos(45),2,0,0,1;
    point = trans * origin;
    std::cout<<point<<std::endl;
    
    Eigen::Vector3f point_P(2.0f, 1.0f, 1.0f);
    Eigen::Matrix3f T_rot, T_shift;
    
    float theta = 45.0f / 180.0f * std::acos(-1);
    T_rot <<
        std::cos(theta),    -std::sin(theta),   0,
        std::sin(theta),    std::cos(theta),    0,
        0,                  0,                  1;
    T_shift <<
        1,                  0,                  1,
        0,                  1,                  2,
        0,                  0,                  1;
    
    std::cout << theta << std::endl;
    std::cout << T_rot << std::endl;
    std::cout << T_shift << std::endl;

    Eigen::Vector3f point_P_transformed;
    point_P_transformed = T_shift * T_rot * point_P;
    std::cout << point_P_transformed << std::endl;

    return 0;
}
