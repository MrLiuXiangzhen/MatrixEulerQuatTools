//
// Created by LiuXiangzhen on 22-9-21.
//
#include "iostream"
#include "MyTF/Quaternion.hpp"
#include "MyTF/MatrixRot.hpp"

using namespace std;

double PI = 3.14;

int main() {
    MyTF::Quaternion q1;
    double roll, pitch, yaw;

    q1.setRPY(PI/2, 0, PI/2);
    cout << "euler to quat:  "  << endl << q1.x() << " " << q1.y() << " " << q1.z() << " " << q1.w() << endl;

    MyTF::Matrix3x3 mat;
    mat.setRotation(q1);
    cout << "quat to Matrix: " << endl
        << mat.getRow(0)[0] << " " << mat.getRow(0)[1] << " " << mat.getRow(0)[2] << endl
        << mat.getRow(1)[0] << " " << mat.getRow(1)[1] << " " << mat.getRow(1)[2] << endl
        << mat.getRow(2)[0] << " " << mat.getRow(2)[1] << " " << mat.getRow(2)[2] << endl;

    MyTF::Matrix3x3(mat.getRow(0)[0], mat.getRow(0)[1], mat.getRow(0)[2],
                    mat.getRow(1)[0], mat.getRow(1)[1], mat.getRow(1)[2],
                    mat.getRow(2)[0], mat.getRow(2)[1], mat.getRow(2)[2]).getRPY(roll, pitch, yaw);
    cout << "Matrix to RPY: " << endl << roll << " " << pitch << " " << yaw << endl;

    return 0;
}
