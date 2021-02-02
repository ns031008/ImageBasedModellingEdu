//
// Created by caoqi on 2018/8/31.
//

//3D:  1.36939, -1.17123, 7.04869
//obs: 0.180123 -0.156584

#include "sfm/bundle_adjustment.h"
/*
 * This function computes the Jacobian entries for the given camera and
 * 3D point pair that leads to one observation.
 *
 * The camera block 'cam_x_ptr' and 'cam_y_ptr' is:
 * - ID 0: Derivative of focal length f
 * - ID 1-2: Derivative of distortion parameters k0, k1
 * - ID 3-5: Derivative of translation t0, t1, t2
 * - ID 6-8: Derivative of rotation w0, w1, w2
 *
 * The 3D point block 'point_x_ptr' and 'point_y_ptr' is:
 * - ID 0-2: Derivative in x, y, and z direction.
 *
 * The function that leads to the observation is given as follows:
 *
 *   u = f * D(x,y) * x  (image observation x coordinate)
 *   v = f * D(x,y) * y  (image observation y coordinate)
 *
 * with the following definitions:
 *
 *   xc = R0 * X + t0  (homogeneous projection)
 *   yc = R1 * X + t1  (homogeneous projection)
 *   zc = R2 * X + t2  (homogeneous projection)
 *   x = xc / zc  (central projection)
 *   y = yc / zc  (central projection)
 *   D(x, y) = 1 + k0 (x^2 + y^2) + k1 (x^2 + y^2)^2  (distortion)
 */

/**
  * /description 给定一个相机参数和一个三维点坐标，求解雅各比矩阵，即公式中的df(theta)/dtheta
  * @param cam       相机参数
  * @param point     三维点坐标
  * @param cam_x_ptr 重投影坐标x 相对于相机参数的偏导数，相机有9个参数： [0] 焦距f; [1-2] 径向畸变系数k1, k2; [3-5] 平移向量 t1, t2, t3
  *                                                               [6-8] 旋转矩阵（角轴向量）
  * @param cam_y_ptr    重投影坐标y 相对于相机参数的偏导数，相机有9个参数
  * @param point_x_ptr  重投影坐标x 相对于三维点坐标的偏导数
  * @param point_y_ptr  重投影坐标y 相对于三维点坐标的偏导数
  */
void jacobian(sfm::ba::Camera const &cam,
              sfm::ba::Point3D const &point,
              double *cam_x_ptr, double *cam_y_ptr,
              double *point_x_ptr, double *point_y_ptr)
{
    const double f = cam.focal_length;
    const double *R = cam.rotation;
    const double *t = cam.translation;
    const double *X = point.pos;
    const double k0 = cam.distortion[0];
    const double k1 = cam.distortion[1];

    // 相机焦距的偏导数
    double X_c[3], xy[2], uv[2];
    X_c[0] = R[0] * X[0] + R[1] * X[1] + R[2] * X[2] + t[0];
    X_c[1] = R[3] * X[0] + R[4] * X[1] + R[5] * X[2] + t[1];
    X_c[2] = R[6] * X[0] + R[7] * X[1] + R[8] * X[2] + t[2];
    xy[0] = X_c[0] / X_c[2];
    xy[1] = X_c[1] / X_c[2];
    double r_2 = xy[0] * xy[0] + xy[1] * xy[1];
    double d = 1 + (k0 + k1 * r_2) * r_2;
    uv[0] = f * d * xy[0];
    uv[1] = f * d * xy[1];

    cam_x_ptr[0] = d * xy[0];
    cam_y_ptr[0] = d * xy[1];

    double deriv_u_d = f * xy[0];
    double deriv_v_d = f * xy[1];
    double deriv_d_k0 = r_2;
    double deriv_d_k1 = r_2 * r_2;

    // 相机径向畸变的偏导数
    cam_x_ptr[1] = f * xy[0] * r_2;
    cam_x_ptr[2] = f * xy[0] * r_2 * r_2;
    cam_y_ptr[1] = f * xy[1] * r_2;
    cam_y_ptr[2] = f * xy[1] * r_2 * r_2;

    double deriv_x_xc = 1 / X_c[2];
    double deriv_x_yc = 0;
    double deriv_x_zc = -xy[0] / X_c[2];

    double deriv_y_xc = 0;
    double deriv_y_yc = 1 / X_c[2];
    double deriv_y_zc = -xy[1] / X_c[2];

    double deriv_u_x = f * d;
    double deriv_v_y = f * d;
    double deriv_d_r2 = k0 + 2 * k1 * r_2;
    double deriv_r2_xc = 2 * xy[0] / X_c[2];
    double deriv_r2_yc = 2 * xy[1] / X_c[2];
    double deriv_r2_zc = -2 * r_2 / X_c[2];
    double deriv_d_xc = deriv_d_r2 * deriv_r2_xc;
    double deriv_d_yc = deriv_d_r2 * deriv_r2_yc;
    double deriv_d_zc = deriv_d_r2 * deriv_r2_zc;

    double deriv_u_xc = deriv_u_d * deriv_d_xc + deriv_u_x * deriv_x_xc;
    double deriv_u_yc = deriv_u_d * deriv_d_yc + deriv_u_x * deriv_x_yc;
    double deriv_u_zc = deriv_u_d * deriv_d_zc + deriv_u_x * deriv_x_zc;
    double deriv_v_xc = deriv_v_d * deriv_d_xc + deriv_v_y * deriv_y_xc;
    double deriv_v_yc = deriv_v_d * deriv_d_yc + deriv_v_y * deriv_y_yc;
    double deriv_v_zc = deriv_v_d * deriv_d_zc + deriv_v_y * deriv_y_zc;

    // 相机平移向量的偏导数
    cam_x_ptr[3] = deriv_u_xc;
    cam_x_ptr[4] = deriv_u_yc;
    cam_x_ptr[5] = deriv_u_zc;
    cam_y_ptr[3] = deriv_v_xc;
    cam_y_ptr[4] = deriv_v_yc;
    cam_y_ptr[5] = deriv_v_zc;

    double r0X = X_c[0] - t[0];
    double r1X = X_c[1] - t[1];
    double r2X = X_c[2] - t[2];
    double deriv_xc_w0 = 0;
    double deriv_xc_w1 = r2X;
    double deriv_xc_w2 = -r1X;
    double deriv_yc_w0 = -r2X;
    double deriv_yc_w1 = 0;
    double deriv_yc_w2 = r0X;
    double deriv_zc_w0 = r1X;
    double deriv_zc_w1 = -r0X;
    double deriv_zc_w2 = 0;

    // 相机旋转矩阵的偏导数
    cam_x_ptr[6] = deriv_u_yc * deriv_yc_w0 + deriv_u_zc * deriv_zc_w0;
    cam_x_ptr[7] = deriv_u_xc * deriv_xc_w1 + deriv_u_zc * deriv_zc_w1;
    cam_x_ptr[8] = deriv_u_xc * deriv_xc_w2 + deriv_u_yc * deriv_yc_w2;
    cam_y_ptr[6] = deriv_v_yc * deriv_yc_w0 + deriv_v_zc * deriv_zc_w0;
    cam_y_ptr[7] = deriv_v_xc * deriv_xc_w1 + deriv_v_zc * deriv_zc_w1;
    cam_y_ptr[8] = deriv_v_xc * deriv_xc_w2 + deriv_v_yc * deriv_yc_w2;

    double deriv_xc_X = R[0];
    double deriv_xc_Y = R[1];
    double deriv_xc_Z = R[2];
    double deriv_yc_X = R[3];
    double deriv_yc_Y = R[4];
    double deriv_yc_Z = R[5];
    double deriv_zc_X = R[6];
    double deriv_zc_Y = R[7];
    double deriv_zc_Z = R[8];

    // 三维点的偏导数
    point_x_ptr[0] = deriv_u_xc * deriv_xc_X + deriv_u_yc * deriv_yc_X + deriv_u_zc * deriv_zc_X;
    point_x_ptr[1] = deriv_u_xc * deriv_xc_Y + deriv_u_yc * deriv_yc_Y + deriv_u_zc * deriv_zc_Y;
    point_x_ptr[2] = deriv_u_xc * deriv_xc_Z + deriv_u_yc * deriv_yc_Z + deriv_u_zc * deriv_zc_Z;
    point_y_ptr[0] = deriv_v_xc * deriv_xc_X + deriv_v_yc * deriv_yc_X + deriv_v_zc * deriv_zc_X;
    point_y_ptr[1] = deriv_v_xc * deriv_xc_Y + deriv_v_yc * deriv_yc_Y + deriv_v_zc * deriv_zc_Y;
    point_y_ptr[2] = deriv_v_xc * deriv_xc_Z + deriv_v_yc * deriv_yc_Z + deriv_v_zc * deriv_zc_Z;
}
int main(int argc, char *argv[])
{

    sfm::ba::Camera cam;
    cam.focal_length = 0.919654;
    cam.distortion[0] = -0.108298;
    cam.distortion[1] = 0.103775;

    cam.rotation[0] = 0.999999;
    cam.rotation[1] = -0.000676196;
    cam.rotation[2] = -0.0013484;
    cam.rotation[3] = 0.000663243;
    cam.rotation[4] = 0.999949;
    cam.rotation[5] = -0.0104095;
    cam.rotation[6] = 0.00135482;
    cam.rotation[7] = 0.0104087;
    cam.rotation[8] = 0.999949;

    cam.translation[0] = 0.00278292;
    cam.translation[1] = 0.0587996;
    cam.translation[2] = -0.127624;

    sfm::ba::Point3D pt3D;
    pt3D.pos[0] = 1.36939;
    pt3D.pos[1] = -1.17123;
    pt3D.pos[2] = 7.04869;

    double cam_x_ptr[9] = {0};
    double cam_y_ptr[9] = {0};
    double point_x_ptr[3] = {0};
    double point_y_ptr[3] = {0};

    jacobian(cam, pt3D, cam_x_ptr, cam_y_ptr, point_x_ptr, point_y_ptr);

    std::cout << "Result is :" << std::endl;
    std::cout << "cam_x_ptr: ";
    for (int i = 0; i < 9; i++)
    {
        std::cout << cam_x_ptr[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "cam_y_ptr: ";
    for (int i = 0; i < 9; i++)
    {

        std::cout << cam_y_ptr[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "point_x_ptr: ";
    std::cout << point_x_ptr[0] << " " << point_x_ptr[1] << " " << point_x_ptr[2] << std::endl;

    std::cout << "point_y_ptr: ";
    std::cout << point_y_ptr[0] << " " << point_y_ptr[1] << " " << point_y_ptr[2] << std::endl;

    std::cout << "\nResult should be :\n"
              << "cam_x_ptr: 0.195942 0.0123983 0.000847141 0.131188 0.000847456 -0.0257388 0.0260453 0.95832 0.164303\n"
              << "cam_y_ptr: -0.170272 -0.010774 -0.000736159 0.000847456 0.131426 0.0223669 -0.952795 -0.0244697 0.179883\n"
              << "point_x_ptr: 0.131153 0.000490796 -0.0259232\n"
              << "point_y_ptr: 0.000964926 0.131652 0.0209965\n";

    return 0;
}
