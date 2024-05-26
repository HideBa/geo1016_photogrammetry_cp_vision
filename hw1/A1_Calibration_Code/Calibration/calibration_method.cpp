/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "calibration.h"
#include "matrix_algo.h"
#include <cassert>

using namespace easy3d;

/**
 * TODO: Finish this function for calibrating a camera from the corresponding
 * 3D-2D point pairs. You may define a few functions for some sub-tasks.
 * @return True on success, otherwise false. On success, the camera parameters
 * are returned by fx, fy, cx, cy, skew, R, and t).
 */
bool Calibration::calibration(
    const std::vector<Vector3D> &points_3d, /// input: An array of 3D points.
    const std::vector<Vector2D>
        &points_2d, /// input: An array of 2D image points.
    double &fx,     /// output: focal length (i.e., K[0][0]).
    double &fy,     /// output: focal length (i.e., K[1][1]).
    double &cx,  /// output: x component of the principal point (i.e., K[0][2]).
    double &cy,  /// output: y component of the principal point (i.e., K[1][2]).
    double &s,   /// output: skew factor (i.e., K[0][1]), which is s = -alpha *
                 /// cot(theta).
    Matrix33 &R, /// output: the 3x3 rotation matrix encoding camera rotation.
    Vector3D &t) /// output：a 3D vector encoding camera translation.
{
  std::cout << "\nTODO: I am going to implement the calibration() function in "
               "the following file:\n"
               "\t    - calibration_method.cpp\n\n";

  std::cout << "[Liangliang]:\n"
               "\tCamera calibration requires computing the SVD and inverse of "
               "matrices.\n"
               "\tIn this assignment, I provide you with a 'Matrix' and a "
               "'Vector' data structures for storing and\n"
               "\tmanipulating matrices and vectors of arbitrary sizes. I also "
               "wrote some code to show you how to:\n"
               "\t    - compute the SVD of a matrix;\n"
               "\t    - compute the inverse of a matrix;\n"
               "\t    - compute the transpose of a matrix.\n\n"
               "\tFeel free to use any of the provided data structures and "
               "functions. The commonly used linear algebra\n"
               "\tfunctions are provided in the following files:\n"
               "\t    - Calibration/matrix.h  Matrices of arbitrary dimensions "
               "and related functions.\n"
               "\t    - Calibration/vector.h  Vectors of arbitrary dimensions "
               "and related functions.\n"
               "\t    - Calibration/matrix_algo.h  Determinant, inverse, SVD, "
               "linear least-squares...\n"
               "\tPlease refer to the above files for a complete list of "
               "useful functions and their usage.\n\n"
               "\tIn your final submission, please\n"
               "\t    - delete ALL unrelated test or debug code and avoid "
               "unnecessary output.\n"
               "\t    - include all the source code (and please do NOT modify "
               "the structure of the directories).\n"
               "\t    - do NOT include the 'build' directory (which contains "
               "the intermediate files in a build step).\n"
               "\t    - make sure your code compiles and can reproduce your "
               "results without ANY modification.\n\n"
            << std::flush;

  std::cout
      << "\n[Liangliang]:\n"
         "\tThe input parameters of this function are:\n"
         "\t\t- points_3d: An array of 3D points (input to this function)\n"
         "\t\t- points_2d: An array of 2D image points (input to this "
         "function)\n"
         "\tThis function must return either 'true' on success or 'false' "
         "otherwise. On success, the camera\n"
         "\tparameters are returned by the following variables:\n"
         "\t\t- fx and fy: the focal lengths\n"
         "\t\t- cx and cy: the principal point\n"
         "\t\t- s: the skew factor, i.e., s = -alpha * cot(theta)\n"
         "\t\t- R: the 3x3 rotation matrix encoding camera orientation\n"
         "\t\t- t: a 3D vector encoding camera location.\n"
         "\tIMPORTANT: don't forget to write your recovered parameters to the "
         "above variables."
      << std::endl;

  // TODO: check if input is valid (e.g., number of correspondences >= 6, sizes
  // of 2D/3D points must match)
  if (end(points_3d) - begin(points_3d) < 6 ||
      end(points_2d) - begin(points_2d) < 6) {
    return false;
  }

  // TODO: construct the P matrix (so P * m = 0).
  Matrix P(2 * points_3d.size(), 12, 0.0);

  for (int i = 0; i < points_3d.size(); i++) {
    double X = points_3d[i].x();
    double Y = points_3d[i].y();
    double Z = points_3d[i].z();
    double x = points_2d[i].x();
    double y = points_2d[i].y();
    P(2 * i, 0) = X;
    P(2 * i, 1) = Y;
    P(2 * i, 2) = Z;
    P(2 * i, 3) = 1;
    P(2 * i, 8) = -x * X;
    P(2 * i, 9) = -x * Y;
    P(2 * i, 10) = -x * Z;
    P(2 * i, 11) = -x;

    P(2 * i + 1, 4) = X;
    P(2 * i + 1, 5) = Y;
    P(2 * i + 1, 6) = Z;
    P(2 * i + 1, 7) = 1;
    P(2 * i + 1, 8) = -y * X;
    P(2 * i + 1, 9) = -y * Y;
    P(2 * i + 1, 10) = -y * Z;
    P(2 * i + 1, 11) = -y;
  }

  std::cout << "p0 " << P.get(0, 0) << " " << P.get(0, 1) << " " << P.get(0, 2)
            << " " << P.get(0, 3) << " " << P.get(0, 4) << " " << P.get(0, 5)
            << " " << P.get(0, 6) << " " << P.get(0, 7) << " " << P.get(0, 8)
            << " " << P.get(0, 9) << " " << P.get(0, 10) << " " << P.get(0, 11)
            << std::endl;
  std::cout << "p0 " << P.get(1, 0) << " " << P.get(1, 1) << " " << P.get(1, 2)
            << " " << P.get(1, 3) << " " << P.get(1, 4) << " " << P.get(1, 5)
            << " " << P.get(1, 6) << " " << P.get(1, 7) << " " << P.get(1, 8)
            << " " << P.get(1, 9) << P.get(1, 11) << " " << P.get(1, 11)
            << std::endl;

  // TODO: solve for M (the whole projection matrix, i.e., M = K * [R, t])
  // using SVD decomposition.
  //   Optional: you can check if your M is correct by applying M on the 3D
  //   points. If correct, the projected point
  //             should be very close to your input images points.
  Matrix U(2 * points_3d.size(), 2 * points_3d.size(), 0.0);
  Matrix D(2 * points_3d.size(), 12, 0.0);
  Matrix V(12, 12, 0.0);
  svd_decompose(P, U, D, V);
  Vector vec_M(V.get_column(11));

  Matrix M(3, 4, vec_M.data());

  //  TODO: check svd
  double sum_error = 0;
  for (size_t i = 0; i < points_3d.size(); i++) {
    Vector3D point = points_3d[i];
    Vector2D point2d = points_2d[i];
    Vector3D projected2d_h = M * point.homogeneous();
    Vector2D projected2d = projected2d_h.cartesian();
    auto diff = (point2d - projected2d).length2();
    //    double x_res = M.get(0, 0) * point.x() + M.get(0, 1) * point.y() +
    //                   M.get(0, 2) * point.z() + M.get(0, 3) -
    //                   M.get(2, 0) * point2d.x() * point.x() -
    //                   M.get(2, 1) * point2d.x() * point.y() -
    //                   M.get(2, 2) * point2d.x() * point.z() - M.get(2, 3);
    //    double y_res = M.get(1, 0) * point.x() + M.get(1, 1) * point.y() +
    //                   M.get(1, 2) * point.z() + M.get(1, 3) -
    //                   M.get(2, 0) * point2d.y() * point.x() -
    //                   M.get(2, 1) * point2d.y() * point.y() -
    //                   M.get(2, 2) * point2d.y() * point.z() - M.get(2, 3);
    sum_error += diff;
  }
  double mean_reprojection_error = sum_error / points_3d.size();
  std::cout << "Mean Reprojection Error: " << mean_reprojection_error
            << std::endl;

  // decompose M into A and b M = [A b]
  Matrix33 A;
  A.set_row(0, Vector3D(M.get(0, 0), M.get(0, 1), M.get(0, 2)));
  A.set_row(1, Vector3D(M.get(1, 0), M.get(1, 1), M.get(1, 2)));
  A.set_row(2, Vector3D(M.get(2, 0), M.get(2, 1), M.get(2, 2)));
  Vector3D b(M.get(0, 3), M.get(1, 3), M.get(2, 3));

  // TODO: extract intrinsic parameters from M.
  //  TODO: check if a3 is not zero
  //  MEMO: sign of rho should be positive so that object comes in front of the
  //  camera
  auto p = 1 / norm(A.get_row(2));

  auto a1 = A.get_row(0);
  auto a2 = A.get_row(1);
  auto a3 = A.get_row(2);
  auto cos_theta = (dot(cross(a1, a3), cross(a2, a3))) /
                   (norm(cross(a1, a3)) * norm(cross(a2, a3)));
  if (norm(cross(a1, a3)) == 0 || norm(cross(a2, a3)) == 0) {
    std::cout << "cross product is zero" << std::endl;
    return false;
  }
  auto sin_theta = sqrt(1 - cos_theta * cos_theta);
  auto alpha = p * p * norm(cross(a1, a3)) * sin_theta;
  auto beta = p * p * norm(cross(a2, a3)) * sin_theta;

  fx = alpha;
  fy = beta / sin_theta;
  s = -alpha * (cos_theta / sin_theta);
  cx = p * p * dot(a1, a3);
  cy = p * p * dot(a2, a3);
  // print out intrinsic parameters
  std::cout << "p (rho) : " << p << std::endl;
  std::cout << "fx: " << fx << " fy: " << fy << " cx: " << cx << " cy: " << cy
            << " s: " << s << " cos theta: " << cos_theta
            << " sin theta: " << sin_theta << std::endl;

  // TODO: extract extrinsic parameters from M.
  if (norm(cross(a2, a3)) == 0) {
    std::cout << "cross product is zero" << std::endl;
    return false;
  }
  auto r1 = cross(a2, a3) / norm(cross(a2, a3));
  auto r3 = p * a3;
  auto r2 = cross(r3, r1);
  R.set_row(0, r1);
  R.set_row(1, r2);
  R.set_row(2, r3);

  Matrix33 K(alpha, s, cx, 0, fy, cy, 0, 0, 1);
  t = p * inverse(K) * b;
  // print out extrinsic parameters
  std::cout << "R: \n" << R << std::endl;
  std::cout << "t: " << t << std::endl;

  std::cout << "\n\tTODO: After you implement this function, please return "
               "'true' - this will trigger the viewer to\n"
               "\t\tupdate the rendering using your recovered camera "
               "parameters. This can help you to visually check\n"
               "\t\tif your calibration is successful or not.\n\n"
            << std::flush;
  return true;
}
