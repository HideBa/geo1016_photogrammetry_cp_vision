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

#include "matrix_algo.h"
#include "triangulation.h"
#include <easy3d/optimizer/optimizer_lm.h>

using namespace easy3d;

Matrix33
normalisation_transformation_matrix(const std::vector<Vector2D> &points);
std::vector<Vector2D> normalize_points(const std::vector<Vector2D> points,
                                       Matrix33 T);

Matrix33 estimate_fundamental_matrix(const Matrix &W);

std::tuple<Matrix33, Matrix33, Vector3D> recover_extrinsic(const Matrix33 &E);
Vector3D reconstruct(const Matrix34 &M, const Matrix34 &M_prime,
                     const Vector2D &p0, const Vector2D &p1);
Vector3D reconstruct_nonlinear(const Matrix34 &M, const Matrix34 &M_prime,
                               const Vector2D &p0, const Vector2D &p1);
bool is_in_front_of_camera(const Vector3D &point, const Matrix33 &R,
                           const Vector3D &t);

double check_matrix(const Matrix33 &Mat, std::vector<Vector2D> points_0,
                    std::vector<Vector2D> points_1);

double evaluate_3dcood(const Matrix34 &M, std::vector<Vector2D> original_points,
                       std::vector<Vector3D> reconstructed_points);
/**
 * TODO: Finish this function for reconstructing 3D geometry from
 * corresponding image points.
 * @return True on success, otherwise false. On success, the reconstructed
 * 3D points must be written to 'points_3d' and the recovered relative pose
 * must be written to R and t.
 */
bool Triangulation::triangulation(
    double fx,
    double fy, /// input: the focal lengths (same for both cameras)
    double cx,
    double cy, /// input: the principal point (same for both cameras)
    double s,  /// input: the skew factor (same for both cameras)
    const std::vector<Vector2D>
        &points_0, /// input: 2D image points in the 1st image.
    const std::vector<Vector2D>
        &points_1, /// input: 2D image points in the 2nd image.
    std::vector<Vector3D> &points_3d, /// output: reconstructed 3D points
    Matrix33 &R, /// output: 3 by 3 matrix, which is the recovered rotation
                 /// of the 2nd camera
    Vector3D &t  /// output: 3D vector, which is the recovered translation of
                 /// the 2nd camera
) const {
  /// NOTE: there might be multiple workflows for reconstructing 3D geometry
  /// from corresponding image points.
  ///       This assignment uses the commonly used one explained in our
  ///       lecture. It is advised to define a function for the sub-tasks.
  ///       This way you have a clean and well-structured implementation,
  ///       which also makes testing and debugging easier. You can put your
  ///       other functions above triangulation(), or put them in one or
  ///       multiple separate files.

  std::cout << "\nTODO: I am going to implement the triangulation() function "
               "in the following file:"
            << std::endl
            << "\t    - triangulation_method.cpp\n\n";

  std::cout << "[Liangliang]:\n"
               "\tFeel free to use any provided data structures and functions. "
               "For your convenience, the\n"
               "\tfollowing three files implement basic linear algebra data "
               "structures and operations:\n"
               "\t    - Triangulation/matrix.h  Matrices of arbitrary "
               "dimensions and related functions.\n"
               "\t    - Triangulation/vector.h  Vectors of arbitrary "
               "dimensions and related functions.\n"
               "\t    - Triangulation/matrix_algo.h  Determinant, inverse, "
               "SVD, linear least-squares...\n"
               "\tPlease refer to the above files for a complete list of "
               "useful functions and their usage.\n\n"
               "\tIf you choose to implement the non-linear method for "
               "triangulation (optional task). Please\n"
               "\trefer to 'Tutorial_NonlinearLeastSquares/main.cpp' for an "
               "example and some explanations.\n\n"
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

  /// For more functions of Matrix and Vector, please refer to 'matrix.h' and
  /// 'vector.h'

  // TODO: delete all above example code in your final submission

  //--------------------------------------------------------------------------------------------------------------
  // implementation starts ...

  // TODO: check if the input is valid (always good because you never known
  // how others will call your function).
  if (fx < 0 || fy < 0) {
    std::cout << "invalid fx and fy" << std::endl;
    return false;
  }
  for (auto p : points_0) {
    assert(p.x() >= 0);
    assert(p.y() >= 0);
  }
  for (auto p : points_1) {
    assert(p.x() >= 0);
    assert(p.y() >= 0);
  }
  assert(points_0.size() == points_1.size() && points_0.size() >= 8);

  // TODO: Estimate relative pose of two views. This can be subdivided into
  //      - estimate the fundamental matrix F;
  //      - compute the essential matrix E;
  //      - recover rotation R and t.

  // Normalise points
  auto T = normalisation_transformation_matrix(points_0);
  auto T_prime = normalisation_transformation_matrix(points_1);
  auto normalised_points_0 = normalize_points(points_0, T);

  auto normalised_points_1 = normalize_points(points_1, T_prime);

  Matrix Wq(points_0.size(), 9, 0.0);
  for (size_t i = 0; i < normalised_points_0.size(); i++) {
    auto p = normalised_points_0[i];
    auto p_prime = normalised_points_1[i];
    auto u = p.x(), v = p.y();
    auto u_prime = p_prime.x(), v_prime = p_prime.y();
    Wq.set_row(i, {u * u_prime, v * u_prime, u_prime, u * v_prime, v * v_prime,
                   v_prime, u, v, 1});
  }

  Matrix33 Fq = estimate_fundamental_matrix(Wq);

  //  Debugging--------
  auto fq_error = check_matrix(Fq, normalised_points_0, normalised_points_1);
  std::cout << "fq error: " << fq_error << std::endl;
  //  Debugging--------

  // Fundamental matrix
  Matrix33 F = transpose(T_prime) * Fq * T;
  F /= F.get(2, 2);

  //  Debugging--------
  auto f_error = check_matrix(F, points_0, points_1);
  std::cout << "f_error: " << f_error << std::endl;
  //  Debugging--------

  // Essential matrix
  Matrix33 K(fx, s, cx, 0, fy, cy, 0, 0, 1);
  Matrix33 E = transpose(K) * F * K;

  auto e_error = check_matrix(E, points_0, points_1);

  Matrix33 R1;
  Matrix33 R2;
  Vector3D temp_t;
  std::tie(R1, R2, temp_t) = recover_extrinsic(E);
  std::vector<std::pair<Matrix33, Vector3D>> R_t_pairs(
      {{R1, temp_t}, {R1, -temp_t}, {R2, temp_t}, {R2, -temp_t}});

  std::vector<int> positive_z_count;
  std::vector<std::vector<Vector3D>> reconstructed_points_list;
  std::vector<std::pair<Matrix34, Matrix34>> M_pairs;
  for (auto pair : R_t_pairs) {
    auto R = pair.first;
    auto t = pair.second;
    Matrix34 extrinsic(1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
                       0.0);
    Matrix34 extrinsic_prime(R.get(0, 0), R.get(0, 1), R.get(0, 2), t.x(),
                             R.get(1, 0), R.get(1, 1), R.get(1, 2), t.y(),
                             R.get(2, 0), R.get(2, 1), R.get(2, 2), t.z());
    Matrix34 M = K * extrinsic;
    Matrix34 M_prime = K * extrinsic_prime;
    M_pairs.push_back(std::pair<Matrix34, Matrix34>(M, M_prime));
    std::vector<Vector3D> reconstructed_points;
    for (size_t i = 0; i < points_0.size(); i++) {
      //      auto point_3d = reconstruct(M, M_prime, points_0[i], points_1[i]);
      auto point_3d =
          reconstruct_nonlinear(M, M_prime, points_0[i], points_1[i]);
      reconstructed_points.push_back(point_3d);
    }

    reconstructed_points_list.push_back(reconstructed_points);
    // calculate the number of points in front of the camera
    auto positive_z = 0;
    for (auto p : reconstructed_points) {
      if (is_in_front_of_camera(p, R, t) && p.z() > 0) {
        positive_z++;
      }
    }
    std::cout << "num of positive z " << positive_z << " t: " << t << std::endl;
    positive_z_count.push_back(positive_z);
  }

  int best_pair_index = std::distance(
      positive_z_count.begin(),
      std::max_element(positive_z_count.begin(), positive_z_count.end()));
  std::tie(R, t) = R_t_pairs[best_pair_index];
  points_3d = reconstructed_points_list[best_pair_index];

  auto best_M_pair = M_pairs[best_pair_index];

  auto mean_error1 = evaluate_3dcood(best_M_pair.first, points_0, points_3d);
  auto mean_error2 = evaluate_3dcood(best_M_pair.second, points_1, points_3d);

  std::cout << "best R: " << R << "\n";
  std::cout << "best t: " << t << "\n";
  std::cout << "error of first image: " << mean_error1 << "\n";
  std::cout << "error of second image: " << mean_error2 << "\n";

  // TODO: Reconstruct 3D points. The main task is
  //      - triangulate a pair of image points (i.e., compute the 3D
  //      coordinates for each corresponding point pair)

  // TODO: Don't forget to
  //          - write your recovered 3D points into 'points_3d' (so the
  //          viewer can visualize the 3D points for you);
  //          - write the recovered relative pose into R and t (the view
  //          will be updated as seen from the 2nd camera,
  //            which can help you check if R and t are correct).
  //       You must return either 'true' or 'false' to indicate whether the
  //       triangulation was successful (so the viewer will be notified to
  //       visualize the 3D points and update the view). There are a few
  //       cases you should return 'false' instead, for example:
  //          - function not implemented yet;
  //          - input not valid (e.g., not enough points, point numbers
  //          don't match);
  //          - encountered failure in any step.
  return points_3d.size() > 0;
}

std::vector<Vector2D> normalize_points(const std::vector<Vector2D> points,
                                       Matrix33 T) {
  std::vector<Vector2D> normalised_points;
  for (const auto p : points) {
    auto homogeneous_p = p.homogeneous();
    Vector3D normalised_p = T * homogeneous_p;
    normalised_points.push_back(normalised_p.cartesian());
  }
  return normalised_points;
}

Matrix33
normalisation_transformation_matrix(const std::vector<Vector2D> &points) {
  double p_x_sum = 0.0, p_y_sum = 0.0;
  int count = 0;
  for (auto p : points) {
    count++;
    p_x_sum += p.x();
    p_y_sum += p.y();
  }
  Vector2D centroid(p_x_sum / points.size(), p_y_sum / points.size());

  double dist_sum = 0.0;
  for (auto p : points) {
    double dx = p.x() - centroid.x();
    double dy = p.y() - centroid.y();
    dist_sum += sqrt(dx * dx + dy * dy);
  }
  double mean_dist = dist_sum / points.size();

  double scale = sqrt(2) / mean_dist;

  Matrix33 T(0.0);
  T.set(0, 0, scale);
  T.set(1, 1, scale);
  T.set(0, 2, -scale * centroid.x());
  T.set(1, 2, -scale * centroid.y());
  T.set(2, 2, 1);
  return T;
};

Matrix33 estimate_fundamental_matrix(const Matrix &W) { // W is Nx9 matrix
  Matrix U(W.rows(), W.rows(), 0.0);
  Matrix D(W.rows(), 9, 0.0);
  Matrix V(9, 9, 0.0);
  svd_decompose(W, U, D, V);
  Vector vec_F_hat(V.get_column(8));
  Matrix F_hat(3, 3, vec_F_hat.data());

  U = Matrix33(0.0);
  D = Matrix33(0.0);
  V = Matrix33(0.0);
  svd_decompose(F_hat, U, D, V);
  D.set(2, 2, 0.0); // This is for the best rank-2 approximation
  Matrix33 F(0.0);
  F = U * D * transpose(V);

  return F;
}

std::tuple<Matrix33, Matrix33, Vector3D> recover_extrinsic(const Matrix33 &E) {
  Matrix33 U(0.0);
  Matrix33 D(0.0);
  Matrix33 V(0.0);
  svd_decompose(E, U, D, V);
  Matrix33 W(0, -1, 0, 1, 0, 0, 0, 0, 1);
  Matrix33 Z(0, 1, 0, -1, 0, 0, 0, 0, 0);
  Matrix33 R1 = U * W * transpose(V);
  Matrix33 R2 = U * transpose(W) * transpose(V);

  if (determinant(R1) < 0) {
    R1 = -R1;
  }
  if (determinant(R2) < 0) {
    R2 = -R2;
  }

  auto t = U.get_column(2);
  return std::tuple<Matrix33, Matrix33, Vector3D>(R1, R2, t);
}

Vector3D reconstruct(const Matrix34 &M, const Matrix34 &M_prime,
                     const Vector2D &p0, const Vector2D &p1) {

  Matrix44 A;
  A.set_row(0, p0.x() * M.get_row(2) - M.get_row(0));
  A.set_row(1, p0.y() * M.get_row(2) - M.get_row(1));
  A.set_row(2, p1.x() * M_prime.get_row(2) - M_prime.get_row(0));
  A.set_row(3, p1.y() * M_prime.get_row(2) - M_prime.get_row(1));

  Matrix U(4, 4);
  Matrix D(4, 4);
  Matrix V(4, 4);
  svd_decompose(A, U, D, V);
  Vector4D h_P = Vector4D(V.get_column(3));
  Vector3D P = h_P.cartesian();
  return P;
}

bool is_in_front_of_camera(const Vector3D &point, const Matrix33 &R,
                           const Vector3D &t) {
  Vector3D point_camera = R * point + t;
  return point_camera.z() > 0;
}

double check_matrix(const Matrix33 &Mat, std::vector<Vector2D> points_0,
                    std::vector<Vector2D> points_1) {
  auto sum_error = 0.0;
  for (size_t i = 0; i < points_0.size(); i++) {
    auto point0 = points_0[i];
    auto point1 = points_1[i];
    auto point0_h = point0.homogeneous();
    auto point1_h = point1.homogeneous();
    Matrix point0_h_m(3, 1, {point0_h.x(), point0_h.y(), point0_h.z()});
    Matrix point1_h_m(3, 1, {point1_h.x(), point1_h.y(), point1_h.z()});
    auto error = norm((transpose(point1_h_m) * Mat * point0_h_m));
    sum_error += error;
  }
  auto ave_error = sum_error / points_0.size();
  return ave_error;
}
double evaluate_3dcood(const Matrix34 &M, std::vector<Vector2D> original_points,
                       std::vector<Vector3D> reconstructed_points) {
  double sum_error = 0;
  for (size_t i = 0; i < reconstructed_points.size(); i++) {
    auto p = reconstructed_points[i];
    Vector3D reprojected_h = M * p.homogeneous();
    auto reprojected = reprojected_h.cartesian();
    auto diff = norm(original_points[i] - reprojected);
    sum_error += diff;
  }
  double ave_error = sum_error / reconstructed_points.size();
  return ave_error;
}

class Objective : public Objective_LM {
public:
  Objective(int num_func, int num_var, const Matrix34 &M,
            const Matrix34 &M_prime, const Vector2D &point_0,
            const Vector2D &point_1)
      : Objective_LM(num_func, num_var), M(M), M_prime(M_prime),
        point_0(point_0), point_1(point_1) {}

protected:
  const Matrix34 &M;
  const Matrix34 &M_prime;
  const Vector2D &point_0;
  const Vector2D &point_1;

  int evaluate(const double *x, double *fvec) {
    Vector3D estimated_point(x[0], x[1], x[2]);
    Vector3D reprojected0_h = this->M * estimated_point.homogeneous();
    Vector2D reprojected0 = reprojected0_h.cartesian();
    Vector3D reprojected1_h = this->M_prime * estimated_point.homogeneous();
    Vector2D reprojected1 = reprojected1_h.cartesian();

    auto diff_x_1 = reprojected0.x() - point_0.x();
    auto diff_x_2 = reprojected1.x() - point_1.x();
    auto diff_y_1 = reprojected0.y() - point_0.y();
    auto diff_y_2 = reprojected1.y() - point_1.y();

    fvec[0] = diff_x_1 * diff_x_1 + diff_y_1 + diff_y_1;
    fvec[1] = diff_x_2 * diff_x_2 + diff_y_2 * diff_y_2;
    return 0;
  }
};

Vector3D reconstruct_nonlinear(const Matrix34 &M, const Matrix34 &M_prime,
                               const Vector2D &p0, const Vector2D &p1) {
  Objective obj(4, 3, M, M_prime, p0, p1);
  Optimizer_LM lm;
  std::vector<double> x = {4.0, -4.0, 4.0};
  bool status = lm.optimize(&obj, x);
  std::cout << "status: " << status << std::endl;
  /// retrieve the result.
  std::cout << "the solution is:     " << x[0] << ", " << x[1] << ", " << x[2]
            << std::endl;

  return Vector3D(x[0], x[1], x[2]);
}
