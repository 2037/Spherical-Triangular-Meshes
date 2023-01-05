// Liuming Wang 111164541
#include "manip.hpp"
#include <algorithm>
#include <cmath>

using namespace std;
namespace ams562_final {

// cross product
void crossProduct(const array<double, 3> A, const array<double, 3> B,
                  array<double, 3> &AxB) {
  AxB[0] = A[1] * B[2] - A[2] * B[1];
  AxB[1] = A[2] * B[0] - A[0] * B[2];
  AxB[2] = A[0] * B[1] - A[1] * B[0];
};

void arrayadd(const array<double, 3> vec1, const array<double, 3> vec2,
              array<double, 3> &res) {
  for (int i = 0; i < 3; i++) {
    res[i] = vec1[i] + vec2[i];
  }
};

void arraysub(const array<double, 3> vec1, const array<double, 3> vec2,
              array<double, 3> &res) {
  for (int i = 0; i < 3; i++) {
    res[i] = vec1[i] - vec2[i];
  }
};

template <class T>
inline std::ostream &operator<<(std::ostream &o_str,
                                const std::array<T, 3> adj) {
  for (unsigned i = 0; i < adj.size(); i++) {
    o_str << to_string(adj[i]) << '\n';
  }
  o_str << "\n";
  return o_str;
}

void compute_n2e_adj(const unsigned n, const Triangles &conn,
                     std::vector<std::vector<int>> &adj) {
  // resize adj to n
  adj.resize(n);
  for_each(adj.begin(), adj.end(), [](vector<int> &i) { i.reserve(10); });

  // FIRST GLOBAL LOOPS on triangles
  // index is the triangle ID of e
  int index = 0;
  for_each(conn.to_vector().begin(), conn.to_vector().end(),
           [&index, &adj](array<int, 3> e) {
             for_each(e.begin(), e.end(),[&index, &adj](int v) { 
                adj[v].push_back(index); 
             });
             index++;
           });

  for_each(adj.begin(), adj.end(), [](vector<int> &i) { i.shrink_to_fit(); });
}

void compute_avg_normals(const SphCo &points, const Triangles &conn,
                         const std::vector<std::vector<int>> &n2e_adj,
                         SphCo &nrms) {
  // resize the nrms
  nrms.resize(points.npoints());

  // initialization
  int n = points.npoints();
  int m = conn.ntris();

  // // Compute facial normal vectors
  // for (signed e = 0; e < m; e++) {
  //   A = points[conn[e][0]];
  //   B = points[conn[e][1]];
  //   C = points[conn[e][2]];
  //   // AC = C-A; AB = B-A
  //   for (int coor = 0; coor < 3; coor++) {
  //     AC[coor] = C[coor] - A[coor];
  //     AB[coor] = B[coor] - A[coor];
  //   }
  //   crossProduct(AB, AC, cp);
  //   ABxAC[e][0] = cp[0];
  //   ABxAC[e][1] = cp[1];
  //   ABxAC[e][2] = cp[2];
  // }

  // SECOND GLOBAL LOOPS on triangles
  array<double, 3> AC, AB, cp;
  SphCo ABxAC(m);
  int index = 0;
  for_each(conn.to_vector().begin(), conn.to_vector().end(),
           [points, &index, &AC, &AB, &cp, &ABxAC](array<int, 3> e) {
             arraysub(points[e[2]], points[e[0]], AC);
             arraysub(points[e[1]], points[e[0]], AB);
             crossProduct(AB, AC, cp);
             ABxAC.to_vector()[index] = cp;
             index++;
           });

  // Normalize facial normal vector n^f
  ABxAC.normalize();

  // Perform averaging
  int k = 0;
  array<double, 3> z = {0};
  int index2 = 0;
  double d = 0;
  // point on conn/triangles at node v
  for (signed v = 0; v < n; v++) {
    // since we used shrink_to_fit in compute_n2e_adj
    k = n2e_adj[v].size();

    // point on n^f
    // sum over all adjacent normals and average them
    z = {0, 0, 0};
    for (signed j = 0; j < k; j++) {
      index2 = n2e_adj[v][j];
      arrayadd(z, ABxAC[index2], z);
    }

    // z = sum/k
    // conpute distance d
    d = 0;
    for_each(z.begin(), z.end(), [&d, k](double &z_i) {
      z_i /= k;
      d += z_i * z_i;
    });
    d = sqrt(d);

    for_each(z.begin(), z.end(), [d](double &z_i) { z_i /= d; });

    // Assign z to nrms_v
    nrms[v] = z;
  }
}

void compute_errors(const SphCo &exact, const SphCo &num,
                    std::vector<double> &acos_theta) {
  // resize the error array
  acos_theta.resize(num.npoints());
  int n = num.npoints();

  const auto dot_product = [](const array<double, 3> a,
                              const array<double, 3> b) {
    return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
  };

  double angle = 0;
  for (int i = 0; i < n; i++) {
    angle = dot_product(exact[i], num[i]);
    acos_theta[i] = acos(angle);
  }
}

}  // namespace ams562_final
// Thu Dec  3 22:53:22 EST 2020