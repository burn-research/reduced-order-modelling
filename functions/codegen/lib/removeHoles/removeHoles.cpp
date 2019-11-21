//
// File: removeHoles.cpp
//
// MATLAB Coder version            : 3.2
// C/C++ source code generated on  : 01-Dec-2017 14:55:34
//

// Include Files
#include "rt_nonfinite.h"
#include "removeHoles.h"

// Function Definitions

//
// Description
//  This function removes the missing number in a series of integers.
//  Example: [1 2 4 5] becomes [1 2 3 4].
//  Input:
//  - x: (N x 1) array of integers
//  - y: (N x 1) array of integers
// Arguments    : double x_data[]
//                int x_size[2]
//                double y_data[]
//                int y_size[2]
// Return Type  : void
//
void removeHoles(double x_data[], int x_size[2], double y_data[], int y_size[2])
{
  int na;
  int n;
  int qEnd;
  int kEnd;
  int idx_data[7];
  int iwork_data[7];
  double u_data[7];
  int k;
  boolean_T p;
  int i;
  int q;
  int i2;
  int j;
  int pEnd;
  int b_p;
  double x;
  int exitg1;
  double m;
  int exponent;
  boolean_T I_data[7];
  double b_x_data[7];

  //  Main
  na = x_size[1];
  n = x_size[1] + 1;
  qEnd = x_size[1];
  for (kEnd = 0; kEnd < qEnd; kEnd++) {
    idx_data[kEnd] = 0;
  }

  if (x_size[1] != 0) {
    for (k = 1; k <= n - 2; k += 2) {
      if ((x_data[k - 1] <= x_data[k]) || rtIsNaN(x_data[k])) {
        p = true;
      } else {
        p = false;
      }

      if (p) {
        idx_data[k - 1] = k;
        idx_data[k] = k + 1;
      } else {
        idx_data[k - 1] = k + 1;
        idx_data[k] = k;
      }
    }

    if ((x_size[1] & 1) != 0) {
      idx_data[x_size[1] - 1] = x_size[1];
    }

    i = 2;
    while (i < n - 1) {
      i2 = i << 1;
      j = 1;
      for (pEnd = 1 + i; pEnd < n; pEnd = qEnd + i) {
        b_p = j;
        q = pEnd - 1;
        qEnd = j + i2;
        if (qEnd > n) {
          qEnd = n;
        }

        k = 0;
        kEnd = qEnd - j;
        while (k + 1 <= kEnd) {
          if ((x_data[idx_data[b_p - 1] - 1] <= x_data[idx_data[q] - 1]) ||
              rtIsNaN(x_data[idx_data[q] - 1])) {
            p = true;
          } else {
            p = false;
          }

          if (p) {
            iwork_data[k] = idx_data[b_p - 1];
            b_p++;
            if (b_p == pEnd) {
              while (q + 1 < qEnd) {
                k++;
                iwork_data[k] = idx_data[q];
                q++;
              }
            }
          } else {
            iwork_data[k] = idx_data[q];
            q++;
            if (q + 1 == qEnd) {
              while (b_p < pEnd) {
                k++;
                iwork_data[k] = idx_data[b_p - 1];
                b_p++;
              }
            }
          }

          k++;
        }

        for (k = 0; k + 1 <= kEnd; k++) {
          idx_data[(j + k) - 1] = iwork_data[k];
        }

        j = qEnd;
      }

      i = i2;
    }
  }

  for (k = 0; k + 1 <= na; k++) {
    u_data[k] = x_data[idx_data[k] - 1];
  }

  k = 0;
  while ((k + 1 <= na) && rtIsInf(u_data[k]) && (u_data[k] < 0.0)) {
    k++;
  }

  q = k;
  k = x_size[1];
  while ((k >= 1) && rtIsNaN(u_data[k - 1])) {
    k--;
  }

  qEnd = x_size[1] - k;
  while ((k >= 1) && rtIsInf(u_data[k - 1]) && (u_data[k - 1] > 0.0)) {
    k--;
  }

  pEnd = (x_size[1] - k) - qEnd;
  b_p = -1;
  if (q > 0) {
    b_p = 0;
  }

  i2 = (q + k) - q;
  while (q + 1 <= i2) {
    x = u_data[q];
    do {
      exitg1 = 0;
      q++;
      if (q + 1 > i2) {
        exitg1 = 1;
      } else {
        m = std::abs(x / 2.0);
        if ((!rtIsInf(m)) && (!rtIsNaN(m))) {
          if (m <= 2.2250738585072014E-308) {
            m = 4.94065645841247E-324;
          } else {
            frexp(m, &exponent);
            m = std::ldexp(1.0, exponent - 53);
          }
        } else {
          m = rtNaN;
        }

        if ((std::abs(x - u_data[q]) < m) || (rtIsInf(u_data[q]) && rtIsInf(x) &&
             ((u_data[q] > 0.0) == (x > 0.0)))) {
          p = true;
        } else {
          p = false;
        }

        if (!p) {
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);

    b_p++;
    u_data[b_p] = x;
  }

  if (pEnd > 0) {
    b_p++;
    u_data[b_p] = u_data[i2];
  }

  q = i2 + pEnd;
  for (j = 1; j <= qEnd; j++) {
    b_p++;
    u_data[b_p] = u_data[(q + j) - 1];
  }

  if (1 > b_p + 1) {
    kEnd = 0;
  } else {
    kEnd = b_p + 1;
  }

  //  Number of unique elements in x
  //  Loop over values in u
  b_p = (int)((1.0 + (-1.0 - ((double)kEnd - 1.0))) / -1.0);
  for (q = 0; q < b_p; q++) {
    i2 = (kEnd - q) - 2;

    //  u(ii) must equal u(ii+1) - 1
    if (u_data[i2] + 1.0 != u_data[i2 + 1]) {
      m = u_data[i2 + 1] - u_data[i2];

      //  Hole's size
      qEnd = x_size[0] * x_size[1];
      for (pEnd = 0; pEnd < qEnd; pEnd++) {
        I_data[pEnd] = (x_data[pEnd] > u_data[i2]);
      }

      //  Values in x who are greater than u(ii)
      i2 = x_size[1] - 1;
      qEnd = 0;
      for (i = 0; i <= i2; i++) {
        if (I_data[i]) {
          qEnd++;
        }
      }

      pEnd = 0;
      for (i = 0; i <= i2; i++) {
        if (I_data[i]) {
          idx_data[pEnd] = i + 1;
          pEnd++;
        }
      }

      for (pEnd = 0; pEnd < qEnd; pEnd++) {
        b_x_data[pEnd] = (x_data[idx_data[pEnd] - 1] - m) + 1.0;
      }

      for (pEnd = 0; pEnd < qEnd; pEnd++) {
        x_data[idx_data[pEnd] - 1] = b_x_data[pEnd];
      }

      //  Get rid of the hole (but add 1)
    }
  }

  //  Output
  y_size[0] = 1;
  y_size[1] = x_size[1];
  qEnd = x_size[0] * x_size[1];
  for (kEnd = 0; kEnd < qEnd; kEnd++) {
    y_data[kEnd] = x_data[kEnd];
  }
}

//
// File trailer for removeHoles.cpp
//
// [EOF]
//
