#ifndef KALMAN_FILTER_H
#define KALMAN_FILTER_H

#include <BasicLinearAlgebra.h>
#include "mecotron.h" // Include MECOTRON header

void TimeUpdate(const Matrix<1> &u, Matrix<1> &xhat, Matrix<1,1> &Phat);
void MeasurementUpdate(const Matrix<1> &y, Matrix<1> &xhat, Matrix<1,1> &Phat, Matrix<1> &nu, Matrix<1,1> &S);

#endif // KALMAN_FILTER_H
