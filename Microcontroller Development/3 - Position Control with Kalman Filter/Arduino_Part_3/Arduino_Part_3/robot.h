#ifndef ROBOT_H
#define ROBOT_H

/*
 * ROBOT Class
 *
 * Class incorporating the robot. This class is used to define state machines,
 * control algorithms, sensor readings,...
 * It should be interfaced with the communicator to send data to the world.
 *
 */

#include "mecotron.h" // Include MECOTRON header
#include <BasicLinearAlgebra.h> // Include BasicLinearAlgebra to make matrix manipulations easier
#include "kalman_filter.h" // Include template to make Kalman filter implementation easier

class Robot : public MECOtron {
  private:

    // Member variables
      // Remember the 2 previous errors
      float errorA[2] = {0.0, 0.0};
      float errorB[2] = {0.0, 0.0};

      // Remember the 2 previous control signals
      float controlA[2] = {0.0, 0.0};
      float controlB[2] = {0.0, 0.0};

      // Velocity Controller parameters
      const float Ti = 0.05;
      const float Ts = TSAMPLE;
      const float K_v = 0.2506;

    // Class variables

    // Kalman filter
    Matrix<1> _xhat;      // state estimate vector
    Matrix<1,1> _Phat;    // state estimate covariance
    Matrix<1> _nu;        // innovation vector
    Matrix<1,1> _S;    // innovation covariance

    // Position controller
    Matrix<1> xref;       // reference state
    Matrix<1,1> K;        // state feedback gain

  public:
    // Constructor
    Robot() { }

    void control();

    // General functions
    bool init();  // Set up the robot

    bool controlEnabled();
    bool KalmanFilterEnabled();

    void resetController();
    void resetKalmanFilter();

    void button0callback();
    void button1callback();

};

#endif // ROBOT_H
