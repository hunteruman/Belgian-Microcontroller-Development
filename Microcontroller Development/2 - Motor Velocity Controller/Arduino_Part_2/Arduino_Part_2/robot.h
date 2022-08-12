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

class Robot : public MECOtron {
  private:
    // Member variables
      // Remember the 2 previous errors
      float errorA[2] = {0.0, 0.0};
      float errorB[2] = {0.0, 0.0};

      // Remember the 2 previous control signals
      float controlA[2] = {0.0, 0.0};
      float controlB[2] = {0.0, 0.0};

      // Controller parameters
      const float Ti = 0.05;
      const float Ts = TSAMPLE;
      const float K = 0.2506;
      
  public:
    // Constructor
    Robot() { }

    void control();

    // General functions
    bool init();  // Set up the robot

    bool controlEnabled();

    void button0callback();
    void button1callback();

    // Controller related functions
    void resetController(); 
};

#endif // ROBOT_H
