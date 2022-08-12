/*
 * KALMAN FILTER TEMPLATE
 *
 * This is a template to get you started with the implementation of the Kalman filter
 * on your own cart.
 *
 */

#include "robot.h"

bool Robot::init() {
  MECOtron::init(); // Initialize the MECOtron

  return true;
}

int tp = 0;       //time counter
float p1 = 0.0;   //position input
void Robot::control() {

  float encC_value = getSpeedMotorA();
  float encD_value = getSpeedMotorB();

  float volt_A = 0.0;
  float volt_B = 0.0;
  Matrix<1> desired_velocity = 0; // defined scalar value for the input of the velocity controller

      tp += 1;

    //step input function
    //int ii = tp%420;
    //if( ii < 150){
     // p1 = readValue(0);
    //}
    //if (ii < 200){
      p1 = readValue(0);
    //}
    //else if (ii < 450){
    //  p1 = readValue(0);
    //}
    //else {
    //  p1 = readValue(0)-0.2;
    //}

  if(controlEnabled()) {   // only do this if controller is enabled (triggered by pushing 'Button 0' in QRoboticsCenter)

     // UNCOMMENT AND COMPLETE LINES BELOW TO IMPLEMENT POSITION CONTROLLER
    float desired_position = p1;                // use channel 0 to provide the constant position reference
    xref(0) = desired_position ;                // transform desired_position to the state reference (make sure units are consistent)
    K(0) = 2 ;                                  // state feedback gain K, to design
    desired_velocity = K * (xref - _xhat);      // calculate the state feedback signal, (i.e. the input for the velocity controller)

    // UNCOMMENT AND COMPLETE LINES BELOW TO IMPLEMENT VELOCITY CONTROLLER

    float vi = desired_velocity(0)/(0.033); //  converting linear velocity to angular for the velocity control
    float r = vi;                              //  reference angular velocity (in radians/second)
    float eA = r-encC_value;                   //  calculate the velocity error of motor A (in radians/second)
    float eB = r-encD_value;                   //  calculate the velocity error of motor B (in radians/second)

    // PI Control Algorithm
    volt_A = controlA[0] + (K_v*Ti+K_v*Ts)*(eA)/(Ti)- K_v*errorA[0];
    volt_B = controlB[0] + (K_v*Ti+K_v*Ts)*(eB)/(Ti)- K_v*errorB[0];

    // store the new errors and new control signals in a member variable
    for(int k=0; k<1; k++) {
      errorA[k+1] = errorA[k];        // shift the memorized errors of motor A with 1 sample
      errorB[k+1] = errorB[k];        // shift the memorized errors of motor B with 1 sample
      controlA[k+1] = controlA[k];    // shift the memorized control signals of motor A with 1 sample
      controlB[k+1] = controlB[k];    // shift the memorized control signals of motor B with 1 sample
    }
    errorA[0] = eA; errorB[0] = eB; controlA[0] = volt_A; controlB[0] = volt_B;    // append the new values


    // Send wheel speed command
    setVoltageMotorA(volt_A);
    setVoltageMotorB(volt_B);
    
  }
  else                      // do nothing since control is disables
  {
    desired_velocity(0) = 0.0;
    setVoltageMotorA(0.0);
    setVoltageMotorB(0.0);
  }

  // Kalman filtering
  if(KalmanFilterEnabled()) {

    // Prediction step
    TimeUpdate(desired_velocity, _xhat, _Phat);                         // do the prediction step -> update _xhat and _Phat

    // Correction step
    Matrix<1> distance_measurement;                                     // define a vector of length 1
    distance_measurement(0) = getFrontDistance();                       // front distance
    MeasurementUpdate(distance_measurement, _xhat, _Phat, _nu, _S);     // do the correction step -> update _xhat, _Phat, _nu, _S

  }

  // Send useful outputs to QRC
  writeValue(0, p1); //u
  writeValue(1, volt_B); //u
  writeValue(2, desired_velocity(0));
  writeValue(3, getPositionMotorA());
  ///writeValue(5, 0.0003); //R
  writeValue(4, getPositionMotorB());
  writeValue(5, getSpeedMotorA());
  writeValue(6, getSpeedMotorB());
  writeValue(7, getFrontDistance()); //y
  writeValue(8, _xhat(0)); //x
  writeValue(9, _Phat(0));
  writeValue(10, _nu(0));
  //float L = _Phat(0)*-1/_S(0); // Computation of L_k+1 as in the Measurement Update
  //writeValue(10, L);
  writeValue(11, _S(0));
  //writeValue(11, p1); 

}

void Robot::resetController(){
  // Set all errors and control signals in the memory back to 0
  for(int k=0; k<2; k++) {
    errorA[k] = 0.0;
    errorB[k] = 0.0;
    controlA[k] = 0.0;
    controlB[k] = 0.0;
  }
}

void Robot::resetKalmanFilter() {
   // UNCOMMENT AND MODIFIES LINES BELOW TO IMPLEMENT THE RESET OF THE KALMAN FILTER
   // Initialize state covariance matrix
   _Phat.Fill(0);            // Initialize the covariance matrix
   _Phat(0,0) = 0.000003;     // Fill the initial covariance matrix, you can change this according to your experiments
  
   // Initialize state estimate
   _xhat(0) = -0.2;     // Change this according to your experiments
}

bool Robot::controlEnabled() {
  return _button_states[0];       // The control is enabled if the state of button 0 is true
}

bool Robot::KalmanFilterEnabled() {
  return _button_states[1];
}

void Robot::button0callback() {
  if(toggleButton(0)) {           // Switches the state of button 0 and checks if the new state is true
    resetController();
    message("Controller reset and enabled.");    // Display a message in the status bar of QRoboticsCenter
  }
  else {
    message("Control disabled.");
  }
}

void Robot::button1callback() {
  if(toggleButton(1)){
      init();                         // Reset the MECOtron and reinitialize the Robot object
      resetKalmanFilter();            // Reset the Kalman filter
      message("Kalman filter reset and enabled.");
  }
  else
  {
    message("Kalman filter disabled.");
  }
}
