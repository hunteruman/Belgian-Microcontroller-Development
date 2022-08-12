/*
 * MECOTRON TUTORIAL
 *
 * This is a template to get you started in the course of the tutorial on the
 * control theory platforms, a.k.a. the MECOtrons.s
 * The tasks of the tutorial session will guide you through this template and
 * ask you to make use of the platform's capabilities in a step-by-step fashion.
 *
 * Every function in this template comes with an opening comment that describes
 * its purpose and functionality. Please also pay attention to the remarks that
 * are made in comment blocks.
 *
 */

#include "robot.h"

bool Robot::init() {
  MECOtron::init(); // Initialize the MECOtron

  return true;
}

int tp = 0;
float vi = 0;
void Robot::control() {
  
  float encC_value = getSpeedMotorA();
  float encD_value = getSpeedMotorB();

  // Compute update of motor voltages if controller is enabled (triggered by
  // pushing 'Button 0' in QRoboticsCenter)
  if(controlEnabled()) {
    // Fill your control law here to conditionally update the motor voltage...
    LED1(ON);
    LED2(OFF);

    tp += 1;

    //step input function
    int ii = tp%600;
    //if( ii < 150){
    //  vi = readValue(11);
    //}
    if (ii < 300){
      vi = readValue(11);
    }
    //else if (ii < 450){
    //  vi = -1*readValue(11);
    //}
    else {
      vi = -1*readValue(11);
    }
    

    float r = vi;          //  use float channel 0 from QRC as the reference angular velocity (in radians/second)
    float eA = r-encC_value;                 //  calculate the velocity error of motor A (in radians/second)
    float eB = r-encD_value;                 //  calculate the velocity error of motor B (in radians/second)

    //PI Control Algorithm
    float uA = controlA[0] + (K*Ti+K*Ts)*(eA)/(Ti)- K*errorA[0];
    float uB = controlB[0] + (K*Ti+K*Ts)*(eB)/(Ti)- K*errorB[0];

    // store the new errors and new control signals in a member variable
    for(int k=0; k<1; k++) {
      errorA[k+1] = errorA[k];        // shift the memorized errors of motor A with 1 sample
      errorB[k+1] = errorB[k];        // shift the memorized errors of motor B with 1 sample
      controlA[k+1] = controlA[k];    // shift the memorized control signals of motor A with 1 sample
      controlB[k+1] = controlB[k];    // shift the memorized control signals of motor B with 1 sample
    }
    errorA[0] = eA; errorB[0] = eB; controlA[0] = uA; controlB[0] = uB;    // append the new values

    // apply the control signal
    setVoltageMotorA(uA);
    setVoltageMotorB(uB);

    // send errors and control signals to QRC
    writeValue(2, r);     // reference
    writeValue(3, eA);    // should go to zero
    writeValue(4, eB);    // should go to zero
    writeValue(5, uA);    // should (preferably) remain between -12V and 12V
    writeValue(6, uB);    // should (preferably) remain between -12V and 12V
    
  } else {
    // If the controller is disabled, you might want to do something else...
    LED1(OFF);
    LED2(ON);
    setVoltageMotorA(0.0); // Apply 0.0 volts to motor A if the control is disabled
    setVoltageMotorB(0.0); // Apply 0.0 volts to motor B if the control is disabled
  }
  
  writeValue(0, encC_value);
  writeValue(1, encD_value);
  

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

bool Robot::controlEnabled() {
  return _button_states[0];       // The control is enabled if the state of button 0 is true
}

void Robot::button0callback() {
  if(toggleButton(0)) {                          // Switches the state of button 0 and checks if the new state is true
    resetController();
    message("Controller reset and enabled.");    // Display a message in the status bar of QRoboticsCenter
  }
  else {
    message("Robot disabled.");
  }
}

void Robot::button1callback() {
  toggleButton(1);
  init();                         // Reset the MECOtron and reinitialize the Robot object
  message("Reset.");
}
