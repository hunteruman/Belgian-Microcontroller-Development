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

  // Initializing the robot's specific variables
  for(int k=0; k<2; k++){
    x[k] = 0.0;   // Set all components of the vector (float array) x to 0 as initialization
  }
  return true;
}

int tp = 0;
void Robot::control() {
  
  tp += 1;
  
  float vA = readValue(0);
  float vB = readValue(1);
  //float encA_value = getPositionMotorA();
  //float encB_value = getPositionMotorB();
  float encC_value = getSpeedMotorA();
  float encD_value = getSpeedMotorB();

  // Compute update of motor voltages if controller is enabled (triggered by
  // pushing 'Button 0' in QRoboticsCenter)
  if(controlEnabled()) {
    // Fill your control law here to conditionally update the motor voltage...
    LED1(ON);
    LED2(OFF);
    //writeValue(3, encA_value); 
    //writeValue(4, encB_value);
    writeValue(5, encC_value); 
    writeValue(6, encD_value);
    writeValue(9, getPendulumAngle());
    int ii = tp%650;
    if( ii < 150){
        //setVoltageMotorA(vA); // Apply vA volts in the first phase of this period
        //setVoltageMotorB(vB); // Apply vB "
        setVoltageMotorA(vA); // One of two superimposed test functions
        setVoltageMotorB(vB); // Test function
        writeValue(7, vA);
        writeValue(8, vB);
    }
    else if (ii < 300){
        setVoltageMotorA(0.0); // One of two superimposed test functions
        setVoltageMotorB(0.0); // Test functions
        writeValue(7, 0.0);
        writeValue(8, 0.0);
    }
    else if (ii < 450){
        setVoltageMotorA(-vA); // One of two superimposed test functions
        setVoltageMotorB(-vB); // Test functions
        writeValue(7, -vA);
        writeValue(8, -vB);
    }
    //else if (ii < 100){
    //    setVoltageMotorA(-vA); // One of two superimposed test functions
    //    setVoltageMotorB(-vB); // Test functions
    //    writeValue(7, -vA);
    //    writeValue(8, -vB);
    //}
    //else if (ii < 125){
    //    setVoltageMotorA(-vA*2); // One of two superimposed test functions
    //    setVoltageMotorB(-vB*2); // Test functions
    //    writeValue(7, -vA*2);
    //    writeValue(8, -vB*2);
    //}
    //else if (ii < 150){
    //    //setVoltageMotorA(0.0); // Apply 0.0 volts to motor A if the control is disabled
    //    //setVoltageMotorB(0.0); // Apply 0.0 volts to motor B if the control is disabled
    //    setVoltageMotorA(-vA); // One of two superimposed test functions
    //    setVoltageMotorB(-vB); // Test functions
    //    writeValue(7, -vA);
    //    writeValue(8, -vB);
    //}
    //else if (ii < 190){
    //    //setVoltageMotorA(-vA); // Apply vA volts in the third phase of this period
    //    //setVoltageMotorB(-vB); // Apply vB "
    //    setVoltageMotorA(0.0); // One of two superimposed test functions
    //    setVoltageMotorB(0.0); // Test functions
    //    writeValue(7, 0.0);
    //   writeValue(8, 0.0);
    //}
    else {
      setVoltageMotorA(0.0); // Apply 0.0 volts to motor A if the control is disabled
      setVoltageMotorB(0.0); // Apply 0.0 volts to motor B if the control is disabled
      writeValue(7, 0.0);
      writeValue(8, 0.0);
    }
  } else {
    // If the controller is disabled, you might want to do something else...
    LED1(OFF);
    LED2(ON);
    setVoltageMotorA(0.0); // Apply 0.0 volts to motor A if the control is disabled
    setVoltageMotorB(0.0); // Apply 0.0 volts to motor B if the control is disabled
  }

  float va = getSpeedMotorA();    // Get the wheel speed of motor A (in radians/second)
  x[1] = x[0]; x[0] = va;         // Memorize the last two samples of the speed of motor A (in fact, a shift register)

  float k = readValue(0); // Read the value you set on QRoboticsCenter's channel 0
  writeValue(0, k);       // Send the value 0.01 to QRoboticsCenter's channel 0
  

}

bool Robot::controlEnabled() {
  return _button_states[0];       // The control is enabled if the state of button 0 is true
}

void Robot::button0callback() {
  if(toggleButton(0)) {           // Switches the state of button 0 and checks if the new state is true
    message("Robot enabled.");    // Display a message in the status bar of QRoboticsCenter
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
