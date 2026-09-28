// Mars lander simulator
// Version 1.11
// Mechanical simulation functions
// Gabor Csanyi and Andrew Gee, August 2019

// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation, to make use of it
// for non-commercial purposes, provided that (a) its original authorship
// is acknowledged and (b) no modified versions of the source code are
// published. Restriction (b) is designed to protect the integrity of the
// exercise for future generations of students. The authors would be happy
// to receive any suggested modifications by private correspondence to
// ahg@eng.cam.ac.uk and gc121@eng.cam.ac.uk.

#include "lander.h"
#include <cmath>
#include <vector>
using namespace std;

void autopilot (void)
  // Autopilot to adjust the engine throttle, parachute and attitude control
{//define variables
double Kh = 0.020;
double Kp = 0.2;
double delta = 0.35;
static ofstream autopilot_log;
static bool autopilot_log_initialized = false;
vector3d e_r = position.norm();
double h = position.abs() - MARS_RADIUS;
double v_r = velocity * e_r;//vector3d dot product re-defined.
double v_target = -(0.5 + Kh * h);
double e = -(0.5 + Kh * h + v_r);
double P_out = Kp * e;

// Minimal-intrusion logging for Assignment 5 analysis in Python
if (!autopilot_log_initialized || simulation_time <= 0.0) {
  if (autopilot_log.is_open()) autopilot_log.close();
  autopilot_log.open("autopilot_descent_log.txt");
  if (autopilot_log) {
    autopilot_log << "time_s altitude_m v_r_mps v_target_mps throttle\n";
  }
  autopilot_log_initialized = true;
}

//强制简化条件：
stabilized_attitude = true;
parachute_status = NOT_DEPLOYED;

//define throttle function
  if (P_out <= -delta) {
    throttle = 0;
  } else if (P_out >= 1-delta) {
    throttle = 1;
  } else {
    throttle = delta + P_out;
  }
  if (throttle < 0.0) throttle = 0.0;
  if (throttle > 1.0) throttle = 1.0;

  if (autopilot_log && (scenario == 1 || scenario == 5)) {
    autopilot_log << simulation_time << " "
                  << h << " "
                  << v_r << " "
                  << v_target << " "
                  << throttle << "\n";
  }
}


void numerical_dynamics (void)
  // This is the function that performs the numerical integration to update the
  // lander's pose. The time step is delta_t (global variable).
{
  const bool use_verlet = true; // true: Verlet, false: Euler.const 说明本次函数调用中它不会被改。
  static bool verlet_initialized = false;//static 可以记住上一次调用时的值，下一次调用时会继续使用上一次的值。
  //第一次调用之后就变成true，之后每次都是true。
  static vector3d previous_position;
  //declare variables
  double mass = UNLOADED_LANDER_MASS + fuel * FUEL_CAPACITY * FUEL_DENSITY;
  double r = position.abs();
  double v_abs = velocity.abs();
  double rho = atmospheric_density(position);
  vector3d r_hat = position.norm();
  double A_lander = M_PI * LANDER_SIZE * LANDER_SIZE;
  //SMALL_NUM 是很小的阈值（接近 0 的值），用于避免导数不稳定
  vector3d v_hat = (v_abs > SMALL_NUM) ? velocity.norm() : vector3d(0.0, 0.0, 0.0);
  vector3d a_gravity = -(GRAVITY * MARS_MASS / (r * r)) * r_hat;
  vector3d a_thrust = thrust_wrt_world() / mass;
  vector3d F_drag_lander = -0.5 * rho * DRAG_COEF_LANDER  * A_lander * velocity.abs2() * v_hat;
  vector3d F_drag_chute(0.0, 0.0, 0.0);
  if (parachute_status == DEPLOYED) {
    double A_chute = 5.0 * (2.0 * LANDER_SIZE) * (2.0 * LANDER_SIZE);
    F_drag_chute = -0.5 * rho * DRAG_COEF_CHUTE  * A_chute * velocity.abs2() * v_hat;
  }
  vector3d a_drag = (F_drag_lander + F_drag_chute) / mass;
  //total acceleration
  vector3d a_total = a_gravity + a_thrust + a_drag;

  // Reset Verlet memory on simulation/scenario restart
  if (simulation_time <= 0.0) {
    verlet_initialized = false;
  }

  //Euler or Verlet:
  if (use_verlet) {//Verlet:
    if (!verlet_initialized) {
      // 2nd-order back-step initialization:
      // x(t-dt) = x(t) - v(t)dt + 0.5*a(t)*dt^2
      previous_position = position - velocity * delta_t + 0.5 * a_total * delta_t * delta_t;
      verlet_initialized = true;
    }
    vector3d next_position = 2.0 * position - previous_position + a_total * delta_t * delta_t;
    vector3d next_velocity = (next_position - previous_position) / (2.0 * delta_t);
    previous_position = position;
    position = next_position;
    velocity = next_velocity;
  } else {//Euler:
  position = position +velocity * delta_t;
  velocity = velocity + a_total * delta_t;
  }

  // Here we can apply an autopilot to adjust the thrust, parachute and attitude
  if (autopilot_enabled) autopilot();

  // Here we can apply 3-axis stabilization to ensure the base is always pointing downwards
  if (stabilized_attitude) attitude_stabilization();
}

void initialize_simulation (void)
  // Lander pose initialization - selects one of 10 possible scenarios
{
  // The parameters to set are:
  // position - in Cartesian planetary coordinate system (m)
  // velocity - in Cartesian planetary coordinate system (m/s)
  // orientation - in lander coordinate system (xyz Euler angles, degrees)
  // delta_t - the simulation time step
  // boolean state variables - parachute_status, stabilized_attitude, autopilot_enabled
  // scenario_description - a descriptive string for the help screen

  scenario_description[0] = "circular orbit";
  scenario_description[1] = "descent from 10km";
  scenario_description[2] = "elliptical orbit, thrust changes orbital plane";
  scenario_description[3] = "polar launch at escape velocity (but drag prevents escape)";
  scenario_description[4] = "elliptical orbit that clips the atmosphere and decays";
  scenario_description[5] = "descent from 200km";
  scenario_description[6] = "areostationary orbit";
  scenario_description[7] = "";
  scenario_description[8] = "";
  scenario_description[9] = "";

  switch (scenario) {

  case 0:
    // a circular equatorial orbit
    position = vector3d(1.2*MARS_RADIUS, 0.0, 0.0);
    velocity = vector3d(0.0, -3247.087385863725, 0.0);
    orientation = vector3d(0.0, 90.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 1:
    // a descent from rest at 10km altitude
    position = vector3d(0.0, -(MARS_RADIUS + 10000.0), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = false;
    break;

  case 2:
    // an elliptical polar orbit
    position = vector3d(0.0, 0.0, 1.2*MARS_RADIUS);
    velocity = vector3d(3500.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 3:
    // polar surface launch at escape velocity (but drag prevents escape)
    position = vector3d(0.0, 0.0, MARS_RADIUS + LANDER_SIZE/2.0);
    velocity = vector3d(0.0, 0.0, 5027.0);
    orientation = vector3d(0.0, 0.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 4:
    // an elliptical orbit that clips the atmosphere each time round, losing energy
    position = vector3d(0.0, 0.0, MARS_RADIUS + 100000.0);
    velocity = vector3d(4000.0, 0.0, 0.0);
    orientation = vector3d(0.0, 90.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 5:
    // a descent from rest at the edge of the exosphere
    position = vector3d(0.0, -(MARS_RADIUS + EXOSPHERE), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = false;
    break;

  case 6:
    // areostationary orbit
    {
      double T_stationary = MARS_DAY; // orbital period (s)
      double r_stationary = cbrt((GRAVITY * MARS_MASS * T_stationary * T_stationary) /
                                 (4.0 * M_PI * M_PI));
      double v_stationary = sqrt((GRAVITY * MARS_MASS) / r_stationary);
      // Equatorial orbit in x-y plane
      position = vector3d(r_stationary, 0.0, 0.0);
      velocity = vector3d(0.0, -v_stationary, 0.0);
    }
    orientation = vector3d(0.0, 90.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 7:
    break;

  case 8:
    break;

  case 9:
    break;

  }
}
//Run:cd "C:\Users\13410\Desktop\Mars lander\lander"
//& "C:\msys64\mingw64\bin\g++.exe" lander.cpp lander_graphics.cpp -O2 -o lander.exe -lfreeglut -lglu32 -lopengl32
//.\lander.exe

//修改完之后要重新编译！