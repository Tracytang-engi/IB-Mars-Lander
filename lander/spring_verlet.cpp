#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main() {

  // declare variables
  double m, k, x0, v0, t_max, dt, t, a_curr, x_prev, x_curr, x_next, v_curr;
  int n_steps;
  vector<double> t_list, x_list, v_list;

  // mass, spring constant, initial position and velocity
  m = 1;
  k = 1;
  x0 = 0;
  v0 = 1;

  // simulation time and timestep
  t_max = 100;
  dt = 0.1;

  // Verlet initialization
  x_prev = x0;                                // x(0)
  a_curr = -k * x_prev / m;                  // a(0)
  x_curr = x_prev + dt * v0 + 0.5 * dt * dt * a_curr; // x(dt)

  // store t=0 state
  t_list.push_back(0.0);
  x_list.push_back(x_prev);
  v_list.push_back(v0);

  // position-Verlet integration
  for (t = dt; t<=t_max; t+=dt) {
    a_curr = -k * x_curr / m;
    x_next = 2.0 * x_curr - x_prev + dt * dt * a_curr;
    v_curr = (x_next - x_prev) / (2.0 * dt);

    t_list.push_back(t);
    x_list.push_back(x_curr);
    v_list.push_back(v_curr);

    x_prev = x_curr;
    x_curr = x_next;
  }

  // Write the trajectories to file
  ofstream fout;
  fout.open("trajectories_verlet.txt");
  if (fout) {
    for (int i = 0; i < t_list.size(); i = i + 1) {
      fout << t_list[i] << ' ' << x_list[i] << ' ' << v_list[i] << endl;
    }
  } else {
    cout << "Could not open trajectory file for writing" << endl;
  }
}
