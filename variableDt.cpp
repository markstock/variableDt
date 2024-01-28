//
// variableDtAB.cpp
//
// test code for variable time step size Adams-Bashforth
//
// (c)2024 Mark J. Stock <markjstock@gmail.com>
//

#include <iostream>
#include <iomanip>
#include <array>
#include <cmath>
#include <numbers>
#include <random>

// time is always double-precision, solution and vel are templated
template<typename T>
struct CosineSystem {
  double tstart = 0.0;
  double tend = std::numbers::pi_v<double>;
  T solution_at_time(const double _time) {
    return std::sin(_time);
  }
  T vel_at_time(const double _time) {
    return std::cos(_time);
  }
};

// times are always in double, computed parts use template type
template<typename T>
T run_ab2 (CosineSystem<T>& _system, const double _basedt,
           const bool _is_rand_dt, std::mt19937& _gen, const bool _silent) {

  // each step needs a dt, value (x), and derivative
  // dt[i] is the time step size for the step just before x[i]
  std::array<T,2> dt, x, dxdt;
  double init_time = _system.tstart;
  for (int32_t i=1; i>-1; --i) {
    x[i]    = _system.solution_at_time(init_time);
    dxdt[i] = _system.vel_at_time(init_time);
    dt[i]   = _basedt;
    init_time -= dt[i];
  }

  std::uniform_real_distribution<T> onetotwo(1.0, 2.0);

  int32_t istep = 0;
  double time = _system.tstart;
  while (time < _system.tend - 1.e-10) {

    // determine "velocity" at this time
    dxdt[1] = _system.vel_at_time(time);

    // take a new step forward with this dt
    T this_dt = _basedt;
    if (_is_rand_dt) this_dt *= onetotwo(_gen);
    if (time + this_dt > _system.tend) this_dt = _system.tend - time;
    // constant-dt AB2:
    //const T new_x = x[1] + 0.5*this_dt*(3.0*dxdt[1]-dxdt[0]);
    // variable-dt AB2:
    const T dtratio = this_dt / dt[1];
    const T new_x = x[1] + (T)0.5*this_dt*(((T)2.0+dtratio)*dxdt[1]-dtratio*dxdt[0]);
    time += this_dt;

    // update the arrays
    dt[0] = dt[1];	dt[1] = this_dt;
    x[0] = x[1];	x[1] = new_x;
    dxdt[0] = dxdt[1];
    ++istep;

    // debug
    //std::cout << istep << "\t" << std::setw(15) << std::setprecision(12) << time << "\t" << x[1] << "\n";
  }
  if (not _silent) std::cout << istep << "\t" << std::setw(15) << std::setprecision(12) << time << "\t" << x[1] << "\n";

  // return absolute error
  return std::abs(x[1] - _system.solution_at_time(_system.tend));
}

template<typename T>
T run_ab3 (CosineSystem<T>& _system, const double _basedt,
           const bool _is_rand_dt, std::mt19937& _gen, const bool _silent) {

  // each step needs a dt, value (x), and derivative
  // dt[i] is the time step size for the step just before x[i]
  std::array<T,3> dt, x, dxdt;
  double init_time = _system.tstart;
  for (int32_t i=2; i>-1; --i) {
    x[i]    = _system.solution_at_time(init_time);
    dxdt[i] = _system.vel_at_time(init_time);
    dt[i]   = _basedt;
    init_time -= dt[i];
  }

  std::uniform_real_distribution<T> onetotwo(1.0, 2.0);

  int32_t istep = 0;
  double time = _system.tstart;
  while (time < _system.tend - 1.e-10) {

    // determine "velocity" at this time
    dxdt[2] = _system.vel_at_time(time);

    // take a new step forward with this dt
    T this_dt = _basedt;
    if (_is_rand_dt) this_dt *= onetotwo(_gen);
    if (time + this_dt > _system.tend) this_dt = _system.tend - time;

    // constant-dt AB3:
    //const T cfactor = ((T)23.*dxdt[2]-(T)16.*dxdt[1]+(T)5.*dxdt[0]) / (T)12.;
    // variable-dt AB3:
    const T dtr1 = this_dt / dt[2];							// 1
    const T dtr2 = this_dt / (this_dt+dt[2]);				// 0.5
    const T dtr3 = dtr1 * (this_dt+dt[2]) / (dt[2]+dt[1]);	// 1
    const T dtr4 = dt[2] / dt[1];							// 1
    const T dtr5 = (T)0.5 * ((T)1.-dtr2/(T)3.);
    const T factor = dxdt[2]
                   + (T)0.5 * dtr1 * (dxdt[2]-dxdt[1]) 
                   + dtr5 * dtr3 * (dxdt[2] - dxdt[1] - dtr4*(dxdt[1]-dxdt[0]));
    //std::cout << istep << "\t" << std::setw(15) << std::setprecision(12) << cfactor << "\t" << factor << "\n";
    const T new_x = x[2] + this_dt*factor;
    time += this_dt;

    // update the arrays
    dt[0] = dt[1];	dt[1] = dt[2];  dt[2] = this_dt;
     x[0] = x[1];	 x[1] = x[2];    x[2] = new_x;
    dxdt[0] = dxdt[1]; dxdt[1] = dxdt[2];
    ++istep;

    // debug
    //std::cout << istep << "\t" << std::setw(15) << std::setprecision(12) << time << "\t" << x[2] << "\n";
  }
  if (not _silent) std::cout << istep << "\t" << std::setw(15) << std::setprecision(12) << time << "\t" << x[2] << "\n";

  return std::abs(x[2] - _system.solution_at_time(_system.tend));
}

template<typename T>
T run_ab4 (CosineSystem<T>& _system, const double _basedt,
           const bool _is_rand_dt, std::mt19937& _gen, const bool _silent) {

  // each step needs a dt, value (x), and derivative
  // dt[i] is the time step size for the step just before x[i]
  std::array<T,4> dt, x, dxdt;
  double init_time = 0.0;
  for (int32_t i=3; i>-1; --i) {
    x[i]    = _system.solution_at_time(init_time);
    dxdt[i] = _system.vel_at_time(init_time);
    dt[i]   = _basedt;
    init_time -= dt[i];
  }

  std::uniform_real_distribution<T> onetotwo(1.0, 2.0);

  int32_t istep = 0;
  double time = _system.tstart;

  while (time < _system.tend - 1.e-10) {

    // determine "velocity" at this time
    dxdt[3] = _system.vel_at_time(time);

    // take a new step forward with this dt
    T this_dt = _basedt;
    if (_is_rand_dt) this_dt *= onetotwo(_gen);
    if (time + this_dt > _system.tend) this_dt = _system.tend - time;

    // constant-dt AB4 (8 flops):
    const double cfactor = ((T)55.*dxdt[3]-(T)59.*dxdt[2]+(T)37.*dxdt[1]-(T)9.*dxdt[0]) / (T)24.0;
    // variable-dt AB4 (49 flops):
    const T dtr1 = this_dt / dt[3];							// 1
    const T dtr2 = this_dt / (this_dt+dt[3]);				// 0.5
    const T dtr3 = dtr1 * (this_dt+dt[3]) / (dt[3]+dt[2]);	// 1
    const T dtr4 = dt[3] / dt[2];							// 1
    const T dtr5 = (T)0.5 * ((T)1.-dtr2/(T)3.);				// 5/12
    const T acc1 = dxdt[3] - dxdt[2];						// 3 - 2
    const T acc2 = acc1 - dtr4*(dxdt[2]-dxdt[1]);			// 3 - 2*2 + 1
    const T factor = dxdt[3]
                   + (T)0.5 * dtr1 * acc1
                   + dtr5 * dtr3 * acc2;					// 3 + 0.5*3 - 0.5*2 + (5/12)*(3 - 2*2 + 1)
															// (47/12)*3 - (16/12)*2 + (5/12)*1
    const T dtr6 = ((T)1.-dtr2/(T)2.) / (T)6.;				// 3/24
    const T dtr7 = this_dt / (this_dt+dt[3]+dt[2]);			// 1/3
    const T dtr8 = dtr3 * (this_dt+dt[3]+dt[2]) / (dt[3]+dt[2]+dt[1]);	// 1
    const T dtr9 = (dt[3]+dt[2]) / (dt[2]+dt[1]);			// 1
    const T acc3 = acc2 - dtr4*dtr9*(dxdt[2] - dxdt[1] - (dt[2]/dt[1])*(dxdt[1]-dxdt[0]));	// 3 - 3*2 + 3*1 - 0
    const T dtr0 = (dtr5 - dtr6*dtr7) * dtr8;				// 5/12 - 3/72 = 27/72
    const T factor2 = factor + dtr0 * acc3;					// 
    //std::cout << istep << "\t" << std::setw(15) << std::setprecision(12) << cfactor << "\t" << factor2 << "\n";
    const double new_x = x[3] + this_dt*factor2;
    time += this_dt;

    // update the arrays
    dt[0] = dt[1];	dt[1] = dt[2];  dt[2] = dt[3];  dt[3] = this_dt;
     x[0] = x[1];	 x[1] = x[2];    x[2] = x[3];  x[3] = new_x;
    dxdt[0] = dxdt[1]; dxdt[1] = dxdt[2]; dxdt[2] = dxdt[3];
    ++istep;

    // debug
    //std::cout << istep << "\t" << std::setw(15) << std::setprecision(12) << time << "\t" << x[3] << "\n";
  }
  if (not _silent) std::cout << istep << "\t" << std::setw(15) << std::setprecision(12) << time << "\t" << x[3] << "\n";

  return std::abs(x[3] - _system.solution_at_time(_system.tend));
}

int main(int argc, char *argv[]) {

  // random number stuff
  std::mt19937 gen(12345);

  // a simple dynamic system
  CosineSystem<double> sys;

  // uniform dt

  std::cout << "Constant dt, AB2:\n";
  double basedt = sys.tend / 10.0;
  double lasterr = 0.0;
  for (int32_t isim=0; isim<6; ++isim) {
    double thiserr = run_ab2<double>(sys, basedt, false, gen, true);
    std::cout << std::setw(12) << std::setprecision(8) << int32_t(0.5+sys.tend/basedt) << "\t" << thiserr;
    if (isim != 0) std::cout << "\t" << lasterr/thiserr;
    std::cout << "\n";
    lasterr = thiserr;
    basedt *= 0.5;
  }

  std::cout << "Constant dt, AB3:\n";
  basedt = sys.tend / 10.0;
  lasterr = 0.0;
  for (int32_t isim=0; isim<6; ++isim) {
    double thiserr = run_ab3<double>(sys, basedt, false, gen, true);
    std::cout << std::setw(12) << std::setprecision(8) << int32_t(0.5+sys.tend/basedt) << "\t" << thiserr;
    if (isim != 0) std::cout << "\t" << lasterr/thiserr;
    std::cout << "\n";
    lasterr = thiserr;
    basedt *= 0.5;
  }

  std::cout << "Constant dt, AB4:\n";
  basedt = sys.tend / 10.0;
  lasterr = 0.0;
  for (int32_t isim=0; isim<6; ++isim) {
    double thiserr = run_ab4<double>(sys, basedt, false, gen, true);
    std::cout << std::setw(12) << std::setprecision(8) << int32_t(0.5+sys.tend/basedt) << "\t" << thiserr;
    if (isim != 0) std::cout << "\t" << lasterr/thiserr;
    std::cout << "\n";
    lasterr = thiserr;
    basedt *= 0.5;
  }

  // variable dt
  const int32_t maxsims = 1000;

  std::cout << "Variable dt, AB2:\n";
  basedt = sys.tend / 10.0;
  lasterr = 0.0;
  for (int32_t isim=0; isim<6; ++isim) {
    double errsum = 0.0;
    for (int32_t istep=0; istep<maxsims; ++istep) {
      errsum += run_ab2<double>(sys, basedt, true, gen, true);
    }
    double thiserr = errsum/maxsims;
    std::cout << std::setw(12) << std::setprecision(8) << int32_t(0.5+sys.tend/basedt) << "-" << int32_t(0.5+2.*sys.tend/basedt) << "  \t" << thiserr;
    if (isim != 0) std::cout << "\t" << lasterr/thiserr;
    std::cout << "\n";
    lasterr = thiserr;
    basedt *= 0.5;
  }

  std::cout << "Variable dt, AB3:\n";
  basedt = sys.tend / 10.0;
  lasterr = 0.0;
  for (int32_t isim=0; isim<6; ++isim) {
    double errsum = 0.0;
    for (int32_t istep=0; istep<maxsims; ++istep) {
      errsum += run_ab3<double>(sys, basedt, true, gen, true);
    }
    double thiserr = errsum/maxsims;
    std::cout << std::setw(12) << std::setprecision(8) << int32_t(0.5+sys.tend/basedt) << "-" << int32_t(0.5+2.*sys.tend/basedt) << "  \t" << thiserr;
    if (isim != 0) std::cout << "\t" << lasterr/thiserr;
    std::cout << "\n";
    lasterr = thiserr;
    basedt *= 0.5;
  }

  std::cout << "Variable dt, AB4:\n";
  basedt = sys.tend / 10.0;
  lasterr = 0.0;
  for (int32_t isim=0; isim<6; ++isim) {
    double errsum = 0.0;
    for (int32_t istep=0; istep<maxsims; ++istep) {
      errsum += run_ab4<double>(sys, basedt, true, gen, true);
    }
    double thiserr = errsum/maxsims;
    std::cout << std::setw(12) << std::setprecision(8) << int32_t(0.5+sys.tend/basedt) << "-" << int32_t(0.5+2.*sys.tend/basedt) << "  \t" << thiserr;
    if (isim != 0) std::cout << "\t" << lasterr/thiserr;
    std::cout << "\n";
    lasterr = thiserr;
    basedt *= 0.5;
  }
}
