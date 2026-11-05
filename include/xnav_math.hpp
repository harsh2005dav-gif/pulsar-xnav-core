#ifndef XNAV_MATH_HPP
#define XNAV_MATH_HPP

#include <cmath>
#include <string>

namespace xnav {

// Physical Constants (SI Units)
static constexpr double SPEED_OF_LIGHT = 299792458.0;       // m/s
static constexpr double G_CONST        = 6.67430e-11;       // m^3 kg^-1 s^-2
static constexpr double SUN_MASS       = 1.98847e30;        // kg
static constexpr double MU_SUN         = 1.32712440018e20;  // m^3 / s^2

struct Vector3D {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    double dot(const Vector3D& other) const {
        return (x * other.x) + (y * other.y) + (z * other.z);
    }

    double norm() const {
        return std::sqrt((x * x) + (y * y) + (z * z));
    }

    Vector3D normalized() const {
        double n = norm();
        if (n > 1e-12) {
            return Vector3D{x / n, y / n, z / n};
        }
        return Vector3D{0.0, 0.0, 0.0};
    }
};

struct PulsarBeacon {
    std::string name;
    Vector3D unit_direction;
    double period;
    double p_dot;
};

// 1. Romer Geometric Delay
inline double roemer_delay(const Vector3D& r_sc, const Vector3D& n_pulsar) {
    return r_sc.dot(n_pulsar) / SPEED_OF_LIGHT;
}

// 2. Relativistic Shapiro Delay around Sun
inline double shapiro_delay(const Vector3D& r_sc, const Vector3D& n_pulsar) {
    double r_mag = r_sc.norm();
    double cos_theta = r_sc.dot(n_pulsar) / r_mag;
    double arg = 1.0 + cos_theta;
    if (arg <= 1e-12) {
        arg = 1e-12;
    }
    return -(2.0 * G_CONST * SUN_MASS / std::pow(SPEED_OF_LIGHT, 3)) * std::log(arg);
}

// 3. Pulse Phase Calculation
inline double compute_pulse_phase(double t_sc, const Vector3D& r_sc, const PulsarBeacon& pulsar) {
    double dt_r = roemer_delay(r_sc, pulsar.unit_direction);
    double dt_s = shapiro_delay(r_sc, pulsar.unit_direction);
    double t_ssb = t_sc + dt_r + dt_s;
    return std::fmod(t_ssb / pulsar.period, 1.0);
}

} // namespace xnav

#endif