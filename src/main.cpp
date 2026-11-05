#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include "../include/xnav_math.hpp"

using namespace xnav;

int main() {
    std::cout << "=========================================================\n";
    std::cout << "   PULSAR-XNAV-CORE: Deep-Space Relativistic Navigator   \n";
    std::cout << "=========================================================\n";

    // High-stability X-ray millisecond pulsars
    std::vector<PulsarBeacon> beacons = {
        {"PSR B1937+21",  Vector3D{0.183, 0.892, 0.413}.normalized(),  0.0015578, 1.05e-19},
        {"PSR B1821-24",  Vector3D{-0.245, -0.781, -0.574}.normalized(), 0.0030543, 1.61e-18},
        {"PSR J0437-4715",Vector3D{0.412, -0.321, -0.852}.normalized(), 0.0057574, 5.73e-20},
        {"Crab B0531+21", Vector3D{-0.851, 0.320, 0.417}.normalized(),  0.0333924, 4.21e-13}
    };

    // True Initial Spacecraft State (Earth-Mars transfer trajectory)
    Vector3D true_pos{1.496e11, 2.5e10, 1.0e10}; 
    Vector3D true_vel{-5000.0, 29000.0, 1200.0}; 

    // Initial Estimated State: 50 km offset
    Vector3D est_pos{true_pos.x + 35000.0, true_pos.y - 25000.0, true_pos.z + 15000.0};
    Vector3D est_vel{true_vel.x, true_vel.y, true_vel.z};

    double p_pos = 2.5e9;    // Initial variance (~50 km 1-sigma)
    double r_noise_std = 1e-9; // Sub-nanosecond photon arrival precision

    std::ofstream out_file("telemetry.csv");
    out_file << "step,true_x,true_y,true_z,est_x,est_y,est_z,pos_error_km,sigma_km\n";

    double dt = 60.0; // 60 seconds per epoch
    std::mt19937_64 rng(42);
    std::normal_distribution<double> sensor_noise(0.0, r_noise_std);

    std::cout << "[SIM] Filtering relativistic pulsar TOA observables...\n";

    for (int step = 0; step < 200; ++step) {
        double t_current = step * dt;

        // Dynamics Propagation (Sun gravity)
        double r3_true = std::pow(true_pos.norm(), 3);
        true_vel.x -= (MU_SUN * true_pos.x / r3_true) * dt;
        true_vel.y -= (MU_SUN * true_pos.y / r3_true) * dt;
        true_vel.z -= (MU_SUN * true_pos.z / r3_true) * dt;
        true_pos.x += true_vel.x * dt;
        true_pos.y += true_vel.y * dt;
        true_pos.z += true_vel.z * dt;

        double r3_est = std::pow(est_pos.norm(), 3);
        est_vel.x -= (MU_SUN * est_pos.x / r3_est) * dt;
        est_vel.y -= (MU_SUN * est_pos.y / r3_est) * dt;
        est_vel.z -= (MU_SUN * est_pos.z / r3_est) * dt;
        est_pos.x += est_vel.x * dt;
        est_pos.y += est_vel.y * dt;
        est_pos.z += est_vel.z * dt;

        // Process noise growth
        p_pos += 25.0 * dt;

        // Ingest pulsar observations
        for (const auto& beacon : beacons) {
            // Relativistic delay difference
            double true_toa = roemer_delay(true_pos, beacon.unit_direction) + 
                              shapiro_delay(true_pos, beacon.unit_direction) + sensor_noise(rng);
            double pred_toa = roemer_delay(est_pos, beacon.unit_direction) + 
                              shapiro_delay(est_pos, beacon.unit_direction);
            
            double delta_t = true_toa - pred_toa;

            // Jacobian component along sightline H = n_i / c
            double h_scale = 1.0 / SPEED_OF_LIGHT;
            double r_cov = r_noise_std * r_noise_std;
            
            // Kalman gain
            double s_innov = (h_scale * h_scale) * p_pos + r_cov;
            double k = (p_pos * h_scale) / s_innov;

            // Position correction
            est_pos.x += k * beacon.unit_direction.x * delta_t;
            est_pos.y += k * beacon.unit_direction.y * delta_t;
            est_pos.z += k * beacon.unit_direction.z * delta_t;

            // Covariance reduction
            p_pos = (1.0 - k * h_scale) * p_pos;
        }

        double dx = est_pos.x - true_pos.x;
        double dy = est_pos.y - true_pos.y;
        double dz = est_pos.z - true_pos.z;
        double error_km = std::sqrt(dx*dx + dy*dy + dz*dz) / 1000.0;
        double sigma_km = std::sqrt(p_pos) / 1000.0;

        out_file << step << "," << true_pos.x << "," << true_pos.y << "," << true_pos.z << ","
                 << est_pos.x << "," << est_pos.y << "," << est_pos.z << ","
                 << error_km << "," << sigma_km << "\n";
    }

    out_file.close();
    std::cout << "[SUCCESS] Convergence telemetry logged to 'telemetry.csv'.\n";
    return 0;
}