#ifndef TRAJECTORY_H_
#define TRAJECTORY_H_

#include <fstream>
#include <vector>
#include <string>
#include <filesystem>

#include "Vector2.h"

namespace fs = std::filesystem;
struct Trajectory {
    std::vector<AVO::Vector2> positions;
    std::vector<AVO::Vector2> velocities;
    std::vector<AVO::Vector2> accelerations;
};

inline void saveTrajectoriesYAML(
    const std::vector<Trajectory>& trajectories,
    const std::string& filename)
  {
    // std::ofstream out(filename);
    std::ofstream out(filename.c_str());

    out << "result:\n";
    for (const auto& traj : trajectories) {

    out << "  -\n";
    out << "    num_states: "
        << traj.positions.size() << "\n";

    out << "    states:\n";

    for (size_t i = 0; i < traj.positions.size(); ++i) {

      const auto& pos = traj.positions[i];
      const auto& vel = traj.velocities[i];

      out << "      - ["
          << pos.x_ << ", "
          << pos.y_ << ", "
          << vel.x_ << ", "
          << vel.y_ << "]\n";
    }
    out << "    num_actions: "
    << traj.accelerations.size() << "\n";
    out << "    actions:\n";
    for (size_t i = 0; i < traj.accelerations.size(); ++i) {
      const auto& acc = traj.accelerations[i];
      out << "      - ["
          << acc.x_ << ", "
          << acc.y_ << "]\n";
    }
  }
}

// given sequence of states, compute accelerations
void computeAccelerations(std::vector<Trajectory>& trajectories, double dt)
{
    for (auto& traj : trajectories) {

        size_t N = traj.velocities.size();
        if (N < 2) {
            traj.accelerations.clear();
            continue;
        }

        traj.accelerations.resize(N - 1);

        for (size_t i = 0; i + 1 < N; ++i) {
            traj.accelerations[i].x_ =
                (traj.velocities[i + 1].x_ - traj.velocities[i].x_) / dt;

            traj.accelerations[i].y_ =
                (traj.velocities[i + 1].y_ - traj.velocities[i].y_) / dt;
        }
    }
}
// sanity check for the dynamics, inter-robot collisions
bool sanityCheck(
    const std::vector<Trajectory>& trajectories,
    double dt,
    double eps = 1e-6)
{
    constexpr double robot_radius = 0.25;
    constexpr double min_dist = 2.0 * robot_radius;

    // Check trajectory dimensions and dynamics
    for (const auto& traj : trajectories) {

      size_t N = traj.positions.size();

      if (traj.velocities.size() != N ||
          traj.accelerations.size() + 1 != N)
      {
          return false;
      }

      for (size_t i = 0; i + 1 < N; ++i) {

        const auto& p = traj.positions[i];
        const auto& p_next = traj.positions[i + 1];

        const auto& v = traj.velocities[i];
        const auto& v_next = traj.velocities[i + 1];
        const auto& a = traj.accelerations[i];

        // Semi-implicit Euler prediction
        auto v_pred = v + a * dt;
        auto p_pred = p + v_pred * dt;

        AVO::Vector2 vv = v_pred - v_next;
        AVO::Vector2 pp = p_pred - p_next;

        if (norm(vv) > eps ||
            norm(pp) > eps)
        {
            std::cout << "position error: " << norm(pp) << std::endl;
            std::cout << "velocity error: " << norm(vv) << std::endl;
            return false;
        }
      }
    }

    // Check inter-robot collisions
    if (trajectories.empty())
        return true;

    const size_t N = trajectories[0].positions.size();

    for (size_t t = 0; t < N; ++t) {

      for (size_t i = 0; i < trajectories.size(); ++i) {
        for (size_t j = i + 1; j < trajectories.size(); ++j) {

          // Make sure both trajectories have this timestep
          if (t >= trajectories[i].positions.size() ||
              t >= trajectories[j].positions.size())
          {
              return false;
          }

          const auto& pi = trajectories[i].positions[t];
          const auto& pj = trajectories[j].positions[t];

          double dist = norm(pi - pj);

          if (dist < min_dist - eps)
          {
            std::cout
                << "INTER-ROBOT COLLISION: "
                << "robots " << i << " and " << j
                << ", timestep " << t
                << ", distance = " << dist
                << ", required >= " << min_dist
                << std::endl;

            return false;
          }
        }
      }
    }

    return true;
}
// for stats
double compute_cost(const std::vector<Trajectory>& trajectories)
{
    double cost = 0.0;

    for (const auto& traj : trajectories) {
        cost += static_cast<double>(traj.accelerations.size());
    }

    return cost * 0.1;
}

size_t compute_makespan(const std::vector<Trajectory>& trajectories)
{
    size_t makespan = 0;

    for (const auto& traj : trajectories) {
        makespan = std::max(makespan, traj.accelerations.size());
    }

    return makespan * 0.1;
}

void create_dir_if_necessary(const std::string& file_path)
{

  fs::path path(file_path);
  if (path.has_parent_path()) {
      fs::create_directories(path.parent_path());
  }
}

#endif