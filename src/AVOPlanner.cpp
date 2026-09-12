#include <cmath>
#include <cstddef>
#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>
#include <yaml-cpp/yaml.h>
#include <chrono>
#include <AVO.h>
#include <Trajectory.h>

static std::vector<AVO::Vector2> goals; // velocity always 0
// later use with input.yaml
void setupScenario(AVO::Simulator *sim)
{
	/* Specify the global time step of the simulation. */
	sim->setTimeStep(0.1f);
    // neighDist, maxNeighb, timeHorizon, radius, maxSpeed, maxAccel, accelInt.
    sim->setAgentDefaults(15.0F, 10U, 10.0F, 0.5, 5.0F, 2.0F, 1.0F);
}

void setPreferredVelocities(AVO::Simulator* sim)
{
    for (size_t i = 0; i < sim->getNumAgents(); ++i) {
        AVO::Vector2 goalVector = goals[i] - sim->getAgentPosition(i);

        if (AVO::absSq(goalVector) > 1.0F) {
            goalVector = normalize(goalVector); // velocity (vx,vy) towards the goal
        }
        if(i==5)
            std::cout << "your robot" << std::endl;
        sim->setAgentPrefVelocity(i, goalVector);
    }
}

bool haveReachedGoals(AVO::Simulator* sim)
{
    double goal_threshold = 0.25;
    for (size_t i = 0; i < sim->getNumAgents(); ++i) {
        if (AVO::absSq(sim->getAgentPosition(i) - goals[i]) > goal_threshold) {
            return false;
        }
    }
    return true;
}

void runAVO(const std::string& input_yaml,
            const std::string& output_yaml,
            const std::string& stats_yaml)
{
    auto start = std::chrono::steady_clock::now();
	  AVO::Simulator *sim = new AVO::Simulator();
    setupScenario(sim);
    YAML::Node config = YAML::LoadFile(input_yaml);
    auto robots = config["robots"];

    if (!robots) {
      throw std::runtime_error("Missing 'robots' field in YAML");
    }

    create_dir_if_necessary(stats_yaml);
    std::ofstream stats(stats_yaml, std::ios::app);
    if (!stats)
    {
      std::cerr << "Failed to open stats.yaml file.\n";
      return;
    }

    std::vector<Trajectory> trajectories;

    goals.clear();
    trajectories.clear();

    // init. robots
    for (size_t i = 0; i < robots.size(); ++i) {
        auto robot = robots[i];

        auto start = robot["start"];
        auto goal  = robot["goal"];

        AVO::Vector2 start_pos(
            start[0].as<float>(),
            start[1].as<float>()
        );

        AVO::Vector2 goal_pos(
            goal[0].as<float>(),
            goal[1].as<float>()
        );

        sim->addAgent(start_pos);
        goals.push_back(goal_pos);
        trajectories.emplace_back();
    }
    // motion planning
    do {
        for (size_t i = 0; i < sim->getNumAgents(); ++i) {
          trajectories[i].positions.push_back(sim->getAgentPosition(i));
          trajectories[i].velocities.push_back(sim->getAgentVelocity(i));
        }
        setPreferredVelocities(sim);
        sim->doStep();

    } while (!haveReachedGoals(sim));
    // compute accelerations
    computeAccelerations(trajectories, /*dt*/sim->getTimeStep());
    if(!sanityCheck(trajectories, /*dt*/sim->getTimeStep())){
      std::cout << "Dynamics are violated" << std::endl;
      return;
    }
    auto end = std::chrono::steady_clock::now();
    double elapsed_sec =
      std::chrono::duration<double>(end - start).count();
    std::cout << "Time: " << elapsed_sec << " sec\n";
    // save the output
    saveTrajectoriesYAML(trajectories, output_yaml);
    std::cout << "AVO simulation finished with "
              << sim->getNumAgents()
              << " agents." << std::endl;
    // save stats
    double cost = compute_cost(trajectories);
    double makespan = compute_makespan(trajectories);
    // save stats
    stats << "stats: " << "\n";
    stats << "  - t: " << elapsed_sec << "\n";
    stats << "    cost: " << cost << "\n";
    stats << "    makespan: " << makespan << "\n";
    stats.flush();
}

int main(int argc, char** argv)
{
    if (argc < 4) {
        std::cerr << "Usage: AVOPlanner <input.yaml> <output.yaml> <stats.yaml>\n";
        return 1;
    }

    try {
        runAVO(argv[1], argv[2], argv[3]);
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}