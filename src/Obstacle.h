#ifndef AVO_OBSTACLE_H_
#define AVO_OBSTACLE_H_

#include "Vector2.h"

namespace AVO {
	/**
	 * \brief      Defines static obstacles in the simulation.
	 */
	class Obstacle {
	private:
		/**
		 * \brief      Constructs a static obstacle instance.
		 */
		Obstacle();

		bool isConvex_;
		Obstacle *nextObstacle_;
		Vector2 point_;
		Obstacle *prevObstacle_;
		Vector2 unitDir_;

		size_t id_;

		friend class Agent;
		friend class KdTree;
		friend class Simulator;
	};
}

#endif