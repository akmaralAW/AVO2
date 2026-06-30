#include "Obstacle.h"
#include "Simulator.h"

namespace AVO {
	Obstacle::Obstacle() : isConvex_(false), nextObstacle_(NULL), prevObstacle_(NULL), id_(0) { }
}