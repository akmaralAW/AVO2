/*
 * KdTree.cc
 * AVO2 Library
 *
 * SPDX-FileCopyrightText: 2010 University of North Carolina at Chapel Hill
 * SPDX-License-Identifier: Apache-2.0
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Please send all bug reports to <geom@cs.unc.edu>.
 *
 * The authors may be contacted via:
 *
 * Jur van den Berg, Jamie Snape, Stephen J. Guy, and Dinesh Manocha
 * Dept. of Computer Science
 * 201 S. Columbia St.
 * Frederick P. Brooks, Jr. Computer Science Bldg.
 * Chapel Hill, N.C. 27599-3175
 * United States of America
 *
 * <https://gamma.cs.unc.edu/AVO/>
 */

#include "KdTree.h"

#include <algorithm>
#include <utility>

#include "Agent.h"
#include "Simulator.h"
#include "Vector2.h"

namespace AVO {
namespace {
const std::size_t AVO_MAX_LEAF_SIZE = 10U;
}  // namespace

class KdTree::AgentTreeNode {
 public:
  AgentTreeNode();

  std::size_t begin;
  std::size_t end;
  std::size_t left;
  std::size_t right;
  float maxX;
  float maxY;
  float minX;
  float minY;
};

KdTree::AgentTreeNode::AgentTreeNode()
    : begin(0U),
      end(0U),
      left(0U),
      right(0U),
      maxX(0.0F),
      maxY(0.0F),
      minX(0.0F),
      minY(0.0F) {}


KdTree::KdTree(Simulator *simulator) : simulator_(simulator) {}

KdTree::~KdTree() {}

void KdTree::buildAgentTree() {
  if (agents_.size() < simulator_->agents_.size()) {
    for (std::size_t agentNo = agents_.size();
         agentNo < simulator_->agents_.size(); ++agentNo) {
      agents_.push_back(simulator_->agents_[agentNo]);
    }

    agentTree_.resize(2U * agents_.size() - 1U);
  }

  if (!agents_.empty()) {
    buildAgentTreeRecursive(0U, agents_.size(), 0U);
  }
}

void KdTree::buildAgentTreeRecursive(std::size_t begin, std::size_t end,
                                     std::size_t node) {
  agentTree_[node].begin = begin;
  agentTree_[node].end = end;
  agentTree_[node].minX = agentTree_[node].maxX = agents_[begin]->position_.x_;
  agentTree_[node].minY = agentTree_[node].maxY = agents_[begin]->position_.y_;

  for (std::size_t i = begin + 1U; i < end; ++i) {
    agentTree_[node].maxX =
        std::max(agentTree_[node].maxX, agents_[i]->position_.x_);
    agentTree_[node].minX =
        std::min(agentTree_[node].minX, agents_[i]->position_.x_);
    agentTree_[node].maxY =
        std::max(agentTree_[node].maxY, agents_[i]->position_.y_);
    agentTree_[node].minY =
        std::min(agentTree_[node].minY, agents_[i]->position_.y_);
  }

  if (end - begin > AVO_MAX_LEAF_SIZE) {
    // No leaf node.
    const bool isVertical = agentTree_[node].maxX - agentTree_[node].minX >
                            agentTree_[node].maxY - agentTree_[node].minY;
    const float splitValue =
        isVertical ? 0.5F * (agentTree_[node].maxX + agentTree_[node].minX)
                   : 0.5F * (agentTree_[node].maxY + agentTree_[node].minY);

    std::size_t left = begin;
    std::size_t right = end;

    while (left < right) {
      while (left < right &&
             (isVertical ? agents_[left]->position_.x_
                         : agents_[left]->position_.y_) < splitValue) {
        ++left;
      }

      while (right > left &&
             (isVertical ? agents_[right - 1U]->position_.x_
                         : agents_[right - 1U]->position_.y_) >= splitValue) {
        --right;
      }

      if (left < right) {
        std::swap(agents_[left], agents_[right - 1U]);
        ++left;
        --right;
      }
    }

    std::size_t leftSize = left - begin;

    if (leftSize == 0U) {
      ++leftSize;
      ++left;
      ++right;
    }

    agentTree_[node].left = node + 1U;
    agentTree_[node].right = node + 1U + (2U * leftSize - 1U);

    buildAgentTreeRecursive(begin, left, agentTree_[node].left);
    buildAgentTreeRecursive(left, end, agentTree_[node].right);
  }
}

// support for obstacles
void KdTree::computeObstacleNeighbors(Agent *agent, float rangeSq) const
{
  queryObstacleTreeRecursive(agent, rangeSq, obstacleTree_);
}

void KdTree::deleteObstacleTree(ObstacleTreeNode *node)
{
  if (node != NULL) {
    deleteObstacleTree(node->left);
    deleteObstacleTree(node->right);
    delete node;
  }
}

void KdTree::queryObstacleTreeRecursive(Agent *agent, float rangeSq, const ObstacleTreeNode *node) const
{
  if (node == NULL) {
    return;
  }
  else {
    const Obstacle *const obstacle1 = node->obstacle;
    const Obstacle *const obstacle2 = obstacle1->nextObstacle_;

    const float agentLeftOfLine = leftOfObs(obstacle1->point_, obstacle2->point_, agent->position_);

    queryObstacleTreeRecursive(agent, rangeSq, (agentLeftOfLine >= 0.0f ? node->left : node->right));

    const float distSqLine = agentLeftOfLine * agentLeftOfLine / absSq(obstacle2->point_ - obstacle1->point_);

    if (distSqLine < rangeSq) {
      if (agentLeftOfLine < 0.0f) {
        /*
          * Try obstacle at this node only if agent is on right side of
          * obstacle (and can see obstacle).
          */
        agent->insertObstacleNeighbor(node->obstacle, rangeSq);
      }

      /* Try other side of line. */
      queryObstacleTreeRecursive(agent, rangeSq, (agentLeftOfLine >= 0.0f ? node->right : node->left));

    }
  }
}

void KdTree::buildObstacleTree()
{
  deleteObstacleTree(obstacleTree_);

  std::vector<Obstacle *> obstacles(simulator_->obstacles_.size());

  for (size_t i = 0; i < simulator_->obstacles_.size(); ++i) {
    obstacles[i] = simulator_->obstacles_[i];
  }

  obstacleTree_ = buildObstacleTreeRecursive(obstacles);
}

KdTree::ObstacleTreeNode *KdTree::buildObstacleTreeRecursive(const std::vector<Obstacle *> &obstacles)
{
  if (obstacles.empty()) {
    return NULL;
  }
  else {
    ObstacleTreeNode *const node = new ObstacleTreeNode;

    size_t optimalSplit = 0;
    size_t minLeft = obstacles.size();
    size_t minRight = obstacles.size();

    for (size_t i = 0; i < obstacles.size(); ++i) {
      size_t leftSize = 0;
      size_t rightSize = 0;

      const Obstacle *const obstacleI1 = obstacles[i];
      const Obstacle *const obstacleI2 = obstacleI1->nextObstacle_;

      /* Compute optimal split node. */
      for (size_t j = 0; j < obstacles.size(); ++j) {
        if (i == j) {
          continue;
        }

        const Obstacle *const obstacleJ1 = obstacles[j];
        const Obstacle *const obstacleJ2 = obstacleJ1->nextObstacle_;

        const float j1LeftOfI = leftOfObs(obstacleI1->point_, obstacleI2->point_, obstacleJ1->point_);
        const float j2LeftOfI = leftOfObs(obstacleI1->point_, obstacleI2->point_, obstacleJ2->point_);

        if (j1LeftOfI >= -AVO_EPSILON && j2LeftOfI >= -AVO_EPSILON) {
          ++leftSize;
        }
        else if (j1LeftOfI <= AVO_EPSILON && j2LeftOfI <= AVO_EPSILON) {
          ++rightSize;
        }
        else {
          ++leftSize;
          ++rightSize;
        }

        if (std::make_pair(std::max(leftSize, rightSize), std::min(leftSize, rightSize)) >= std::make_pair(std::max(minLeft, minRight), std::min(minLeft, minRight))) {
          break;
        }
      }

      if (std::make_pair(std::max(leftSize, rightSize), std::min(leftSize, rightSize)) < std::make_pair(std::max(minLeft, minRight), std::min(minLeft, minRight))) {
        minLeft = leftSize;
        minRight = rightSize;
        optimalSplit = i;
      }
    }

    /* Build split node. */
    std::vector<Obstacle *> leftObstacles(minLeft);
    std::vector<Obstacle *> rightObstacles(minRight);

    size_t leftCounter = 0;
    size_t rightCounter = 0;
    const size_t i = optimalSplit;

    const Obstacle *const obstacleI1 = obstacles[i];
    const Obstacle *const obstacleI2 = obstacleI1->nextObstacle_;

    for (size_t j = 0; j < obstacles.size(); ++j) {
      if (i == j) {
        continue;
      }

      Obstacle *const obstacleJ1 = obstacles[j];
      Obstacle *const obstacleJ2 = obstacleJ1->nextObstacle_;

      const float j1LeftOfI = leftOfObs(obstacleI1->point_, obstacleI2->point_, obstacleJ1->point_);
      const float j2LeftOfI = leftOfObs(obstacleI1->point_, obstacleI2->point_, obstacleJ2->point_);

      if (j1LeftOfI >= -AVO_EPSILON && j2LeftOfI >= -AVO_EPSILON) {
        leftObstacles[leftCounter++] = obstacles[j];
      }
      else if (j1LeftOfI <= AVO_EPSILON && j2LeftOfI <= AVO_EPSILON) {
        rightObstacles[rightCounter++] = obstacles[j];
      }
      else {
        /* Split obstacle j. */
        const float t = det(obstacleI2->point_ - obstacleI1->point_, obstacleJ1->point_ - obstacleI1->point_) / det(obstacleI2->point_ - obstacleI1->point_, obstacleJ1->point_ - obstacleJ2->point_);

        const Vector2 splitpoint = obstacleJ1->point_ + t * (obstacleJ2->point_ - obstacleJ1->point_);

        Obstacle *const newObstacle = new Obstacle();
        newObstacle->point_ = splitpoint;
        newObstacle->prevObstacle_ = obstacleJ1;
        newObstacle->nextObstacle_ = obstacleJ2;
        newObstacle->isConvex_ = true;
        newObstacle->unitDir_ = obstacleJ1->unitDir_;

        newObstacle->id_ = simulator_->obstacles_.size();

        simulator_->obstacles_.push_back(newObstacle);

        obstacleJ1->nextObstacle_ = newObstacle;
        obstacleJ2->prevObstacle_ = newObstacle;

        if (j1LeftOfI > 0.0f) {
          leftObstacles[leftCounter++] = obstacleJ1;
          rightObstacles[rightCounter++] = newObstacle;
        }
        else {
          rightObstacles[rightCounter++] = obstacleJ1;
          leftObstacles[leftCounter++] = newObstacle;
        }
      }
    }

    node->obstacle = obstacleI1;
    node->left = buildObstacleTreeRecursive(leftObstacles);
    node->right = buildObstacleTreeRecursive(rightObstacles);
    return node;
  }
}

void KdTree::computeAgentNeighbors(Agent *agent, float &rangeSq) const {
  queryAgentTreeRecursive(agent, rangeSq, 0U);
  }

void KdTree::queryAgentTreeRecursive(Agent *agent, float &rangeSq,
                                     std::size_t node) const {
  if (agentTree_[node].end - agentTree_[node].begin <= AVO_MAX_LEAF_SIZE) {
    for (std::size_t i = agentTree_[node].begin; i < agentTree_[node].end;
         ++i) {
      agent->insertAgentNeighbor(agents_[i], rangeSq);
    }
  } else {
    const float distLeftMinX = std::max(
        0.0F, agentTree_[agentTree_[node].left].minX - agent->position_.x_);
    const float distLeftMaxX = std::max(
        0.0F, agent->position_.x_ - agentTree_[agentTree_[node].left].maxX);
    const float distLeftMinY = std::max(
        0.0F, agentTree_[agentTree_[node].left].minY - agent->position_.y_);
    const float distLeftMaxY = std::max(
        0.0F, agent->position_.y_ - agentTree_[agentTree_[node].left].maxY);

    const float distSqLeft =
        distLeftMinX * distLeftMinX + distLeftMaxX * distLeftMaxX +
        distLeftMinY * distLeftMinY + distLeftMaxY * distLeftMaxY;

    const float distRightMinX = std::max(
        0.0F, agentTree_[agentTree_[node].right].minX - agent->position_.x_);
    const float distRightMaxX = std::max(
        0.0F, agent->position_.x_ - agentTree_[agentTree_[node].right].maxX);
    const float distRightMinY = std::max(
        0.0F, agentTree_[agentTree_[node].right].minY - agent->position_.y_);
    const float distRightMaxY = std::max(
        0.0F, agent->position_.y_ - agentTree_[agentTree_[node].right].maxY);

    const float distSqRight =
        distRightMinX * distRightMinX + distRightMaxX * distRightMaxX +
        distRightMinY * distRightMinY + distRightMaxY * distRightMaxY;

    if (distSqLeft < distSqRight) {
      if (distSqLeft < rangeSq) {
        queryAgentTreeRecursive(agent, rangeSq, agentTree_[node].left);

        if (distSqRight < rangeSq) {
          queryAgentTreeRecursive(agent, rangeSq, agentTree_[node].right);
        }
      }
    } else {
      if (distSqRight < rangeSq) {
        queryAgentTreeRecursive(agent, rangeSq, agentTree_[node].right);

        if (distSqLeft < rangeSq) {
          queryAgentTreeRecursive(agent, rangeSq, agentTree_[node].left);
        }
      }
    }
  }
}
}  // namespace AVO
