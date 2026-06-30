/*
 * Agent.cc
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
#include <iostream>
#include "Agent.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "KdTree.h"

namespace AVO {
namespace {
const std::size_t AVO_MAX_STEPS = 100U;

Vector2 tangentPoint(const Vector2 &p, float r) {
  const float pLengthSq = absSq(p);
  const float l = std::sqrt(pLengthSq - r * r);
  return Vector2(l * p.x_ - r * p.y_, r * p.x_ + l * p.y_) * (l / pLengthSq);
}

Vector2 boundaryAVO(const Vector2 &p, const Vector2 &v, float r, float delta,
                    float t) {
  const float e = std::exp(-t / delta) - 1.0F;
  const float temp = 1.0F / (t + delta * e);

  const Vector2 dc = (-e * p - (delta * e + (e + 1.0F) * t) * v) * temp * temp;
  const Vector2 rdrdc = (delta + t / e) * dc;

  const Vector2 c = (delta * e * v - p) * temp;
  const float radius = r * temp;

  return c - rdrdc + tangentPoint(rdrdc, radius) - v;
}

Vector2 centerAVO(const Vector2 &p, const Vector2 &v, float delta, float t) {
  const float temp = delta * (std::exp(-t / delta) - 1.0F);
  return (temp * v - p) / (t + temp) - v;
}

std::pair<Vector2, Vector2> circleCircleIntersection(float R, float r,
                                                     const Vector2 &p) {
  const float RSq = R * R;
  const float rSq = r * r;
  const float pxSq = p.x_ * p.x_;
  const float pySq = p.y_ * p.y_;
  const float squareRoot = std::sqrt(
      -pxSq * (pxSq * pxSq + pySq * pySq + 2.0F * pxSq * (pySq - rSq - RSq) +
               (rSq - RSq) * (rSq - RSq) - 2.0F * pySq * (rSq + RSq)));

  const float xNum = pxSq * pxSq + pxSq * pySq - pxSq * rSq + pxSq * RSq;
  const float yNum = pxSq * p.y_ + p.y_ * pySq - p.y_ * rSq + p.y_ * RSq;
  const float yDenom = 2.0F * (pxSq + pySq);
  const float xDenom = yDenom * p.x_;

  return std::make_pair(Vector2((xNum - p.y_ * squareRoot) / xDenom,
                                (yNum + squareRoot) / yDenom),
                        Vector2((xNum + p.y_ * squareRoot) / xDenom,
                                (yNum - squareRoot) / yDenom));
}

float distSqPointLineSegment(const Vector2 &vector1, const Vector2 &vector2,
                             const Vector2 &vector3) {
  const float r =
      ((vector3 - vector1) * (vector2 - vector1)) / absSq(vector2 - vector1);

  if (r < 0.0F) {
    return absSq(vector3 - vector1);
  }
  if (r > 1.0F) {
    return absSq(vector3 - vector2);
  }
  return absSq(vector3 - (vector1 + r * (vector2 - vector1)));
}

float leftOf(const Vector2 &vector1, const Vector2 &vector2,
             const Vector2 &vector3) {
  return det(vector1 - vector3, vector2 - vector1);
}

float radiusAVO(float r, float delta, float t) {
  return r / (t + delta * (std::exp(-t / delta) - 1.0F));
}

bool linearProgram1(const std::vector<Line> &lines, std::size_t lineNo,
                    float radius, const Vector2 &optVelocity, bool directionOpt,
                    Vector2 &result) {  // NOLINT(runtime/references)
  const float dotProduct = lines[lineNo].point * lines[lineNo].direction;
  const float discriminant =
      dotProduct * dotProduct + radius * radius - absSq(lines[lineNo].point);

  if (discriminant < 0.0F) {
    // Max speed circle fully invalidates line lineNo.
    return false;
  }

  const float sqrtDiscriminant = std::sqrt(discriminant);
  float tLeft = -dotProduct - sqrtDiscriminant;
  float tRight = -dotProduct + sqrtDiscriminant;

  for (std::size_t i = 0U; i < lineNo; ++i) {
    const float denominator = det(lines[lineNo].direction, lines[i].direction);
    const float numerator =
        det(lines[i].direction, lines[lineNo].point - lines[i].point);

    if (std::fabs(denominator) <= AVO_EPSILON) {
      // Lines lineNo and i are (almost) parallel.
      if (numerator < 0.0F) {
        return false;
      }
      continue;
    }

    const float t = numerator / denominator;

    if (denominator >= 0.0F) {
      // Line i bounds line lineNo on the right.
      tRight = std::min(tRight, t);
    } else {
      // Line i bounds line lineNo on the left.
      tLeft = std::max(tLeft, t);
    }

    if (tLeft > tRight) {
      return false;
    }
  }

  if (directionOpt) {
    // Optimize direction.
    if (optVelocity * lines[lineNo].direction > 0.0F) {
      // Take right extreme.
      result = lines[lineNo].point + tRight * lines[lineNo].direction;
    } else {
      // Take left extreme.
      result = lines[lineNo].point + tLeft * lines[lineNo].direction;
    }
  } else {
    // Optimize closest point.
    const float t =
        lines[lineNo].direction * (optVelocity - lines[lineNo].point);

    if (t < tLeft) {
      result = lines[lineNo].point + tLeft * lines[lineNo].direction;
    } else if (t > tRight) {
      result = lines[lineNo].point + tRight * lines[lineNo].direction;
    } else {
      result = lines[lineNo].point + t * lines[lineNo].direction;
    }
  }

  return true;
}

std::size_t linearProgram2(const std::vector<Line> &lines, float radius,
                           const Vector2 &optVelocity, bool directionOpt,
                           Vector2 &result) {  // NOLINT(runtime/references)
  if (directionOpt) {
    // Optimize direction. The optimization velocity is of unit length in this
    // case.
    result = optVelocity * radius;
  } else if (absSq(optVelocity) > radius * radius) {
    // Optimize closest point and outside circle.
    result = normalize(optVelocity) * radius;
  } else {
    // Optimize closest point and inside circle.
    result = optVelocity;
  }

  for (std::size_t lineNo = 0U; lineNo < lines.size(); ++lineNo) {
    if (det(lines[lineNo].direction, lines[lineNo].point - result) > 0.0F) {
      // Result does not satisfy constraint i. Compute new optimal result.
      const Vector2 tempResult = result;

      if (!linearProgram1(lines, lineNo, radius, optVelocity, directionOpt,
                          result)) {
        result = tempResult;
        return lineNo;
      }
    }
  }

  return lines.size();
}

void linearProgram3(const std::vector<Line> &lines, std::size_t numObstLines,
                    float radius,
                    Vector2 &result) {  // NOLINT(runtime/references)
  float distance = 0.0F;

  for (std::size_t lineNo = numObstLines; lineNo < lines.size(); ++lineNo) {
    if (det(lines[lineNo].direction, lines[lineNo].point - result) > distance) {
      // Result does not satisfy constraint of line.
      std::vector<Line> projLines(
          lines.begin(),
          lines.begin() + static_cast<std::ptrdiff_t>(numObstLines));

      for (std::size_t j = numObstLines; j < lineNo; ++j) {
        Line line;

        float determinant = det(lines[lineNo].direction, lines[j].direction);

        if (std::fabs(determinant) <= AVO_EPSILON) {
          // Line i and line j are parallel.
          if (lines[lineNo].direction * lines[j].direction > 0.0F) {
            // Line i and line j point in the same direction.
            continue;
          }
          // Line i and line j point in opposite direction.
          line.point = 0.5F * (lines[lineNo].point + lines[j].point);

        } else {
          line.point =
              lines[lineNo].point +
              (det(lines[j].direction, lines[lineNo].point - lines[j].point) /
               determinant) *
                  lines[lineNo].direction;
        }

        line.direction =
            normalize(lines[j].direction - lines[lineNo].direction);
        projLines.push_back(line);
      }

      const Vector2 tempResult = result;

      if (linearProgram2(
              projLines, radius,
              Vector2(-lines[lineNo].direction.y_, lines[lineNo].direction.x_),
              true, result) < projLines.size()) {
        // This should in principle not happen. The result is by definition
        // already in the feasible region of this linear program. If it fails,
        // it is due to small floating point error, and the current result is
        // kept.
        result = tempResult;
      }

      distance = det(lines[lineNo].direction, lines[lineNo].point - result);
    }
  }
}
}  // namespace

Agent::Agent()
    : id_(0U),
      maxNeighbors_(0U),
      accelInterval_(0.0F),
      maxAccel_(0.0F),
      maxSpeed_(0.0F),
      neighborDist_(0.0F),
      radius_(0.0F),
      timeHorizon_(0.0F),
      timeHorizonObst_(0.0F) {}

Agent::~Agent() {}

void Agent::computeNeighbors(const KdTree *kdTree) {
  obstacleNeighbors_.clear();
  float rangeSq = compute_sqr(timeHorizonObst_ * maxSpeed_ + radius_);
  kdTree->computeObstacleNeighbors(this, rangeSq);

  agentNeighbors_.clear();

  if (maxNeighbors_ > 0U) {
    float rangeSq = neighborDist_ * neighborDist_;
    kdTree->computeAgentNeighbors(this, rangeSq);
  }
}

void Agent::computeNewVelocity(float timeStep) {
  boundary_.clear();
  orcaLines_.clear();

	const float invTimeHorizonObst = 1.0f / timeHorizonObst_;
  float speed = abs(velocity_);
  Vector2 direction = velocity_ / speed;

  Line speedLine;

  // Insert maximum speed line.
  speedLine.point = (maxSpeed_ - speed) * direction;
  speedLine.direction = Vector2(-direction.y_, direction.x_);
  orcaLines_.push_back(speedLine);

  // Insert minimum speed line.
  speedLine.point = (maxSpeed_ - speed) * -direction;
  speedLine.direction = Vector2(direction.y_, -direction.x_);
  orcaLines_.push_back(speedLine);

  // const std::size_t numObstLines = orcaLines_.size();

  Line line;

  // Create obstacle ORCA lines
  for (size_t i = 0; i < obstacleNeighbors_.size(); ++i) {

    const Obstacle *obstacle1 = obstacleNeighbors_[i].second;
    const Obstacle *obstacle2 = obstacle1->nextObstacle_;

    const Vector2 relativePosition1 = obstacle1->point_ - position_;
    const Vector2 relativePosition2 = obstacle2->point_ - position_;

    /*
    * Check if velocity obstacle of obstacle is already taken care of by
    * previously constructed obstacle ORCA lines.
    */
    for (size_t j = 0; j < orcaLines_.size(); ++j) {
      if (det(invTimeHorizonObst * relativePosition1 - orcaLines_[j].point, orcaLines_[j].direction) - invTimeHorizonObst * radius_ >= -AVO_EPSILON && det(invTimeHorizonObst * relativePosition2 - orcaLines_[j].point, orcaLines_[j].direction) - invTimeHorizonObst * radius_ >=  -AVO_EPSILON) {
        break;
      }
    }

    /* Not yet covered. Check for collisions. */

    const float distSq1 = absSq(relativePosition1);
    const float distSq2 = absSq(relativePosition2);

    const float radiusSq = compute_sqr(radius_);

    const Vector2 obstacleVector = obstacle2->point_ - obstacle1->point_;
    const float s = (-relativePosition1 * obstacleVector) / absSq(obstacleVector);
    const float distSqLine = absSq(-relativePosition1 - s * obstacleVector);
    Line line;

    if (s < 0.0f && distSq1 <= radiusSq) {
      /* Collision with left vertex. Ignore if non-convex. */
      if (obstacle1->isConvex_) {
        line.point = Vector2(0.0f, 0.0f);
        line.direction = normalize(Vector2(-relativePosition1.y_, relativePosition1.x_));
        orcaLines_.push_back(line);
      }

      continue;
    }
    else if (s > 1.0f && distSq2 <= radiusSq) {
      /* Collision with right vertex. Ignore if non-convex
        * or if it will be taken care of by neighoring obstace */
      if (obstacle2->isConvex_ && det(relativePosition2, obstacle2->unitDir_) >= 0.0f) {
        line.point = Vector2(0.0f, 0.0f);
        line.direction = normalize(Vector2(-relativePosition2.y_, relativePosition2.x_));
        orcaLines_.push_back(line);
      }

      continue;
    }
    else if (s >= 0.0f && s < 1.0f && distSqLine <= radiusSq) {
      /* Collision with obstacle segment. */
      line.point = Vector2(0.0f, 0.0f);
      line.direction = -obstacle1->unitDir_;
      orcaLines_.push_back(line);
      continue;
    }

    /*
      * No collision.
      * Compute legs. When obliquely viewed, both legs can come from a single
      * vertex. Legs extend cut-off line when nonconvex vertex.
      */

    Vector2 leftLegDirection, rightLegDirection;

    if (s < 0.0f && distSqLine <= radiusSq) {
      /*
        * Obstacle viewed obliquely so that left vertex
        * defines velocity obstacle.
        */
      if (!obstacle1->isConvex_) {
        /* Ignore obstacle. */
        continue;
      }

      obstacle2 = obstacle1;

      const float leg1 = std::sqrt(distSq1 - radiusSq);
      leftLegDirection = Vector2(relativePosition1.x_ * leg1 - relativePosition1.y_ * radius_, relativePosition1.x_ * radius_ + relativePosition1.y_ * leg1) / distSq1;
      rightLegDirection = Vector2(relativePosition1.x_ * leg1 + relativePosition1.y_ * radius_, -relativePosition1.x_ * radius_ + relativePosition1.y_ * leg1) / distSq1;
    }
    else if (s > 1.0f && distSqLine <= radiusSq) {
      /*
        * Obstacle viewed obliquely so that
        * right vertex defines velocity obstacle.
        */
      if (!obstacle2->isConvex_) {
        /* Ignore obstacle. */
        continue;
      }

      obstacle1 = obstacle2;

      const float leg2 = std::sqrt(distSq2 - radiusSq);
      leftLegDirection = Vector2(relativePosition2.x_ * leg2 - relativePosition2.y_ * radius_, relativePosition2.x_ * radius_ + relativePosition2.y_ * leg2) / distSq2;
      rightLegDirection = Vector2(relativePosition2.x_ * leg2 + relativePosition2.y_ * radius_, -relativePosition2.x_ * radius_ + relativePosition2.y_ * leg2) / distSq2;
    }
    else {
      /* Usual situation. */
      if (obstacle1->isConvex_) {
        const float leg1 = std::sqrt(distSq1 - radiusSq);
        leftLegDirection = Vector2(relativePosition1.x_ * leg1 - relativePosition1.y_ * radius_, relativePosition1.x_ * radius_ + relativePosition1.y_ * leg1) / distSq1;
      }
      else {
        /* Left vertex non-convex; left leg extends cut-off line. */
        leftLegDirection = -obstacle1->unitDir_;
      }

      if (obstacle2->isConvex_) {
        const float leg2 = std::sqrt(distSq2 - radiusSq);
        rightLegDirection = Vector2(relativePosition2.x_ * leg2 + relativePosition2.y_ * radius_, -relativePosition2.x_ * radius_ + relativePosition2.y_ * leg2) / distSq2;
      }
      else {
        /* Right vertex non-convex; right leg extends cut-off line. */
        rightLegDirection = obstacle1->unitDir_;
      }
    }

    /*
      * Legs can never point into neighboring edge when convex vertex,
      * take cutoff-line of neighboring edge instead. If velocity projected on
      * "foreign" leg, no constraint is added.
      */

    const Obstacle *const leftNeighbor = obstacle1->prevObstacle_;

    bool isLeftLegForeign = false;
    bool isRightLegForeign = false;

    if (obstacle1->isConvex_ && det(leftLegDirection, -leftNeighbor->unitDir_) >= 0.0f) {
      /* Left leg points into obstacle. */
      leftLegDirection = -leftNeighbor->unitDir_;
      isLeftLegForeign = true;
    }

    if (obstacle2->isConvex_ && det(rightLegDirection, obstacle2->unitDir_) <= 0.0f) {
      /* Right leg points into obstacle. */
      rightLegDirection = obstacle2->unitDir_;
      isRightLegForeign = true;
    }

    /* Compute cut-off centers. */
    const Vector2 leftCutoff = invTimeHorizonObst * (obstacle1->point_ - position_);
    const Vector2 rightCutoff = invTimeHorizonObst * (obstacle2->point_ - position_);
    const Vector2 cutoffVec = rightCutoff - leftCutoff;

    /* Project current velocity on velocity obstacle. */

    /* Check if current velocity is projected on cutoff circles. */
    const float t = (obstacle1 == obstacle2 ? 0.5f : ((velocity_ - leftCutoff) * cutoffVec) / absSq(cutoffVec));
    const float tLeft = ((velocity_ - leftCutoff) * leftLegDirection);
    const float tRight = ((velocity_ - rightCutoff) * rightLegDirection);

    if ((t < 0.0f && tLeft < 0.0f) || (obstacle1 == obstacle2 && tLeft < 0.0f && tRight < 0.0f)) {
      /* Project on left cut-off circle. */
      const Vector2 unitW = normalize(velocity_ - leftCutoff);

      line.direction = Vector2(unitW.y_, -unitW.x_);
      line.point = leftCutoff + radius_ * invTimeHorizonObst * unitW;
      orcaLines_.push_back(line);
      continue;
    }
    else if (t > 1.0f && tRight < 0.0f) {
      /* Project on right cut-off circle. */
      const Vector2 unitW = normalize(velocity_ - rightCutoff);

      line.direction = Vector2(unitW.y_, -unitW.x_);
      line.point = rightCutoff + radius_ * invTimeHorizonObst * unitW;
      orcaLines_.push_back(line);
      continue;
    }

    /*
      * Project on left leg, right leg, or cut-off line, whichever is closest
      * to velocity.
      */
    const float distSqCutoff = ((t < 0.0f || t > 1.0f || obstacle1 == obstacle2) ? std::numeric_limits<float>::infinity() : absSq(velocity_ - (leftCutoff + t * cutoffVec)));
    const float distSqLeft = ((tLeft < 0.0f) ? std::numeric_limits<float>::infinity() : absSq(velocity_ - (leftCutoff + tLeft * leftLegDirection)));
    const float distSqRight = ((tRight < 0.0f) ? std::numeric_limits<float>::infinity() : absSq(velocity_ - (rightCutoff + tRight * rightLegDirection)));

    if (distSqCutoff <= distSqLeft && distSqCutoff <= distSqRight) {
      /* Project on cut-off line. */
      line.direction = -obstacle1->unitDir_;
      line.point = leftCutoff + radius_ * invTimeHorizonObst * Vector2(-line.direction.y_, line.direction.x_);
      orcaLines_.push_back(line);
      continue;
    }
    else if (distSqLeft <= distSqRight) {
      /* Project on left leg. */
      if (isLeftLegForeign) {
        continue;
      }

      line.direction = leftLegDirection;
      line.point = leftCutoff + radius_ * invTimeHorizonObst * Vector2(-line.direction.y_, line.direction.x_);
      orcaLines_.push_back(line);
      continue;
    }
    else {
      /* Project on right leg. */
      if (isRightLegForeign) {
        continue;
      }

      line.direction = -rightLegDirection;
      line.point = rightCutoff + radius_ * invTimeHorizonObst * Vector2(-line.direction.y_, line.direction.x_);
      orcaLines_.push_back(line);
      continue;
    }
  }
  
  std::cout << "obstacleNeighbors_.size() = "
          << obstacleNeighbors_.size()
          << std::endl;

  const std::size_t numObstLines = orcaLines_.size();

  std::cout << "numObstLines = "
          << numObstLines
          << std::endl;
          
  for (size_t i = 0; i < numObstLines; i++) {
    std::cout << "obst line " << i
              << " point=" << orcaLines_[i].point
              << " dir=" << orcaLines_[i].direction
              << std::endl;
  }
  // Create agent ORCA lines.
  for (std::size_t agentNo = 0U; agentNo < agentNeighbors_.size(); ++agentNo) {
    const Agent *const other = agentNeighbors_[agentNo].second;

    const Vector2 relativePosition = position_ - other->position_;
    const Vector2 relativeVelocity = velocity_ - other->velocity_;
    const float combinedRadius = radius_ + other->radius_;
    const float combinedMaxAccel = maxAccel_ + other->maxAccel_;
    const float relativePositionSq = absSq(relativePosition);
    const float combinedRadiusSq = compute_sqr(combinedRadius);
    std::deque<Vector2> boundary;
    bool convex = false;
    float left = 0.0F;  // 1 if left, -1 if right.
    Vector2 cutoffCenter;
    float cutoffRadius = 0.0F;

    // Check if left boundary or right boundary matters.
    const float determinant = det(relativePosition, relativeVelocity);

    if (determinant > 0.0F) {
      // Right boundary matters.
      left = -1.0F;
    } else {
      // Left boundary matters.
      left = 1.0F;
    }

    if (relativePositionSq <= combinedRadiusSq) {
      // Collision.
      const Vector2 centerDt = centerAVO(relativePosition, relativeVelocity,
                                         accelInterval_, timeStep);
      const Vector2 centerTau = centerAVO(relativePosition, relativeVelocity,
                                          accelInterval_, timeHorizon_);
      const float radiusDt =
          radiusAVO(combinedRadius, accelInterval_, timeStep);
      const float radiusTau =
          radiusAVO(combinedRadius, accelInterval_, timeHorizon_);

      const Vector2 centerVec = centerDt - centerTau;
      const float centerVecSq = absSq(centerVec);
      const float diffRadius = left * (radiusDt - radiusTau);
      const float diffRadiusSq = diffRadius * diffRadius;

      if (centerVecSq > diffRadiusSq) {
        // Compute tangent line of circle at time timeStep and circle at time
        // timeHorizon as boundary.
        const float l = std::sqrt(centerVecSq - diffRadiusSq);
        const Vector2 unitDir =
            Vector2(l * centerVec.x_ - diffRadius * centerVec.y_,
                    diffRadius * centerVec.x_ + l * centerVec.y_) /
            centerVecSq;
        boundary.push_back(centerTau + radiusTau * Vector2(-left * unitDir.y_,
                                                           left * unitDir.x_));
        boundary.push_back(centerDt + radiusDt * Vector2(-left * unitDir.y_,
                                                         left * unitDir.x_));

        cutoffCenter = centerTau;
        cutoffRadius = radiusTau;
      } else {
        // Project on circle of time timeStep.
        cutoffCenter = centerDt;
        cutoffRadius = radiusDt;
      }
    } else {
      // No collision. Check if current velocity in VO to see if convex.
      const float discr = combinedRadiusSq * absSq(relativeVelocity) -
                          determinant * determinant;

      if (discr >= 0.0F) {
        convex = std::sqrt(discr) - relativePosition * relativeVelocity >= 0.0F;
      } else {
        convex = false;
      }

      float lastTime = timeStep;

      for (std::size_t i = 0U; i <= AVO_MAX_STEPS; ++i) {
        const float t = timeStep + static_cast<float>(i) *
                                       (timeHorizon_ - timeStep) /
                                       AVO_MAX_STEPS;
        boundary.push_front(boundaryAVO(relativePosition, relativeVelocity,
                                        left * combinedRadius, accelInterval_,
                                        t));
        float boundaryFrontX = boundary.front().x_;

        if (boundaryFrontX != boundaryFrontX) {
          // Boundary not well-defined.
          boundary.pop_front();
          convex = true;
          break;
        }
        if (convex && i >= 2U &&
            left * leftOf(boundary[0], boundary[1], boundary[2]) >
                AVO_EPSILON) {
          // Failed to detect undefined part of boundary, and now running into
          // nonconvex vertex.
          boundary.pop_front();
          boundary.pop_front();
          break;
        }

        lastTime = t;
      }

      cutoffCenter = centerAVO(relativePosition, relativeVelocity,
                               accelInterval_, lastTime);
      cutoffRadius = radiusAVO(combinedRadius, accelInterval_, lastTime);

      if (!convex) {
        // Find intersections of non-convex boundary.
        std::vector<Vector2> intersections;
        const float rSq = accelInterval_ * accelInterval_ * combinedMaxAccel *
                          combinedMaxAccel;

        for (int i = 0; i < static_cast<int>(boundary.size()) - 1; ++i) {
          const Vector2 p = boundary[i];
          const Vector2 v = boundary[i + 1] - boundary[i];
          const float vSq = absSq(v);
          const float detPv = det(p, v);
          const float discr2 = rSq * vSq - detPv * detPv;

          const float distSq1 = absSq(boundary[i]);
          const float distSq2 = absSq(boundary[i + 1]);

          if (discr2 > 0.0F && ((distSq1 < rSq && distSq2 > rSq) ||
                                (distSq1 >= rSq && distSq2 <= rSq) ||
                                (distSq1 >= rSq && distSq2 > rSq &&
                                 absSq(p - (p * v / vSq) * v) < rSq))) {
            const float discrSqrt = std::sqrt(discr2);
            const float t1 = -(p * v + discrSqrt) / vSq;
            const float t2 = -(p * v - discrSqrt) / vSq;

            if (t1 >= 0.0F && t1 < 1.0F) {
              // Segment intersects disc.
              intersections.push_back(p + t1 * v);
            }

            if (t2 >= 0.0F && t2 < 1.0F) {
              // Segment intersects disc.
              intersections.push_back(p + t2 * v);
            }
          }
        }

        const float velocityPlusCutoffRadius =
            accelInterval_ * combinedMaxAccel + cutoffRadius;
        const bool cutoffInDisc =
            absSq(cutoffCenter) <
                velocityPlusCutoffRadius * velocityPlusCutoffRadius &&
            (absSq(boundary[0]) < rSq ||
             left * leftOf(Vector2(), cutoffCenter, boundary[0]) > 0.0F);

        if (!cutoffInDisc) {
          if (!intersections.empty()) {
            // Line between first and last left intersection.
            line.direction =
                left * normalize(intersections.back() - intersections.front());
            line.point = 0.5F * intersections.front();
            orcaLines_.push_back(line);
            continue;
          }
          // Nothing in disc.
          continue;
        }
        // Cutoff arc intersects disc. Connect last intersection to disc and
        // replace boundary.
        if (intersections.empty() && absSq(boundary[0]) <= rSq) {
          // Boundary too short to intersect disc.
          intersections.push_back(boundary.back());
        } else if (!intersections.empty() &&
                   absSq(boundary[0] - intersections.back()) < AVO_EPSILON) {
          // Boundary starts exactly at disc.
          intersections.clear();
        }

        if (!intersections.empty()) {
          Vector2 q = tangentPoint(cutoffCenter - intersections.back(),
                                   -left * cutoffRadius) +
                      intersections.back();

          if (absSq(q) > rSq) {
            // Compute intersection point of cutoff arc and disc and create
            // line.
            std::pair<Vector2, Vector2> p = circleCircleIntersection(
                accelInterval_ * combinedMaxAccel, cutoffRadius, cutoffCenter);

            if (-left * leftOf(intersections.back(), cutoffCenter, p.first) >
                0.0F) {
              q = p.first;
            } else {
              q = p.second;
            }

            line.direction =
                static_cast<float>(left) * normalize(intersections.back() - q);
            line.point = 0.5F * q;
            orcaLines_.push_back(line);
            continue;
          }
          boundary[0] = q;
          boundary[1] = intersections.back();
          boundary.erase(boundary.begin() + 2U, boundary.end());
        } else {
          boundary.clear();
        }
      }
    }

    // Treat as convex obstacle.
    if (boundary.size() <= 1U ||
        left * leftOf(cutoffCenter, boundary[0], Vector2()) >= 0.0F) {
      // Project on cutoff center.
      const float wLength = abs(cutoffCenter);
      const Vector2 unitW = cutoffCenter / -wLength;

      line.direction = Vector2(unitW.y_, -unitW.x_);
      line.point = 0.5F * (cutoffRadius - wLength) * unitW;
      orcaLines_.push_back(line);
      continue;
    }
    // Find closest point on boundary.
    float minDistSq = std::numeric_limits<float>::infinity();
    std::size_t minSegment = 0U;

    for (int i = 0; i < static_cast<int>(boundary.size()) - 1; ++i) {
      const float distSq =
          distSqPointLineSegment(boundary[i], boundary[i + 1], Vector2());

      if (distSq < minDistSq) {
        minDistSq = distSq;
        minSegment = i;
      }
    }

    line.direction =
        static_cast<float>(left) *
        normalize(boundary[minSegment + 1U] - boundary[minSegment]);
    line.point = 0.5F * boundary[minSegment];
    orcaLines_.push_back(line);
  }

  if (linearProgram2(orcaLines_, maxAccel_ * accelInterval_,
                     prefVelocity_ - velocity_, false, newVelocity_) == 0U) {
    linearProgram3(orcaLines_, numObstLines, maxAccel_ * accelInterval_,
                   newVelocity_);
  }

  newVelocity_ += velocity_;
}

void Agent::insertAgentNeighbor(const Agent *agent, float &rangeSq) {
  if (this != agent) {
    const float distSq = absSq(position_ - agent->position_);

    if (distSq < rangeSq) {
      if (agentNeighbors_.size() < maxNeighbors_) {
        agentNeighbors_.push_back(std::make_pair(distSq, agent));
      }

      std::size_t agentNo = agentNeighbors_.size() - 1U;

      while (agentNo != 0 && distSq < agentNeighbors_[agentNo - 1U].first) {
        agentNeighbors_[agentNo] = agentNeighbors_[agentNo - 1U];
        --agentNo;
      }

      agentNeighbors_[agentNo] = std::make_pair(distSq, agent);

      if (agentNeighbors_.size() == maxNeighbors_) {
        rangeSq = agentNeighbors_.back().first;
      }
    }
  }
}

void Agent::update(float timeStep) {
  velocity_ = newVelocity_;
  position_ += velocity_ * timeStep;
}

void Agent::insertObstacleNeighbor(const Obstacle *obstacle, float rangeSq)
{
  const Obstacle *const nextObstacle = obstacle->nextObstacle_;

  const float distSq = distSqPointLineSegment(obstacle->point_, nextObstacle->point_, position_);
  if (distSq < rangeSq) {
    obstacleNeighbors_.push_back(std::make_pair(distSq, obstacle));

    size_t i = obstacleNeighbors_.size() - 1;

    while (i != 0 && distSq < obstacleNeighbors_[i - 1].first) {
      obstacleNeighbors_[i] = obstacleNeighbors_[i - 1];
      --i;
    }
    obstacleNeighbors_[i] = std::make_pair(distSq, obstacle);
  }
}


}  // namespace AVO
