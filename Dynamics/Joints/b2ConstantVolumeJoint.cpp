#include <Box2D/Dynamics/Joints/b2ConstantVolumeJoint.h>
#include <Box2D/Dynamics/b2Body.h>
#include <Box2D/Dynamics/b2Fixture.h>
#include <Box2D/Collision/Shapes/b2CircleShape.h>
#include <Box2D/Dynamics/b2TimeStep.h>
#include <Box2D/Dynamics/Joints/b2DistanceJoint.h>
#include <Box2D/Dynamics/b2World.h>

#include <iostream>

b2ConstantVolumeJointDef::b2ConstantVolumeJointDef()
	: frequencyHq(0.0f), dampingRatio(0.0f) {

	type = e_constantVolumeJoint;
}

void b2ConstantVolumeJointDef::AddBody(b2Body* body) {

	bodies.push_back(body);

	if (bodies.size() == 1) {

		bodyA = body;
	}
	else if (bodies.size() == 2) {

		bodyB = body;
	}
}

b2ConstantVolumeJoint::b2ConstantVolumeJoint(b2World* world, const b2ConstantVolumeJointDef* def ): 
	b2Joint(def), 
	world(world),
	rotationalAcceleration(0),
	impulse(0),
	bodies(def->bodies),
	normals(def->bodies.size(), b2Vec2(0, 0)),
	targetLengths(def->bodies.size(), 0),
	deltaPositions(def->bodies.size(), b2Vec2(0, 0)),
	deltaVelocities(def->bodies.size(), b2Vec2(0, 0)) {

	b2Vec2 distance;

	for (unsigned i = 0; i < targetLengths.size(); ++i) {

		const unsigned next = (i == targetLengths.size() - 1) ? 0 : i + 1;

		distance = bodies[i]->GetWorldCenter() - bodies[next]->GetWorldCenter();
		float dist = distance.Length();

		targetLengths[i] = dist;
	}

	targetVolume = GetBodyArea();

	CreateDistanceJoints(def);
}


b2ConstantVolumeJoint::~b2ConstantVolumeJoint() {

	for (unsigned i = 0; i < distanceJoints.size(); ++i) {

		world->DestroyJoint(distanceJoints[i]);
	}

	distanceJoints.clear();
}



void b2ConstantVolumeJoint::CreateDistanceJoints(const b2ConstantVolumeJointDef* def) {

	std::vector<b2Body*> bodies = def->bodies;

	b2DistanceJointDef jd;
	jd.bodyA = bodies[0];
	jd.bodyB = bodies[1];
	jd.dampingRatio = def->dampingRatio;
	jd.frequencyHz = def->frequencyHq;
	jd.collideConnected = def->collideConnected;

	jd.Initialize(bodies[0], bodies.back(), bodies.front()->GetWorldCenter(), bodies.back()->GetWorldCenter());

	distanceJoints.push_back(static_cast<b2DistanceJoint*>(world->CreateJoint(&jd)));

	for (unsigned i = 0; i < bodies.size() - 1; i++) {

		jd.Initialize(bodies[i], bodies[i + 1], bodies[i]->GetWorldCenter(), bodies[i + 1]->GetWorldCenter());
		distanceJoints.push_back(static_cast<b2DistanceJoint*>(world->CreateJoint(&jd)));
	}
}


float b2ConstantVolumeJoint::GetBodyArea() {

	float area = 0;

	area += (bodies[bodies.size() - 1]->GetWorldCenter().x * bodies[0]->GetWorldCenter().y) - (bodies[0]->GetWorldCenter().x * bodies[bodies.size() - 1]->GetWorldCenter().y);

	for(unsigned i = 0; i < bodies.size() - 1; ++i) {

		area += (bodies[i]->GetWorldCenter().x * bodies[i + 1]->GetWorldCenter().y) - (bodies[i + 1]->GetWorldCenter().x * bodies[i]->GetWorldCenter().y);
	}

	return area * 0.5f;
}

float b2ConstantVolumeJoint::GetSolverArea( const b2SolverData& data ) {

	float area = 0;

	int index = bodies[bodies.size() - 1]->m_islandIndex;
	int nextIndex = bodies[0]->m_islandIndex;;

	area += (data.positions[index].c.x * data.positions[nextIndex].c.y) - (data.positions[nextIndex].c.x	* data.positions[index].c.y);

	for(unsigned i = 0; i < bodies.size() - 1; ++i) {

		index = bodies[i]->m_islandIndex;
		nextIndex = bodies[i + 1]->m_islandIndex;

		area += (data.positions[index].c.x * data.positions[nextIndex].c.y)	- (data.positions[nextIndex].c.x * data.positions[index].c.y);
	}

	return area * 0.5f;
}

void b2ConstantVolumeJoint::InitVelocityConstraints( const b2SolverData& data ) {

	for (unsigned int i = 0; i < bodies.size(); ++i) {

		const unsigned prev = (i == 0) ? bodies.size() - 1 : i - 1;
		const unsigned next = (i == bodies.size() - 1) ? 0 : i + 1;

		deltaPositions[i] = data.positions[bodies[next]->m_islandIndex].c;
		deltaPositions[i] -= data.positions[bodies[prev]->m_islandIndex].c;
	}

	if (data.step.warmStarting) {

		impulse *= data.step.dtRatio;

		for (unsigned i = 0; i < bodies.size(); ++i) {

			data.velocities[bodies[i]->m_islandIndex].v.x += bodies[i]->m_invMass * deltaPositions[i].y * 0.5f * impulse;
			data.velocities[bodies[i]->m_islandIndex].v.y += bodies[i]->m_invMass * -deltaPositions[i].x * 0.5f * impulse;
		}
	} 
	else {
		impulse = 0.0f;
	}
}

void b2ConstantVolumeJoint::SolveVelocityConstraints( const b2SolverData& data ) {

	float crossMassSum = 0.0f;
	float dotMassSum = 0.0f;

	for (unsigned i = 0; i < bodies.size(); ++i) {

		const unsigned prev = (i == 0) ? bodies.size() - 1 : i - 1;
		const unsigned next = (i == bodies.size() - 1) ? 0 : i + 1;

		deltaVelocities[i] = data.positions[bodies[next]->m_islandIndex].c;
		deltaVelocities[i] -= data.positions[bodies[prev]->m_islandIndex].c;

		dotMassSum += (deltaVelocities[i].LengthSquared()) / bodies[i]->GetMass();
		crossMassSum += b2Cross(data.velocities[bodies[i]->m_islandIndex].v, deltaVelocities[i]);
	}

	float lambda = -2.0f * crossMassSum / dotMassSum;

	impulse += lambda;

	for (unsigned i = 0; i < bodies.size(); ++i) {

		data.velocities[bodies[i]->m_islandIndex].v.x += bodies[i]->m_invMass * deltaVelocities[i].y * 0.5f * lambda;
		data.velocities[bodies[i]->m_islandIndex].v.y += bodies[i]->m_invMass * -deltaVelocities[i].x * 0.5f * lambda;
	}
}

bool b2ConstantVolumeJoint::SolvePositionConstraints( const b2SolverData& data ) {

	float perimeter = 0.0f;

	for (unsigned i = 0; i < bodies.size(); ++i) {

		const unsigned next = (i == bodies.size() - 1) ? 0 : i + 1;

		float dx = data.positions[bodies[next]->m_islandIndex].c.x - data.positions[bodies[i]->m_islandIndex].c.x;
		float dy = data.positions[bodies[next]->m_islandIndex].c.y - data.positions[bodies[i]->m_islandIndex].c.y;

		float dist = sqrt(dx * dx + dy * dy);

		if (dist < b2_epsilon) {
			dist = 1.0f;
		}

		normals[i].x = dy / dist;
		normals[i].y = -dx / dist;

		perimeter += dist;
	}

	b2Vec2 delta;

	float deltaArea = targetVolume - GetSolverArea(data);
	float toExtrude = 0.5f * deltaArea / perimeter;

	bool done = true;

	for (unsigned i = 0; i < bodies.size(); ++i) {

		const unsigned next = (i == bodies.size() - 1) ? 0 : i + 1;

		delta.Set(toExtrude * (normals[i].x + normals[next].x), toExtrude * (normals[i].y + normals[next].y));

		float norm = delta.Length();

		if (norm > b2_maxLinearCorrection) {
			delta *= (b2_maxLinearCorrection / norm);
		}

		if (norm > b2_linearSlop) {
			done = false;
		}

		data.positions[bodies[next]->m_islandIndex].c += delta;
	}

	return done;
}

void b2ConstantVolumeJoint::Inflate(float32 factor) {
	targetVolume *= factor;
}

b2Vec2 b2ConstantVolumeJoint::GetAnchorA() const {
	return bodies[0]->GetPosition();
}

b2Vec2 b2ConstantVolumeJoint::GetAnchorB() const {
	return bodies[1]->GetPosition();
}

b2Vec2 b2ConstantVolumeJoint::GetReactionForce( float32 inv_dt ) const {
	B2_NOT_USED(inv_dt);
	return b2Vec2(0, 0);
}

float32 b2ConstantVolumeJoint::GetReactionTorque( float32 inv_dt ) const {
	B2_NOT_USED(inv_dt);
	return 0;
}

const vector<b2Body*>& b2ConstantVolumeJoint::GetBodies() const {
	return bodies;
}
