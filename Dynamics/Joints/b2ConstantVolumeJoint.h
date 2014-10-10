#ifndef b2ConstantVolumeJoint_h__
#define b2ConstantVolumeJoint_h__

#include <vector>
#include <Box2D/Dynamics/Joints/b2Joint.h>

class b2World;
class b2DistanceJoint;

using namespace std;

struct b2ConstantVolumeJointDef : public b2JointDef {

public:

	b2ConstantVolumeJointDef();

	void AddBody(b2Body* body);
	
	float32 dampingRatio;
	float32 frequencyHq;

private: 

	friend class b2ConstantVolumeJoint;

	vector<b2Body*> bodies;
};

class b2ConstantVolumeJoint : public b2Joint {

public:

	virtual ~b2ConstantVolumeJoint();

	void Inflate(float32 factor);
	const vector<b2Body*>& GetBodies() const;

	virtual b2Vec2 GetAnchorA() const;
	virtual b2Vec2 GetAnchorB() const;

	virtual b2Vec2 GetReactionForce(float32 inv_dt) const;
	virtual float32 GetReactionTorque(float32 inv_dt) const;

protected:

	friend class b2Joint;

	b2ConstantVolumeJoint(b2World* world, const b2ConstantVolumeJointDef* def);

	virtual void InitVelocityConstraints(const b2SolverData& data);
	virtual void SolveVelocityConstraints(const b2SolverData& data);
	
	virtual bool SolvePositionConstraints(const b2SolverData& data);

private:  // methods

	void CreateDistanceJoints(const b2ConstantVolumeJointDef* def);


	float GetBodyArea();
	float GetSolverArea(const b2SolverData& data);

private: // fields

	b2World* world;

	float impulse;
	float targetVolume;	
	float rotationalAcceleration;

	vector<b2Vec2> normals;
	vector<float> targetLengths;

	vector<b2Vec2> deltaPositions;
	vector<b2Vec2> deltaVelocities;

	vector<b2Body*> bodies;
	vector<b2DistanceJoint*> distanceJoints;

};

#endif // b2ConstantVolumeJoint_h__