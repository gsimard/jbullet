package com.bulletphysics.dynamics.constraintsolver;

import com.bulletphysics.dynamics.RigidBody;
import com.bulletphysics.linearmath.Transform;
import com.bulletphysics.linearmath.VectorUtil;
import cz.advel.stack.Stack;
import javax.vecmath.Matrix3d;
import javax.vecmath.Vector3d;

/**
 * Rope-type constraint between two rigid bodies each with a pivot point.
 *
 * @author simard
 */
public class RopeConstraint extends TypedConstraint {

	private final JacobianEntry[] jac = new JacobianEntry[]/*[3]*/ { new JacobianEntry(), new JacobianEntry(), new JacobianEntry() }; // 3 orthogonal linear constraints

	private final Vector3d pivotInA = new Vector3d();
	private final Vector3d pivotInB = new Vector3d();

    private double distance = 0f;

	public ConstraintSetting setting = new ConstraintSetting();

	public RopeConstraint() {
		super(TypedConstraintType.ROPE_CONSTRAINT_TYPE);
	}

	public RopeConstraint(RigidBody rbA, RigidBody rbB, Vector3d pivotInA, Vector3d pivotInB, double distance) {
		super(TypedConstraintType.ROPE_CONSTRAINT_TYPE, rbA, rbB);
		this.pivotInA.set(pivotInA);
		this.pivotInB.set(pivotInB);
        this.distance = distance;
	}

	public RopeConstraint(RigidBody rbA, Vector3d pivotInA) {
		super(TypedConstraintType.ROPE_CONSTRAINT_TYPE, rbA);
		this.pivotInA.set(pivotInA);
		this.pivotInB.set(pivotInA);
		rbA.getCenterOfMassTransform(Stack.alloc(Transform.class)).transform(this.pivotInB);
	}

	@Override
	public void buildJacobian() {
		appliedImpulse = 0f;

		Vector3d normal = Stack.alloc(Vector3d.class);
		normal.set(0f, 0f, 0f);

		Matrix3d tmpMat1 = Stack.alloc(Matrix3d.class);
		Matrix3d tmpMat2 = Stack.alloc(Matrix3d.class);
		Vector3d tmp1 = Stack.alloc(Vector3d.class);
		Vector3d tmp2 = Stack.alloc(Vector3d.class);
		Vector3d tmpVec = Stack.alloc(Vector3d.class);

		Transform centerOfMassA = rbA.getCenterOfMassTransform(Stack.alloc(Transform.class));
		Transform centerOfMassB = rbB.getCenterOfMassTransform(Stack.alloc(Transform.class));

		for (int i = 0; i < 3; i++) {
			VectorUtil.setCoord(normal, i, 1f);

			tmpMat1.transpose(centerOfMassA.basis);
			tmpMat2.transpose(centerOfMassB.basis);

			tmp1.set(pivotInA);
			centerOfMassA.transform(tmp1);
			tmp1.sub(rbA.getCenterOfMassPosition(tmpVec));

			tmp2.set(pivotInB);
			centerOfMassB.transform(tmp2);
			tmp2.sub(rbB.getCenterOfMassPosition(tmpVec));

			jac[i].init(
					tmpMat1,
					tmpMat2,
					tmp1,
					tmp2,
					normal,
					rbA.getInvInertiaDiagLocal(Stack.alloc(Vector3d.class)),
					rbA.getInvMass(),
					rbB.getInvInertiaDiagLocal(Stack.alloc(Vector3d.class)),
					rbB.getInvMass());
			VectorUtil.setCoord(normal, i, 0f);
		}
	}

	@Override
	public void solveConstraint(double timeStep) {
		Vector3d tmp = Stack.alloc(Vector3d.class);
		Vector3d tmp2 = Stack.alloc(Vector3d.class);
		Vector3d tmpVec = Stack.alloc(Vector3d.class);

		Transform centerOfMassA = rbA.getCenterOfMassTransform(Stack.alloc(Transform.class));
		Transform centerOfMassB = rbB.getCenterOfMassTransform(Stack.alloc(Transform.class));

		Vector3d pivotAInW = Stack.alloc(pivotInA);
		centerOfMassA.transform(pivotAInW);

		Vector3d pivotBInW = Stack.alloc(pivotInB);
		centerOfMassB.transform(pivotBInW);

		Vector3d normal = Stack.alloc(Vector3d.class);
		normal.set(0f, 0f, 0f);

		//btVector3 angvelA = m_rbA.getCenterOfMassTransform().getBasis().transpose() * m_rbA.getAngularVelocity();
		//btVector3 angvelB = m_rbB.getCenterOfMassTransform().getBasis().transpose() * m_rbB.getAngularVelocity();

		for (int i = 0; i < 3; i++) {
			VectorUtil.setCoord(normal, i, 1f);
			double jacDiagABInv = 1f / jac[i].getDiagonal();

			Vector3d rel_pos1 = Stack.alloc(Vector3d.class);
			rel_pos1.sub(pivotAInW, rbA.getCenterOfMassPosition(tmpVec));
			Vector3d rel_pos2 = Stack.alloc(Vector3d.class);
			rel_pos2.sub(pivotBInW, rbB.getCenterOfMassPosition(tmpVec));
			// this jacobian entry could be re-used for all iterations

			Vector3d vel1 = rbA.getVelocityInLocalPoint(rel_pos1, Stack.alloc(Vector3d.class));
			Vector3d vel2 = rbB.getVelocityInLocalPoint(rel_pos2, Stack.alloc(Vector3d.class));
			Vector3d vel = Stack.alloc(Vector3d.class);
			vel.sub(vel1, vel2);

			double rel_vel;
			rel_vel = normal.dot(vel);

			/*
			//velocity error (first order error)
			btScalar rel_vel = m_jac[i].getRelativeVelocity(m_rbA.getLinearVelocity(),angvelA,
			m_rbB.getLinearVelocity(),angvelB);
			 */

			// positional error (zeroth order error)
			tmp.sub(pivotAInW, pivotBInW);

            // GS: For a rope, don't correct anything if the rope is not stretched
            if (tmp.length() < distance) {
                tmp.set(0f,0f,0f);
            } else {
                // If the rope is streched, correct, but begin by
                // substracting rope length to prevent huge
                // unwarranted errors
                Vector3d rope = Stack.alloc(Vector3d.class);
                rope.normalize(tmp);
                rope.scale(distance);

                tmp.sub(rope);
            }

			double depth = -tmp.dot(normal); //this is the error projected on the normal

			//double impulse = depth * setting.tau / timeStep * jacDiagABInv - setting.damping * rel_vel * jacDiagABInv;
			double impulse = depth * setting.tau / timeStep * jacDiagABInv;

            //System.out.println("impulse: " + impulse);

			double impulseClamp = setting.impulseClamp;
			if (impulseClamp > 0f) {
				if (impulse < -impulseClamp) {
					impulse = -impulseClamp;
				}
				if (impulse > impulseClamp) {
					impulse = impulseClamp;
				}
			}

			appliedImpulse += impulse;
			Vector3d impulse_vector = Stack.alloc(Vector3d.class);
			impulse_vector.scale(impulse, normal);
			tmp.sub(pivotAInW, rbA.getCenterOfMassPosition(tmpVec));
			rbA.applyImpulse(impulse_vector, tmp);
			tmp.negate(impulse_vector);
			tmp2.sub(pivotBInW, rbB.getCenterOfMassPosition(tmpVec));
			rbB.applyImpulse(tmp, tmp2);

			VectorUtil.setCoord(normal, i, 0f);
		}
	}

	public void updateRHS(double timeStep) {
	}

	public void setPivotA(Vector3d pivotA) {
		pivotInA.set(pivotA);
	}

	public void setPivotB(Vector3d pivotB) {
		pivotInB.set(pivotB);
	}

	public Vector3d getPivotInA(Vector3d out) {
		out.set(pivotInA);
		return out;
	}

	public Vector3d getPivotInB(Vector3d out) {
		out.set(pivotInB);
		return out;
	}

	////////////////////////////////////////////////////////////////////////////

	public static class ConstraintSetting {
		public double tau = 0.0001f; //GS
		public double damping = 1f;
		public double impulseClamp = 0f;
	}

}
