/*
 * Java port of Bullet (c) 2008 Martin Dvorak <jezek2@advel.cz>
 *
 * This source file is part of GIMPACT Library.
 *
 * For the latest info, see http://gimpact.sourceforge.net/
 *
 * Copyright (c) 2007 Francisco Leon Najera. C.C. 80087371.
 * email: projectileman@yahoo.com
 *
 * This software is provided 'as-is', without any express or implied warranty.
 * In no event will the authors be held liable for any damages arising from
 * the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgment in the product documentation would be
 *    appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 */

package com.bulletphysics.extras.gimpact;

import com.bulletphysics.BulletGlobals;
import com.bulletphysics.linearmath.VectorUtil;
import cz.advel.stack.Stack;
import javax.vecmath.Vector3d;
import javax.vecmath.Vector4d;

/**
 *
 * @author jezek2
 */
class GeometryOperations {

	public static final double PLANEDIREPSILON = 0.0000001f;
	public static final double PARALELENORMALS = 0.000001f;

	public static final double CLAMP(double number, double minval, double maxval) {
		return (number < minval? minval : (number > maxval? maxval : number));
	}

	/**
	 * Calc a plane from a triangle edge an a normal.
	 */
	public static void edge_plane(Vector3d e1, Vector3d e2, Vector3d normal, Vector4d plane) {
		Vector3d planenormal = Stack.alloc(Vector3d.class);
		planenormal.sub(e2, e1);
		planenormal.cross(planenormal, normal);
		planenormal.normalize();

		plane.set(planenormal);
		plane.w = e2.dot(planenormal);
	}

	/**
	 * Finds the closest point(cp) to (v) on a segment (e1,e2).
	 */
	public static void closest_point_on_segment(Vector3d cp, Vector3d v, Vector3d e1, Vector3d e2) {
		Vector3d n = Stack.alloc(Vector3d.class);
		n.sub(e2, e1);
		cp.sub(v, e1);
		double _scalar = cp.dot(n) / n.dot(n);
		if (_scalar < 0.0f) {
			cp = e1;
		}
		else if (_scalar > 1.0f) {
			cp = e2;
		}
		else {
			cp.scaleAdd(_scalar, n, e1);
		}
	}

	/**
	 * Line plane collision.
	 *
	 * @return -0 if the ray never intersects, -1 if the ray collides in front, -2 if the ray collides in back
	 */
	public static int line_plane_collision(Vector4d plane, Vector3d vDir, Vector3d vPoint, Vector3d pout, double[] tparam, double tmin, double tmax) {
		double _dotdir = VectorUtil.dot3(vDir, plane);

		if (Math.abs(_dotdir) < PLANEDIREPSILON) {
			tparam[0] = tmax;
			return 0;
		}

		double _dis = ClipPolygon.distance_point_plane(plane, vPoint);
		int returnvalue = _dis < 0.0f ? 2 : 1;
		tparam[0] = -_dis / _dotdir;

		if (tparam[0] < tmin) {
			returnvalue = 0;
			tparam[0] = tmin;
		}
		else if (tparam[0] > tmax) {
			returnvalue = 0;
			tparam[0] = tmax;
		}
		pout.scaleAdd(tparam[0], vDir, vPoint);
		return returnvalue;
	}

	/**
	 * Find closest points on segments.
	 */
	public static void segment_collision(Vector3d vA1, Vector3d vA2, Vector3d vB1, Vector3d vB2, Vector3d vPointA, Vector3d vPointB) {
		Vector3d AD = Stack.alloc(Vector3d.class);
		AD.sub(vA2, vA1);

		Vector3d BD = Stack.alloc(Vector3d.class);
		BD.sub(vB2, vB1);

		Vector3d N = Stack.alloc(Vector3d.class);
		N.cross(AD, BD);
		double[] tp = new double[] { N.lengthSquared() };

		Vector4d _M = Stack.alloc(Vector4d.class);//plane

		if (tp[0] < BulletGlobals.SIMD_EPSILON)//ARE PARALELE
		{
			// project B over A
			boolean invert_b_order = false;
			_M.x = vB1.dot(AD);
			_M.y = vB2.dot(AD);

			if (_M.x > _M.y) {
				invert_b_order = true;
				//BT_SWAP_NUMBERS(_M[0],_M[1]);
				_M.x = _M.x + _M.y;
				_M.y = _M.x - _M.y;
				_M.x = _M.x - _M.y;
			}
			_M.z = vA1.dot(AD);
			_M.w = vA2.dot(AD);
			// mid points
			N.x = (_M.x + _M.y) * 0.5f;
			N.y = (_M.z + _M.w) * 0.5f;

			if (N.x < N.y) {
				if (_M.y < _M.z) {
					vPointB = invert_b_order ? vB1 : vB2;
					vPointA = vA1;
				}
				else if (_M.y < _M.w) {
					vPointB = invert_b_order ? vB1 : vB2;
					closest_point_on_segment(vPointA, vPointB, vA1, vA2);
				}
				else {
					vPointA = vA2;
					closest_point_on_segment(vPointB, vPointA, vB1, vB2);
				}
			}
			else {
				if (_M.w < _M.x) {
					vPointB = invert_b_order ? vB2 : vB1;
					vPointA = vA2;
				}
				else if (_M.w < _M.y) {
					vPointA = vA2;
					closest_point_on_segment(vPointB, vPointA, vB1, vB2);
				}
				else {
					vPointB = invert_b_order ? vB1 : vB2;
					closest_point_on_segment(vPointA, vPointB, vA1, vA2);
				}
			}
			return;
		}

		N.cross(N, BD);
		_M.set(N.x, N.y, N.z, vB1.dot(N));

		// get point A as the plane collision point
		line_plane_collision(_M, AD, vA1, vPointA, tp, 0f, 1f);

		/*Closest point on segment*/
		vPointB.sub(vPointA, vB1);
		tp[0] = vPointB.dot(BD);
		tp[0] /= BD.dot(BD);
		tp[0] = CLAMP(tp[0], 0.0f, 1.0f);

		vPointB.scaleAdd(tp[0], BD, vB1);
	}

}
