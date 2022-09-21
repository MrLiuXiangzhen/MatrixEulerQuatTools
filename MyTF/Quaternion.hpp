#pragma once
#include "QuadWord.hpp"
#include "Vector3.hpp"

namespace MyTF {

    // 用于欧拉角和四元数之前相互转换
    class Quaternion : public QuadWord {
    public:
        /**@brief No initialization constructor */
        Quaternion() {}


        //		template <typename tfScalar>
        //		explicit Quaternion(const tfScalar *v) : Tuple4<tfScalar>(v) {}
        /**@brief Constructor from scalars */
        Quaternion(const tfScalar &x, const tfScalar &y, const tfScalar &z, const tfScalar &w)
                : QuadWord(x, y, z, w) {}

        /**@brief Axis angle Constructor
         * @param axis The axis which the rotation is around
         * @param angle The magnitude of the rotation around the angle (Radians) */
        Quaternion(const Vector3 &axis, const tfScalar &angle) {
            setRotation(axis, angle);
        }

        /**@brief Constructor from Euler angles
         * @param yaw Angle around Y unless TF_EULER_DEFAULT_ZYX defined then Z
         * @param pitch Angle around X unless TF_EULER_DEFAULT_ZYX defined then Y
         * @param roll Angle around Z unless TF_EULER_DEFAULT_ZYX defined then X */
        Quaternion(const tfScalar &yaw, const tfScalar &pitch, const tfScalar &roll) {
#ifndef TF_EULER_DEFAULT_ZYX
            setEuler(yaw, pitch, roll);
#else
            setRPY(roll, pitch, yaw);
#endif
        }

        /**@brief Set the rotation using axis angle notation
         * @param axis The axis around which to rotate
         * @param angle The magnitude of the rotation in Radians */
        void setRotation(const Vector3 &axis, const tfScalar &angle) {
            tfScalar d = axis.length();
            tfAssert(d != tfScalar(0.0));
            tfScalar s = tfSin(angle * tfScalar(0.5)) / d;
            setValue(axis.x() * s, axis.y() * s, axis.z() * s,
                     tfCos(angle * tfScalar(0.5)));
        }

        /**@brief Set the quaternion using Euler angles
         * @param yaw Angle around Y
         * @param pitch Angle around X
         * @param roll Angle around Z */
        void setEuler(const tfScalar &yaw, const tfScalar &pitch, const tfScalar &roll) {
            tfScalar halfYaw = tfScalar(yaw) * tfScalar(0.5);
            tfScalar halfPitch = tfScalar(pitch) * tfScalar(0.5);
            tfScalar halfRoll = tfScalar(roll) * tfScalar(0.5);
            tfScalar cosYaw = tfCos(halfYaw);
            tfScalar sinYaw = tfSin(halfYaw);
            tfScalar cosPitch = tfCos(halfPitch);
            tfScalar sinPitch = tfSin(halfPitch);
            tfScalar cosRoll = tfCos(halfRoll);
            tfScalar sinRoll = tfSin(halfRoll);
            setValue(cosRoll * sinPitch * cosYaw + sinRoll * cosPitch * sinYaw,
                     cosRoll * cosPitch * sinYaw - sinRoll * sinPitch * cosYaw,
                     sinRoll * cosPitch * cosYaw - cosRoll * sinPitch * sinYaw,
                     cosRoll * cosPitch * cosYaw + sinRoll * sinPitch * sinYaw);
        }

        /**@brief Set the quaternion using fixed axis RPY
         * @param roll Angle around X
         * @param pitch Angle around Y
         * @param yaw Angle around Z*/
        void setRPY(const tfScalar &roll, const tfScalar &pitch, const tfScalar &yaw) {
            tfScalar halfYaw = tfScalar(yaw) * tfScalar(0.5);
            tfScalar halfPitch = tfScalar(pitch) * tfScalar(0.5);
            tfScalar halfRoll = tfScalar(roll) * tfScalar(0.5);
            tfScalar cosYaw = tfCos(halfYaw);
            tfScalar sinYaw = tfSin(halfYaw);
            tfScalar cosPitch = tfCos(halfPitch);
            tfScalar sinPitch = tfSin(halfPitch);
            tfScalar cosRoll = tfCos(halfRoll);
            tfScalar sinRoll = tfSin(halfRoll);
            setValue(sinRoll * cosPitch * cosYaw - cosRoll * sinPitch * sinYaw, //x
                     cosRoll * sinPitch * cosYaw + sinRoll * cosPitch * sinYaw, //y
                     cosRoll * cosPitch * sinYaw - sinRoll * sinPitch * cosYaw, //z
                     cosRoll * cosPitch * cosYaw + sinRoll * sinPitch * sinYaw); //formerly yzx
        }

        /**@brief Set the quaternion using euler angles
         * @param yaw Angle around Z
         * @param pitch Angle around Y
         * @param roll Angle around X */
        void setEulerZYX(const tfScalar &yaw, const tfScalar &pitch, const tfScalar &roll) {
            setRPY(roll, pitch, yaw);
        }

        /**@brief Add two quaternions
         * @param q The quaternion to add to this one */
        Quaternion &operator+=(const Quaternion &q) {
            m_floats[0] += q.x();
            m_floats[1] += q.y();
            m_floats[2] += q.z();
            m_floats[3] += q.m_floats[3];
            return *this;
        }

        /**@brief Sutfract out a quaternion
         * @param q The quaternion to sutfract from this one */
        Quaternion &operator-=(const Quaternion &q) {
            m_floats[0] -= q.x();
            m_floats[1] -= q.y();
            m_floats[2] -= q.z();
            m_floats[3] -= q.m_floats[3];
            return *this;
        }

        /**@brief Scale this quaternion
         * @param s The scalar to scale by */
        Quaternion &operator*=(const tfScalar &s) {
            m_floats[0] *= s;
            m_floats[1] *= s;
            m_floats[2] *= s;
            m_floats[3] *= s;
            return *this;
        }

        /**@brief Multiply this quaternion by q on the right
         * @param q The other quaternion
         * Equivilant to this = this * q */
        Quaternion &operator*=(const Quaternion &q) {
            setValue(m_floats[3] * q.x() + m_floats[0] * q.m_floats[3] + m_floats[1] * q.z() - m_floats[2] * q.y(),
                     m_floats[3] * q.y() + m_floats[1] * q.m_floats[3] + m_floats[2] * q.x() - m_floats[0] * q.z(),
                     m_floats[3] * q.z() + m_floats[2] * q.m_floats[3] + m_floats[0] * q.y() - m_floats[1] * q.x(),
                     m_floats[3] * q.m_floats[3] - m_floats[0] * q.x() - m_floats[1] * q.y() - m_floats[2] * q.z());
            return *this;
        }

        /**@brief Return the dot product between this quaternion and another
         * @param q The other quaternion */
        tfScalar dot(const Quaternion &q) const {
            return m_floats[0] * q.x() + m_floats[1] * q.y() + m_floats[2] * q.z() + m_floats[3] * q.m_floats[3];
        }

        /**@brief Return the length squared of the quaternion */
        tfScalar length2() const {
            return dot(*this);
        }

        /**@brief Return the length of the quaternion */
        tfScalar length() const {
            return tfSqrt(length2());
        }

        /**@brief Normalize the quaternion
         * Such that x^2 + y^2 + z^2 +w^2 = 1 */
        Quaternion &normalize() {
            return *this /= length();
        }

        /**@brief Return a scaled version of this quaternion
         * @param s The scale factor */
        Quaternion
        operator*(const tfScalar &s) const {
            return Quaternion(x() * s, y() * s, z() * s, m_floats[3] * s);
        }


        /**@brief Return an inversely scaled versionof this quaternion
         * @param s The inverse scale factor */
        Quaternion operator/(const tfScalar &s) const {
            tfAssert(s != tfScalar(0.0));
            return *this * (tfScalar(1.0) / s);
        }

        /**@brief Inversely scale this quaternion
         * @param s The scale factor */
        Quaternion &operator/=(const tfScalar &s) {
            tfAssert(s != tfScalar(0.0));
            return *this *= tfScalar(1.0) / s;
        }

        /**@brief Return a normalized version of this quaternion */
        Quaternion normalized() const {
            return *this / length();
        }

        /**@brief Return the ***half*** angle between this quaternion and the other
         * @param q The other quaternion */
        tfScalar angle(const Quaternion &q) const {
            tfScalar s = tfSqrt(length2() * q.length2());
            tfAssert(s != tfScalar(0.0));
            return tfAcos(dot(q) / s);
        }

        /**@brief Return the angle between this quaternion and the other along the shortest path
        * @param q The other quaternion */
        tfScalar angleShortestPath(const Quaternion &q) const {
            tfScalar s = tfSqrt(length2() * q.length2());
            tfAssert(s != tfScalar(0.0));
            if (dot(q) < 0) // Take care of long angle case see http://en.wikipedia.org/wiki/Slerp
                return tfAcos(dot(-q) / s) * tfScalar(2.0);
            else
                return tfAcos(dot(q) / s) * tfScalar(2.0);
        }

        /**@brief Return the angle [0, 2Pi] of rotation represented by this quaternion */
        tfScalar getAngle() const {
            tfScalar s = tfScalar(2.) * tfAcos(m_floats[3]);
            return s;
        }

        /**@brief Return the angle [0, Pi] of rotation represented by this quaternion along the shortest path*/
        tfScalar getAngleShortestPath() const {
            tfScalar s;
            if (m_floats[3] < 0)
                s = tfScalar(2.) * tfAcos(-m_floats[3]);
            else
                s = tfScalar(2.) * tfAcos(m_floats[3]);
            return s;
        }

        /**@brief Return the axis of the rotation represented by this quaternion */
        Vector3 getAxis() const {
            tfScalar s_squared = tfScalar(1.) - tfPow(m_floats[3], tfScalar(2.));
            if (s_squared < tfScalar(10.) * TFSIMD_EPSILON) //Check for divide by zero
                return Vector3(1.0, 0.0, 0.0);  // Arbitrary
            tfScalar s = tfSqrt(s_squared);
            return Vector3(m_floats[0] / s, m_floats[1] / s, m_floats[2] / s);
        }

        /**@brief Return the inverse of this quaternion */
        Quaternion inverse() const {
            return Quaternion(-m_floats[0], -m_floats[1], -m_floats[2], m_floats[3]);
        }

        /**@brief Return the sum of this quaternion and the other
         * @param q2 The other quaternion */
        Quaternion
        operator+(const Quaternion &q2) const {
            const Quaternion &q1 = *this;
            return Quaternion(q1.x() + q2.x(), q1.y() + q2.y(), q1.z() + q2.z(), q1.m_floats[3] + q2.m_floats[3]);
        }

        /**@brief Return the difference between this quaternion and the other
         * @param q2 The other quaternion */
        Quaternion
        operator-(const Quaternion &q2) const {
            const Quaternion &q1 = *this;
            return Quaternion(q1.x() - q2.x(), q1.y() - q2.y(), q1.z() - q2.z(), q1.m_floats[3] - q2.m_floats[3]);
        }

        /**@brief Return the negative of this quaternion
         * This simply negates each element */
        Quaternion operator-() const {
            const Quaternion &q2 = *this;
            return Quaternion(-q2.x(), -q2.y(), -q2.z(), -q2.m_floats[3]);
        }

        /**@todo document this and it's use */
        Quaternion farthest(const Quaternion &qd) const {
            Quaternion diff, sum;
            diff = *this - qd;
            sum = *this + qd;
            if (diff.dot(diff) > sum.dot(sum))
                return qd;
            return (-qd);
        }

        /**@todo document this and it's use */
        Quaternion nearest(const Quaternion &qd) const {
            Quaternion diff, sum;
            diff = *this - qd;
            sum = *this + qd;
            if (diff.dot(diff) < sum.dot(sum))
                return qd;
            return (-qd);
        }


        /**@brief Return the quaternion which is the result of Spherical Linear Interpolation between this and the other quaternion
         * @param q The other quaternion to interpolate with
         * @param t The ratio between this and q to interpolate.  If t = 0 the result is this, if t=1 the result is q.
         * Slerp interpolates assuming constant velocity.  */
        Quaternion slerp(const Quaternion &q, const tfScalar &t) const {
            tfScalar theta = angleShortestPath(q) / tfScalar(2.0);
            if (theta != tfScalar(0.0)) {
                tfScalar d = tfScalar(1.0) / tfSin(theta);
                tfScalar s0 = tfSin((tfScalar(1.0) - t) * theta);
                tfScalar s1 = tfSin(t * theta);
                if (dot(q) < 0) // Take care of long angle case see http://en.wikipedia.org/wiki/Slerp
                    return Quaternion((m_floats[0] * s0 + -q.x() * s1) * d,
                                      (m_floats[1] * s0 + -q.y() * s1) * d,
                                      (m_floats[2] * s0 + -q.z() * s1) * d,
                                      (m_floats[3] * s0 + -q.m_floats[3] * s1) * d);
                else
                    return Quaternion((m_floats[0] * s0 + q.x() * s1) * d,
                                      (m_floats[1] * s0 + q.y() * s1) * d,
                                      (m_floats[2] * s0 + q.z() * s1) * d,
                                      (m_floats[3] * s0 + q.m_floats[3] * s1) * d);

            } else {
                return *this;
            }
        }

        static const Quaternion &getIdentity() {
            static const Quaternion identityQuat(tfScalar(0.), tfScalar(0.), tfScalar(0.), tfScalar(1.));
            return identityQuat;
        }

        const tfScalar &getW() const { return m_floats[3]; }


    };
}
