#pragma once
#include "Scalar.hpp"

namespace MyTF {
#define Vector3Data Vector3DoubleData
#define Vector3DataName "Vector3DoubleData"

/**
 * @class Vector3
 * @brief Vector3 can be used to represent 3D points and vectors.
 * It has an un-used w component to suit 16-byte alignment when Vector3 is stored in containers. This extra component can be used by derived classes (Quaternion?) or by user
 * Ideally, this class should be replaced by a platform optimized SIMD version that keeps the data in registers
 */
    ATTRIBUTE_ALIGNED16(class) Vector3 {
    public:

#if defined (__SPU__) && defined (__CELLOS_LV2__)
        tfScalar	m_floats[4];
public:
     const vec_float4&	get128() const
    {
        return *((const vec_float4*)&m_floats[0]);
    }
public:
#else //__CELLOS_LV2__ __SPU__
#ifdef TF_USE_SSE // _WIN32
        union {
            __m128 mVec128;
            tfScalar	m_floats[4];
        };
            __m128	get128() const
        {
            return mVec128;
        }
            void	set128(__m128 v128)
        {
            mVec128 = v128;
        }
#else
        tfScalar m_floats[4];
#endif
#endif //__CELLOS_LV2__ __SPU__

    public:

        /**@brief No initialization constructor */
        Vector3() {}


        /**@brief Constructor from scalars
         * @param x X value
         * @param y Y value
         * @param z Z value
         */
        Vector3(const tfScalar &x, const tfScalar &y, const tfScalar &z) {
            m_floats[0] = x;
            m_floats[1] = y;
            m_floats[2] = z;
            m_floats[3] = tfScalar(0.);
        }

/**@brief Add a vector to this one 
 * @param The vector to add to this one */
        Vector3 &operator+=(const Vector3 &v) {

            m_floats[0] += v.m_floats[0];
            m_floats[1] += v.m_floats[1];
            m_floats[2] += v.m_floats[2];
            return *this;
        }


        /**@brief Sutfract a vector from this one
         * @param The vector to sutfract */
        Vector3 &operator-=(const Vector3 &v) {
            m_floats[0] -= v.m_floats[0];
            m_floats[1] -= v.m_floats[1];
            m_floats[2] -= v.m_floats[2];
            return *this;
        }

        /**@brief Scale the vector
         * @param s Scale factor */
        Vector3 &operator*=(const tfScalar &s) {
            m_floats[0] *= s;
            m_floats[1] *= s;
            m_floats[2] *= s;
            return *this;
        }

        /**@brief Inversely scale the vector
         * @param s Scale factor to divide by */
        Vector3 &operator/=(const tfScalar &s) {
            tfFullAssert(s != tfScalar(0.0));
            return *this *= tfScalar(1.0) / s;
        }

        /**@brief Return the dot product
         * @param v The other vector in the dot product */
        tfScalar dot(const Vector3 &v) const {
            return m_floats[0] * v.m_floats[0] + m_floats[1] * v.m_floats[1] + m_floats[2] * v.m_floats[2];
        }

        /**@brief Return the length of the vector squared */
        tfScalar length2() const {
            return dot(*this);
        }

        /**@brief Return the length of the vector */
        tfScalar length() const {
            return tfSqrt(length2());
        }

        /**@brief Return the distance squared between the ends of this and another vector
         * This is symantically treating the vector like a point */
        tfScalar distance2(const Vector3 &v) const;

        /**@brief Return the distance between the ends of this and another vector
         * This is symantically treating the vector like a point */
        tfScalar distance(const Vector3 &v) const;

        /**@brief Normalize this vector
         * x^2 + y^2 + z^2 = 1 */
        Vector3 &normalize() {
            return *this /= length();
        }

        /**@brief Return a normalized version of this vector */
        Vector3 normalized() const;

        /**@brief Rotate this vector
         * @param wAxis The axis to rotate about
         * @param angle The angle to rotate by */
        Vector3 rotate(const Vector3 &wAxis, const tfScalar angle) const;

        /**@brief Return the angle between this and another vector
         * @param v The other vector */
        tfScalar angle(const Vector3 &v) const {
            tfScalar s = tfSqrt(length2() * v.length2());
            tfFullAssert(s != tfScalar(0.0));
            return tfAcos(dot(v) / s);
        }

        /**@brief Return a vector will the absolute values of each element */
        Vector3 absolute() const {
            return Vector3(
                    tfFabs(m_floats[0]),
                    tfFabs(m_floats[1]),
                    tfFabs(m_floats[2]));
        }

        /**@brief Return the cross product between this and another vector
         * @param v The other vector */
        Vector3 cross(const Vector3 &v) const {
            return Vector3(
                    m_floats[1] * v.m_floats[2] - m_floats[2] * v.m_floats[1],
                    m_floats[2] * v.m_floats[0] - m_floats[0] * v.m_floats[2],
                    m_floats[0] * v.m_floats[1] - m_floats[1] * v.m_floats[0]);
        }

        tfScalar triple(const Vector3 &v1, const Vector3 &v2) const {
            return m_floats[0] * (v1.m_floats[1] * v2.m_floats[2] - v1.m_floats[2] * v2.m_floats[1]) +
                   m_floats[1] * (v1.m_floats[2] * v2.m_floats[0] - v1.m_floats[0] * v2.m_floats[2]) +
                   m_floats[2] * (v1.m_floats[0] * v2.m_floats[1] - v1.m_floats[1] * v2.m_floats[0]);
        }

        /**@brief Return the axis with the smallest value
         * Note return values are 0,1,2 for x, y, or z */
        int minAxis() const {
            return m_floats[0] < m_floats[1] ? (m_floats[0] < m_floats[2] ? 0 : 2) : (m_floats[1] < m_floats[2] ? 1
                                                                                                                : 2);
        }

        /**@brief Return the axis with the largest value
         * Note return values are 0,1,2 for x, y, or z */
        int maxAxis() const {
            return m_floats[0] < m_floats[1] ? (m_floats[1] < m_floats[2] ? 2 : 1) : (m_floats[0] < m_floats[2] ? 2
                                                                                                                : 0);
        }

        int furthestAxis() const {
            return absolute().minAxis();
        }

        int closestAxis() const {
            return absolute().maxAxis();
        }

        void setInterpolate3(const Vector3 &v0, const Vector3 &v1, tfScalar rt) {
            tfScalar s = tfScalar(1.0) - rt;
            m_floats[0] = s * v0.m_floats[0] + rt * v1.m_floats[0];
            m_floats[1] = s * v0.m_floats[1] + rt * v1.m_floats[1];
            m_floats[2] = s * v0.m_floats[2] + rt * v1.m_floats[2];
            //don't do the unused w component
            //		m_co[3] = s * v0[3] + rt * v1[3];
        }

        /**@brief Return the linear interpolation between this and another vector
         * @param v The other vector
         * @param t The ration of this to v (t = 0 => return this, t=1 => return other) */
        Vector3 lerp(const Vector3 &v, const tfScalar &t) const {
            return Vector3(m_floats[0] + (v.m_floats[0] - m_floats[0]) * t,
                           m_floats[1] + (v.m_floats[1] - m_floats[1]) * t,
                           m_floats[2] + (v.m_floats[2] - m_floats[2]) * t);
        }

        /**@brief Elementwise multiply this vector by the other
         * @param v The other vector */
        Vector3 &operator*=(const Vector3 &v) {
            m_floats[0] *= v.m_floats[0];
            m_floats[1] *= v.m_floats[1];
            m_floats[2] *= v.m_floats[2];
            return *this;
        }

        /**@brief Return the x value */
        const tfScalar &getX() const { return m_floats[0]; }

        /**@brief Return the y value */
        const tfScalar &getY() const { return m_floats[1]; }

        /**@brief Return the z value */
        const tfScalar &getZ() const { return m_floats[2]; }

        /**@brief Set the x value */
        void setX(tfScalar x) { m_floats[0] = x; };

        /**@brief Set the y value */
        void setY(tfScalar y) { m_floats[1] = y; };

        /**@brief Set the z value */
        void setZ(tfScalar z) { m_floats[2] = z; };

        /**@brief Set the w value */
        void setW(tfScalar w) { m_floats[3] = w; };

        /**@brief Return the x value */
        const tfScalar &x() const { return m_floats[0]; }

        /**@brief Return the y value */
        const tfScalar &y() const { return m_floats[1]; }

        /**@brief Return the z value */
        const tfScalar &z() const { return m_floats[2]; }

        /**@brief Return the w value */
        const tfScalar &w() const { return m_floats[3]; }

        // tfScalar&       operator[](int i)       { return (&m_floats[0])[i];	}
        // const tfScalar& operator[](int i) const { return (&m_floats[0])[i]; }
        ///operator tfScalar*() replaces operator[], using implicit conversion. We added operator != and operator == to avoid pointer comparisons.
        operator tfScalar *() { return &m_floats[0]; }

        operator const tfScalar *() const { return &m_floats[0]; }

        bool operator==(const Vector3 &other) const {
            return ((m_floats[3] == other.m_floats[3]) && (m_floats[2] == other.m_floats[2]) &&
                    (m_floats[1] == other.m_floats[1]) && (m_floats[0] == other.m_floats[0]));
        }

        bool operator!=(const Vector3 &other) const {
            return !(*this == other);
        }

        /**@brief Set each element to the max of the current values and the values of another Vector3
      * @param other The other Vector3 to compare with
      */
        void setMax(const Vector3 &other) {
            tfSetMax(m_floats[0], other.m_floats[0]);
            tfSetMax(m_floats[1], other.m_floats[1]);
            tfSetMax(m_floats[2], other.m_floats[2]);
            tfSetMax(m_floats[3], other.w());
        }

        /**@brief Set each element to the min of the current values and the values of another Vector3
         * @param other The other Vector3 to compare with
         */
        void setMin(const Vector3 &other) {
            tfSetMin(m_floats[0], other.m_floats[0]);
            tfSetMin(m_floats[1], other.m_floats[1]);
            tfSetMin(m_floats[2], other.m_floats[2]);
            tfSetMin(m_floats[3], other.w());
        }

        void setValue(const tfScalar &x, const tfScalar &y, const tfScalar &z) {
            m_floats[0] = x;
            m_floats[1] = y;
            m_floats[2] = z;
            m_floats[3] = tfScalar(0.);
        }

        void getSkewSymmetricMatrix(Vector3 *v0, Vector3 *v1, Vector3 *v2) const {
            v0->setValue(0., -z(), y());
            v1->setValue(z(), 0., -x());
            v2->setValue(-y(), x(), 0.);
        }

        void setZero() {
            setValue(tfScalar(0.), tfScalar(0.), tfScalar(0.));
        }

        bool isZero() const {
            return m_floats[0] == tfScalar(0) && m_floats[1] == tfScalar(0) && m_floats[2] == tfScalar(0);
        }

        bool fuzzyZero() const {
            return length2() < TFSIMD_EPSILON;
        }

        void serialize(struct Vector3Data &dataOut) const;

        void deSerialize(const struct Vector3Data &dataIn);

        void serializeFloat(struct Vector3FloatData &dataOut) const;

        void deSerializeFloat(const struct Vector3FloatData &dataIn);

        void serializeDouble(struct Vector3DoubleData &dataOut) const;

        void deSerializeDouble(const struct Vector3DoubleData &dataIn);

    };

/**@brief Return the sum of two vectors (Point symantics)*/
    Vector3
    operator+(const Vector3 &v1, const Vector3 &v2) {
        return Vector3(v1.m_floats[0] + v2.m_floats[0], v1.m_floats[1] + v2.m_floats[1],
                       v1.m_floats[2] + v2.m_floats[2]);
    }

/**@brief Return the elementwise product of two vectors */
    Vector3
    operator*(const Vector3 &v1, const Vector3 &v2) {
        return Vector3(v1.m_floats[0] * v2.m_floats[0], v1.m_floats[1] * v2.m_floats[1],
                       v1.m_floats[2] * v2.m_floats[2]);
    }

/**@brief Return the difference between two vectors */
    Vector3
    operator-(const Vector3 &v1, const Vector3 &v2) {
        return Vector3(v1.m_floats[0] - v2.m_floats[0], v1.m_floats[1] - v2.m_floats[1],
                       v1.m_floats[2] - v2.m_floats[2]);
    }

/**@brief Return the negative of the vector */
    Vector3
    operator-(const Vector3 &v) {
        return Vector3(-v.m_floats[0], -v.m_floats[1], -v.m_floats[2]);
    }

/**@brief Return the vector scaled by s */
    Vector3
    operator*(const Vector3 &v, const tfScalar &s) {
        return Vector3(v.m_floats[0] * s, v.m_floats[1] * s, v.m_floats[2] * s);
    }

/**@brief Return the vector scaled by s */
    Vector3
    operator*(const tfScalar &s, const Vector3 &v) {
        return v * s;
    }

/**@brief Return the vector inversely scaled by s */
    Vector3
    operator/(const Vector3 &v, const tfScalar &s) {
        tfFullAssert(s != tfScalar(0.0));
        return v * (tfScalar(1.0) / s);
    }

/**@brief Return the vector inversely scaled by s */
    Vector3
    operator/(const Vector3 &v1, const Vector3 &v2) {
        return Vector3(v1.m_floats[0] / v2.m_floats[0], v1.m_floats[1] / v2.m_floats[1],
                       v1.m_floats[2] / v2.m_floats[2]);
    }

/**@brief Return the dot product between two vectors */
    tfScalar
    tfDot(const Vector3 &v1, const Vector3 &v2) {
        return v1.dot(v2);
    }


/**@brief Return the distance squared between two vectors */
    tfScalar
    tfDistance2(const Vector3 &v1, const Vector3 &v2) {
        return v1.distance2(v2);
    }


/**@brief Return the distance between two vectors */
    tfScalar
    tfDistance(const Vector3 &v1, const Vector3 &v2) {
        return v1.distance(v2);
    }

/**@brief Return the angle between two vectors */
    tfScalar
    tfAngle(const Vector3 &v1, const Vector3 &v2) {
        return v1.angle(v2);
    }

/**@brief Return the cross product of two vectors */
    Vector3
    tfCross(const Vector3 &v1, const Vector3 &v2) {
        return v1.cross(v2);
    }

    tfScalar
    tfTriple(const Vector3 &v1, const Vector3 &v2, const Vector3 &v3) {
        return v1.triple(v2, v3);
    }

/**@brief Return the linear interpolation between two vectors
 * @param v1 One vector 
 * @param v2 The other vector 
 * @param t The ration of this to v (t = 0 => return v1, t=1 => return v2) */
    Vector3
    lerp(const Vector3 &v1, const Vector3 &v2, const tfScalar &t) {
        return v1.lerp(v2, t);
    }


    tfScalar Vector3::distance2(const Vector3 &v) const {
        return (v - *this).length2();
    }

    tfScalar Vector3::distance(const Vector3 &v) const {
        return (v - *this).length();
    }

    Vector3 Vector3::normalized() const {
        return *this / length();
    }

    Vector3 Vector3::rotate(const Vector3 &wAxis, const tfScalar angle) const {
        // wAxis must be a unit lenght vector

        Vector3 o = wAxis * wAxis.dot(*this);
        Vector3 x = *this - o;
        Vector3 y;

        y = wAxis.cross(*this);

        return (o + x * tfCos(angle) + y * tfSin(angle));
    }

    class tfVector4 : public Vector3 {
    public:

        tfVector4() {}


        tfVector4(const tfScalar &x, const tfScalar &y, const tfScalar &z, const tfScalar &w)
                : Vector3(x, y, z) {
            m_floats[3] = w;
        }


        tfVector4 absolute4() const {
            return tfVector4(
                    tfFabs(m_floats[0]),
                    tfFabs(m_floats[1]),
                    tfFabs(m_floats[2]),
                    tfFabs(m_floats[3]));
        }


        tfScalar getW() const { return m_floats[3]; }


        int maxAxis4() const {
            int maxIndex = -1;
            tfScalar maxVal = tfScalar(-TF_LARGE_FLOAT);
            if (m_floats[0] > maxVal) {
                maxIndex = 0;
                maxVal = m_floats[0];
            }
            if (m_floats[1] > maxVal) {
                maxIndex = 1;
                maxVal = m_floats[1];
            }
            if (m_floats[2] > maxVal) {
                maxIndex = 2;
                maxVal = m_floats[2];
            }
            if (m_floats[3] > maxVal) {
                maxIndex = 3;
                maxVal = m_floats[3];
            }


            return maxIndex;

        }


        int minAxis4() const {
            int minIndex = -1;
            tfScalar minVal = tfScalar(TF_LARGE_FLOAT);
            if (m_floats[0] < minVal) {
                minIndex = 0;
                minVal = m_floats[0];
            }
            if (m_floats[1] < minVal) {
                minIndex = 1;
                minVal = m_floats[1];
            }
            if (m_floats[2] < minVal) {
                minIndex = 2;
                minVal = m_floats[2];
            }
            if (m_floats[3] < minVal) {
                minIndex = 3;
                minVal = m_floats[3];
            }

            return minIndex;

        }


        int closestAxis4() const {
            return absolute4().maxAxis4();
        }




        /**@brief Set x,y,z and zero w
         * @param x Value of x
         * @param y Value of y
         * @param z Value of z
         */


/*		void getValue(tfScalar *m) const 
		{
			m[0] = m_floats[0];
			m[1] = m_floats[1];
			m[2] =m_floats[2];
		}
*/
/**@brief Set the values 
   * @param x Value of x
   * @param y Value of y
   * @param z Value of z
   * @param w Value of w
   */
        void setValue(const tfScalar &x, const tfScalar &y, const tfScalar &z, const tfScalar &w) {
            m_floats[0] = x;
            m_floats[1] = y;
            m_floats[2] = z;
            m_floats[3] = w;
        }


    };


///tfSwapVector3Endian swaps vector endianness, useful for network and cross-platform serialization
    void tfSwapScalarEndian(const tfScalar &sourceVal, tfScalar &destVal) {
        unsigned char *dest = (unsigned char *) &destVal;
        const unsigned char *src = (const unsigned char *) &sourceVal;
        dest[0] = src[7];
        dest[1] = src[6];
        dest[2] = src[5];
        dest[3] = src[4];
        dest[4] = src[3];
        dest[5] = src[2];
        dest[6] = src[1];
        dest[7] = src[0];
    }

///tfSwapVector3Endian swaps vector endianness, useful for network and cross-platform serialization
    void tfSwapVector3Endian(const Vector3 &sourceVec, Vector3 &destVec) {
        for (int i = 0; i < 4; i++) {
            tfSwapScalarEndian(sourceVec[i], destVec[i]);
        }

    }

///tfUnSwapVector3Endian swaps vector endianness, useful for network and cross-platform serialization
    void tfUnSwapVector3Endian(Vector3 &vector) {

        Vector3 swappedVec;
        for (int i = 0; i < 4; i++) {
            tfSwapScalarEndian(vector[i], swappedVec[i]);
        }
        vector = swappedVec;
    }

    void tfPlaneSpace1(const Vector3 &n, Vector3 &p, Vector3 &q) {
        if (tfFabs(n.z()) > TFSIMDSQRT12) {
            // choose p in y-z plane
            tfScalar a = n[1] * n[1] + n[2] * n[2];
            tfScalar k = tfRecipSqrt (a);
            p.setValue(0, -n[2] * k, n[1] * k);
            // set q = n x p
            q.setValue(a * k, -n[0] * p[2], n[0] * p[1]);
        } else {
            // choose p in x-y plane
            tfScalar a = n.x() * n.x() + n.y() * n.y();
            tfScalar k = tfRecipSqrt (a);
            p.setValue(-n.y() * k, n.x() * k, 0);
            // set q = n x p
            q.setValue(-n.z() * p.y(), n.z() * p.x(), a * k);
        }
    }


    struct Vector3FloatData {
        float m_floats[4];
    };

    struct Vector3DoubleData {
        double m_floats[4];

    };

    void Vector3::serializeFloat(struct Vector3FloatData &dataOut) const {
        ///could also do a memcpy, check if it is worth it
        for (int i = 0; i < 4; i++)
            dataOut.m_floats[i] = float(m_floats[i]);
    }

    void Vector3::deSerializeFloat(const struct Vector3FloatData &dataIn) {
        for (int i = 0; i < 4; i++)
            m_floats[i] = tfScalar(dataIn.m_floats[i]);
    }


    void Vector3::serializeDouble(struct Vector3DoubleData &dataOut) const {
        ///could also do a memcpy, check if it is worth it
        for (int i = 0; i < 4; i++)
            dataOut.m_floats[i] = double(m_floats[i]);
    }

    void Vector3::deSerializeDouble(const struct Vector3DoubleData &dataIn) {
        for (int i = 0; i < 4; i++)
            m_floats[i] = tfScalar(dataIn.m_floats[i]);
    }


    void Vector3::serialize(struct Vector3Data &dataOut) const {
        ///could also do a memcpy, check if it is worth it
        for (int i = 0; i < 4; i++)
            dataOut.m_floats[i] = m_floats[i];
    }

    void Vector3::deSerialize(const struct Vector3Data &dataIn) {
        for (int i = 0; i < 4; i++)
            m_floats[i] = dataIn.m_floats[i];
    }

}