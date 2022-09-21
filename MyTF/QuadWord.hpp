#pragma once
#include "Scalar.hpp"

template <class T>
const T& tfMin(const T& a, const T& b) {
  return a < b ? a : b ;
}

template <class T>
const T& tfMax(const T& a, const T& b) {
  return  a > b ? a : b;
}


template <class T>
void tfSetMin(T& a, const T& b) {
    if (b < a) 
	{
		a = b;
	}
}

template <class T>
void tfSetMax(T& a, const T& b) {
    if (a < b) 
	{
		a = b;
	}
}


namespace MyTF {

    class QuadWord {
    protected:
#if defined (__SPU__) && defined (__CELLOS_LV2__)
        union {
            vec_float4 mVec128;
            tfScalar	m_floats[4];
        };
    public:
        vec_float4	get128() const
        {
            return mVec128;
        }
    protected:
#else //__CELLOS_LV2__ __SPU__
        tfScalar m_floats[4];
#endif //__CELLOS_LV2__ __SPU__

    public:


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

        //TFSIMD_FORCE_INLINE tfScalar&       operator[](int i)       { return (&m_floats[0])[i];	}
        //TFSIMD_FORCE_INLINE const tfScalar& operator[](int i) const { return (&m_floats[0])[i]; }
        ///operator tfScalar*() replaces operator[], using implicit conversion. We added operator != and operator == to avoid pointer comparisons.
        operator tfScalar *() { return &m_floats[0]; }

        operator const tfScalar *() const { return &m_floats[0]; }

        bool operator==(const QuadWord &other) const {
            return ((m_floats[3] == other.m_floats[3]) && (m_floats[2] == other.m_floats[2]) &&
                    (m_floats[1] == other.m_floats[1]) && (m_floats[0] == other.m_floats[0]));
        }

        bool operator!=(const QuadWord &other) const {
            return !(*this == other);
        }

        /**@brief Set x,y,z and zero w
         * @param x Value of x
         * @param y Value of y
         * @param z Value of z
         */
        void setValue(const tfScalar &x, const tfScalar &y, const tfScalar &z) {
            m_floats[0] = x;
            m_floats[1] = y;
            m_floats[2] = z;
            m_floats[3] = 0.f;
        }

/*		void getValue(tfScalar *m) const 
		{
			m[0] = m_floats[0];
			m[1] = m_floats[1];
			m[2] = m_floats[2];
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

        /**@brief No initialization constructor */
        QuadWord()
        //	:m_floats[0](tfScalar(0.)),m_floats[1](tfScalar(0.)),m_floats[2](tfScalar(0.)),m_floats[3](tfScalar(0.))
        {
        }

        /**@brief Three argument constructor (zeros w)
         * @param x Value of x
         * @param y Value of y
         * @param z Value of z
         */
        QuadWord(const tfScalar &x, const tfScalar &y, const tfScalar &z) {
            m_floats[0] = x, m_floats[1] = y, m_floats[2] = z, m_floats[3] = 0.0f;
        }

/**@brief Initializing constructor
   * @param x Value of x
   * @param y Value of y
   * @param z Value of z
   * @param w Value of w
   */
        QuadWord(const tfScalar &x, const tfScalar &y, const tfScalar &z, const tfScalar &w) {
            m_floats[0] = x, m_floats[1] = y, m_floats[2] = z, m_floats[3] = w;
        }

        /**@brief Set each element to the max of the current values and the values of another QuadWord
         * @param other The other QuadWord to compare with
         */
        void setMax(const QuadWord &other) {
            tfSetMax(m_floats[0], other.m_floats[0]);
            tfSetMax(m_floats[1], other.m_floats[1]);
            tfSetMax(m_floats[2], other.m_floats[2]);
            tfSetMax(m_floats[3], other.m_floats[3]);
        }

        /**@brief Set each element to the min of the current values and the values of another QuadWord
         * @param other The other QuadWord to compare with
         */
        void setMin(const QuadWord &other) {
            tfSetMin(m_floats[0], other.m_floats[0]);
            tfSetMin(m_floats[1], other.m_floats[1]);
            tfSetMin(m_floats[2], other.m_floats[2]);
            tfSetMin(m_floats[3], other.m_floats[3]);
        }
    };

}