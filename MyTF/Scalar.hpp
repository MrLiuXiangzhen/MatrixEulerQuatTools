#pragma once
#include <cmath>
#include <cstdlib>
#include <cfloat>

#define ATTRIBUTE_ALIGNED16(a) a
#include <assert.h>
#define tfAssert(x)
#define tfFullAssert(x)

namespace MyTF {

//The tfScalar type abstracts floating point numbers, to easily switch between double and single floating point precision.
    typedef double tfScalar;
//this number could be bigger in double precision
#define TF_LARGE_FLOAT 1e30

    inline tfScalar tfSqrt(tfScalar x) { return sqrt(x); }

    inline tfScalar tfFabs(tfScalar x) { return fabs(x); }

    inline tfScalar tfCos(tfScalar x) { return cos(x); }

    inline tfScalar tfSin(tfScalar x) { return sin(x); }

    inline tfScalar tfTan(tfScalar x) { return tan(x); }

    inline tfScalar tfAcos(tfScalar x) {
        if (x < tfScalar(-1)) x = tfScalar(-1);
        if (x > tfScalar(1)) x = tfScalar(1);
        return acos(x);
    }

    inline tfScalar tfAsin(tfScalar x) {
        if (x < tfScalar(-1)) x = tfScalar(-1);
        if (x > tfScalar(1)) x = tfScalar(1);
        return asin(x);
    }

    inline tfScalar tfAtan(tfScalar x) { return atan(x); }

    inline tfScalar tfAtan2(tfScalar x, tfScalar y) { return atan2(x, y); }

    inline tfScalar tfExp(tfScalar x) { return exp(x); }

    inline tfScalar tfLog(tfScalar x) { return log(x); }

    inline tfScalar tfPow(tfScalar x, tfScalar y) { return pow(x, y); }

    inline tfScalar tfFmod(tfScalar x, tfScalar y) { return fmod(x, y); }

#define TFSIMD_2_PI         tfScalar(6.283185307179586232)
#define TFSIMD_PI           (TFSIMD_2_PI * tfScalar(0.5))
#define TFSIMD_HALF_PI      (TFSIMD_2_PI * tfScalar(0.25))
#define TFSIMD_RADS_PER_DEG (TFSIMD_2_PI / tfScalar(360.0))
#define TFSIMD_DEGS_PER_RAD  (tfScalar(360.0) / TFSIMD_2_PI)
#define TFSIMDSQRT12 tfScalar(0.7071067811865475244008443621048490)

#define tfRecipSqrt(x) ((tfScalar)(tfScalar(1.0)/tfSqrt(tfScalar(x))))        /* reciprocal square root */
#define TFSIMD_EPSILON      DBL_EPSILON
#define TFSIMD_INFINITY     DBL_MAX

    tfScalar tfAtan2Fast(tfScalar y, tfScalar x) {
        tfScalar coeff_1 = TFSIMD_PI / 4.0f;
        tfScalar coeff_2 = 3.0f * coeff_1;
        tfScalar abs_y = tfFabs(y);
        tfScalar angle;
        if (x >= 0.0f) {
            tfScalar r = (x - abs_y) / (x + abs_y);
            angle = coeff_1 - coeff_1 * r;
        } else {
            tfScalar r = (x + abs_y) / (abs_y - x);
            angle = coeff_2 - coeff_1 * r;
        }
        return (y < 0.0f) ? -angle : angle;
    }

    bool tfFuzzyZero(tfScalar x) { return tfFabs(x) < TFSIMD_EPSILON; }

    bool tfEqual(tfScalar a, tfScalar eps) {
        return (((a) <= eps) && !((a) < -eps));
    }

    bool tfGreaterEqual(tfScalar a, tfScalar eps) {
        return (!((a) <= eps));
    }


    int tfIsNegative(tfScalar x) {
        return x < tfScalar(0.0) ? 1 : 0;
    }

    tfScalar tfRadians(tfScalar x) { return x * TFSIMD_RADS_PER_DEG; }

    tfScalar tfDegrees(tfScalar x) { return x * TFSIMD_DEGS_PER_RAD; }

#define TF_DECLARE_HANDLE(name) typedef struct name##__ { int unused; } *name

#ifndef tfFsel

    tfScalar tfFsel(tfScalar a, tfScalar b, tfScalar c) {
        return a >= 0 ? b : c;
    }

#endif
#define tfFsels(a, b, c) (tfScalar)tfFsel(a,b,c)


    bool tfMachineIsLittleEndian() {
        long int i = 1;
        const char *p = (const char *) &i;
        if (p[0] == 1)  // Lowest address contains the least significant byte
            return true;
        else
            return false;
    }


///tfSelect avoids branches, which makes performance much better for consoles like Playstation 3 and XBox 360
///Thanks Phil Knight. See also http://www.cellperformance.com/articles/2006/04/more_techniques_for_eliminatin_1.html
    unsigned tfSelect(unsigned condition, unsigned valueIfConditionNonZero, unsigned valueIfConditionZero) {
        // Set testNz to 0xFFFFFFFF if condition is nonzero, 0x00000000 if condition is zero
        // Rely on positive value or'ed with its negative having sign bit on
        // and zero value or'ed with its negative (which is still zero) having sign bit off
        // Use arithmetic shift right, shifting the sign bit through all 32 bits
        unsigned testNz = (unsigned) (((int) condition | -(int) condition) >> 31);
        unsigned testEqz = ~testNz;
        return ((valueIfConditionNonZero & testNz) | (valueIfConditionZero & testEqz));
    }

    int tfSelect(unsigned condition, int valueIfConditionNonZero, int valueIfConditionZero) {
        unsigned testNz = (unsigned) (((int) condition | -(int) condition) >> 31);
        unsigned testEqz = ~testNz;
        return static_cast<int>((valueIfConditionNonZero & testNz) | (valueIfConditionZero & testEqz));
    }

    float tfSelect(unsigned condition, float valueIfConditionNonZero, float valueIfConditionZero) {
#ifdef TF_HAVE_NATIVE_FSEL
        return (float)tfFsel((tfScalar)condition - tfScalar(1.0f), valueIfConditionNonZero, valueIfConditionZero);
#else
        return (condition != 0) ? valueIfConditionNonZero : valueIfConditionZero;
#endif
    }

    template<typename T>
    void tfSwap(T &a, T &b) {
        T tmp = a;
        a = b;
        b = tmp;
    }


//PCK: endian swapping functions
    unsigned tfSwapEndian(unsigned val) {
        return (((val & 0xff000000) >> 24) | ((val & 0x00ff0000) >> 8) | ((val & 0x0000ff00) << 8) |
                ((val & 0x000000ff) << 24));
    }

    unsigned short tfSwapEndian(unsigned short val) {
        return static_cast<unsigned short>(((val & 0xff00) >> 8) | ((val & 0x00ff) << 8));
    }

    unsigned tfSwapEndian(int val) {
        return tfSwapEndian((unsigned) val);
    }

    unsigned short tfSwapEndian(short val) {
        return tfSwapEndian((unsigned short) val);
    }

///tfSwapFloat uses using char pointers to swap the endianness
////tfSwapFloat/tfSwapDouble will NOT return a float, because the machine might 'correct' invalid floating point values
///Not all values of sign/exponent/mantissa are valid floating point numbers according to IEEE 754.
///When a floating point unit is faced with an invalid value, it may actually change the value, or worse, throw an exception.
///In most systems, running user mode code, you wouldn't get an exception, but instead the hardware/os/runtime will 'fix' the number for you.
///so instead of returning a float/double, we return integer/long long integer
    unsigned int tfSwapEndianFloat(float d) {
        unsigned int a = 0;
        unsigned char *dst = (unsigned char *) &a;
        unsigned char *src = (unsigned char *) &d;

        dst[0] = src[3];
        dst[1] = src[2];
        dst[2] = src[1];
        dst[3] = src[0];
        return a;
    }

// unswap using char pointers
    float tfUnswapEndianFloat(unsigned int a) {
        float d = 0.0f;
        unsigned char *src = (unsigned char *) &a;
        unsigned char *dst = (unsigned char *) &d;

        dst[0] = src[3];
        dst[1] = src[2];
        dst[2] = src[1];
        dst[3] = src[0];

        return d;
    }


// swap using char pointers
    void tfSwapEndianDouble(double d, unsigned char *dst) {
        unsigned char *src = (unsigned char *) &d;

        dst[0] = src[7];
        dst[1] = src[6];
        dst[2] = src[5];
        dst[3] = src[4];
        dst[4] = src[3];
        dst[5] = src[2];
        dst[6] = src[1];
        dst[7] = src[0];

    }

// unswap using char pointers
    double tfUnswapEndianDouble(const unsigned char *src) {
        double d = 0.0;
        unsigned char *dst = (unsigned char *) &d;

        dst[0] = src[7];
        dst[1] = src[6];
        dst[2] = src[5];
        dst[3] = src[4];
        dst[4] = src[3];
        dst[5] = src[2];
        dst[6] = src[1];
        dst[7] = src[0];

        return d;
    }

// returns normalized value in range [-TFSIMD_PI, TFSIMD_PI]
    tfScalar tfNormalizeAngle(tfScalar angleInRadians) {
        angleInRadians = tfFmod(angleInRadians, TFSIMD_2_PI);
        if (angleInRadians < -TFSIMD_PI) {
            return angleInRadians + TFSIMD_2_PI;
        } else if (angleInRadians > TFSIMD_PI) {
            return angleInRadians - TFSIMD_2_PI;
        } else {
            return angleInRadians;
        }
    }

///rudimentary class to provide type info
    struct tfTypedObject {
        tfTypedObject(int objectType)
                : m_objectType(objectType) {
        }

        int m_objectType;

        inline int getObjectType() const {
            return m_objectType;
        }
    };

}