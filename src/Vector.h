/******************************************************************************

    Copyright (C) 2013, Mathieu Tauban
    Email: mathieu.tauban@solvay.com

 *******************************************************************************/
#ifndef VECTOR_H
#define VECTOR_H

#include <stdint.h>
#include <iostream>
#include <sstream>
#include <cmath>


typedef long long int steroint;
//! Implementation of a generic sign function

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
static const double _m_pi = 3.14159265358979323846;

template<typename T, typename T1> inline T imod64(T & n1, T1 & N1) {
    return n1 - floor((double) n1 / N1) * N1;
}

template<typename T, typename T1> inline T1 row_major(T n1, T n2, T n3, T1 N1, T1 N2, T1 N3) {
    return (T1) imod64(n1, N1) * N2 * N3 + imod64(n2, N2) * N3 + imod64(n3, N3);
}

template<typename T, typename T1> inline T1 column_major(T n1, T n2, T n3, T1 N1, T1 N2, T1 N3) {
    return (T1) imod64(n1, N1) + imod64(n2, N2) * N1 + imod64(n3, N3) * N2*N1;
}

template<typename T, typename T1> inline T1 row_major(T n1, T n2, T1 N1, T1 N2) {
    return (T1) imod64(n1, N1) * N2 + imod64(n2, N2);
}

template<typename T, typename T1> inline T1 column_major(T n1, T n2, T1 N1, T1 N2) {
    return (T1) imod64(n1, N1) + imod64(n2, N2) * N1;
}

/*******************************************************************************
 * mtnVector Class for maths
 * Author: Mathieu Tauban
 * 
 ******************************************************************************/
class mtnVector {
public:
    typedef double Scalar;

    inline Scalar & operator() (int i) {
        return data[i];
    }

    inline const Scalar & operator() (int i) const {
        return data[i];
    }

    mtnVector(const mtnVector & ref) {
        data[0] = ref(0);
        data[1] = ref(1);
        data[2] = ref(2);
    }

    mtnVector() {
        data[0] = data[1] = data[2] = 0;
    }

    mtnVector(const Scalar x, const Scalar y, const Scalar z) {

        data[0] = x;
        data[1] = y;
        data[2] = z;
    }

    mtnVector(mtnVector & ref) {
        data[0] = ref.data[0];
        data[1] = ref.data[1];
        data[2] = ref.data[2];
    }

    inline Scalar squaredNorm() {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }

    inline void zeros() {
        data[0] = 0;
        data[1] = 0;
        data[2] = 0;
    }

    inline Scalar norm() {
        Scalar sn = this->squaredNorm();
        return sqrt(sn > 0 ? sn : 0.0);
    }

    inline mtnVector operator*=(const Scalar & a) {
        data[0] *= a;
        data[1] *= a;
        data[2] *= a;
        return *this;
    }

    inline mtnVector operator+=(const mtnVector & a) {
        data[0] += a.data[0];
        data[1] += a.data[1];
        data[2] += a.data[2];
        return *this;
    }

    inline const mtnVector operator-=(const mtnVector & a) {
        data[0] -= a.data[0];
        data[1] -= a.data[1];
        data[2] -= a.data[2];
        return *this;
    }

    inline bool operator==(const mtnVector & b) {
        return ((b.data[0] == data[0]) && (b.data[1] == data[1]) && (b.data[2] == data[2]));
    }

    inline bool operator!=(const mtnVector & b) {
        return ((b.data[0] != data[0]) || (b.data[1] != data[1]) || (b.data[2] != data[2]));
    }
protected:
    Scalar data [3];
};
typedef mtnVector mVector;

template <class OutputStreamT>
inline OutputStreamT& operator<<(OutputStreamT& s, const mtnVector&b) {
    for (int i = 0; i < 3; i++) {
        s << b(i) << ' ';
    }
    return s;
}

inline mtnVector operator+(const mtnVector& a, const mtnVector&b) {
    return mtnVector(a(0) + b(0), a(1) + b(1), a(2) + b(2));
}

inline mtnVector operator-(const mtnVector& a, const mtnVector&b) {
    return mtnVector(a(0) - b(0), a(1) - b(1), a(2) - b(2));
}

inline mtnVector operator*(const mtnVector::Scalar a, const mtnVector &b) {
    return mtnVector(a * b(0), a * b(1), a * b(2));
}

inline mtnVector operator*(const mtnVector &b, const mtnVector::Scalar a) {
    return mtnVector(a * b(0), a * b(1), a * b(2));
}

inline mtnVector::Scalar operator*(const mtnVector &b, const mtnVector&a) {
    return a(0) * b(0) + a(1) * b(1) + a(2) * b(2);
}

inline mtnVector cross(const mtnVector & a, const mtnVector&b) {
    return mtnVector(a(1) * b(2) - a(2) * b(1), a(2) * b(0) - a(0) * b(2), a(0) * b(1) - a(1) * b(0));
}

inline mtnVector operator/(const mtnVector &b, const mtnVector::Scalar a) {
    return mtnVector(b(0) / a, b(1) / a, b(2) / a);
}

class quaternion {
public:
    typedef double Scalar;
    Scalar data[4];

    quaternion() {
    }

    quaternion(mVector & _v) {
        data[0] = 0;
        data[1] = _v(0);
        data[2] = _v(1);
        data[3] = _v(2);
    }

    quaternion(Scalar m, mtnVector _v) {
        data[0] = m;
        data[1] = _v(0);
        data[2] = _v(1);
        data[3] = _v(2);
    }

    quaternion(Scalar _a, Scalar _b, Scalar _c, Scalar _d) {
        data[0] = _a;
        data[1] = _b;
        data[2] = _c;
        data[3] = _d;
    }

    inline quaternion & operator+=(const quaternion & a) {
        data[0] += a.data[0];
        data[1] += a.data[1];
        data[2] += a.data[2];
        data[3] += a.data[3];
        return *this;
    }

    inline Scalar norm() {
        return sqrt(data[0] * data[0] + data[1] * data[1] + data[2] * data[2] + data[3] * data[3]);
    }

    inline quaternion conjugate() {
        return quaternion(data[0], -data[1], -data[2], -data[3]);
    }

    inline Scalar & operator() (int i) {
        return data[i];
    }

    inline const Scalar & operator() (int i) const {
        return data[i];
    }

    //    quaternion(double _a, vector3 _uv) {
    //        double C1 = cos(_a / 2);
    //        double S1 = sin(_a / 2);
    //        a = C1;
    //        b = S1 * _uv.x;
    //        c = S1 * _uv.y;
    //        d = S1 * _uv.z;
    //    }


    // Scalar a, b, c, d;
};
typedef quaternion mQuaternion;

template <class OutputStreamT>
inline OutputStreamT& operator<<(OutputStreamT& s, const quaternion& q) {
    s << q(0) << " " << q(1) << " " << q(2) << " " << q(3) << " ";
    return s;
}

inline quaternion operator+(const quaternion& q1, const quaternion & q2) {
    return quaternion(
            q1(0) + q2(0),
            q1(1) + q2(1),
            q1(2) + q2(2),
            q1(3) + q2(3)
            );
}

inline quaternion quatProd(quaternion& q1, quaternion& q2) {
    return quaternion(
            q1(0) * q2(0) - q1(1) * q2(1) - q1(2) * q2(2) - q1(3) * q2(3),
            q1(0) * q2(1) + q1(1) * q2(0) + q1(2) * q2(3) - q1(3) * q2(2),
            q1(0) * q2(2) + q1(2) * q2(0) + q1(3) * q2(1) - q1(1) * q2(3),
            q1(0) * q2(3) + q1(3) * q2(0) + q1(1) * q2(2) - q1(2) * q2(1)
            );
}

inline quaternion quatProd(const quaternion& q1, const quaternion& q2) {
    return quaternion(
            q1(0) * q2(0) - q1(1) * q2(1) - q1(2) * q2(2) - q1(3) * q2(3),
            q1(0) * q2(1) + q1(1) * q2(0) + q1(2) * q2(3) - q1(3) * q2(2),
            q1(0) * q2(2) + q1(2) * q2(0) + q1(3) * q2(1) - q1(1) * q2(3),
            q1(0) * q2(3) + q1(3) * q2(0) + q1(1) * q2(2) - q1(2) * q2(1)
            );
}

inline quaternion operator*(const quaternion& q1, const quaternion& q2) {
    return quaternion(
            q1(0) * q2(0) - q1(1) * q2(1) - q1(2) * q2(2) - q1(3) * q2(3),
            q1(0) * q2(1) + q1(1) * q2(0) + q1(2) * q2(3) - q1(3) * q2(2),
            q1(0) * q2(2) + q1(2) * q2(0) + q1(3) * q2(1) - q1(1) * q2(3),
            q1(0) * q2(3) + q1(3) * q2(0) + q1(1) * q2(2) - q1(2) * q2(1)
            );
}

inline quaternion operator*(const quaternion::Scalar a, const quaternion& q1) {
    return quaternion(
            a * q1(0),
            a * q1(1),
            a * q1(2),
            a * q1(3)
            );
}

inline quaternion operator*(const quaternion& q1, const quaternion::Scalar a) {
    return quaternion(
            a * q1(0),
            a * q1(1),
            a * q1(2),
            a * q1(3)
            );
}

inline quaternion operator/(const quaternion& q1, const quaternion::Scalar a) {
    return quaternion(
            q1(0) / a,
            q1(1) / a,
            q1(2) / a,
            q1(3) / a
            );
}


//inline quaternion operator+(const vector3&u, const quaternion::DataType a) {
//    return quaternion(a, u.x, u.y, u.z);
//}
//
//inline quaternion operator+(const quaternion::DataType a, const vector3& u) {
//    return u + a;
//}

inline quaternion conjugate(const quaternion&q) {
    return quaternion(q(0), -q(1), -q(2), -q(3));
}


//inline vector3 rotate(const vector3&u, const quaternion&q) {
//    quaternion::DataType t2 = q(0) * q(1);
//    quaternion::DataType t3 = q(0) * q(2);
//    quaternion::DataType t4 = q(0) * q(3);
//    quaternion::DataType t5 = -q(1) * q(1);
//    quaternion::DataType t6 = q(1) * q(2);
//    quaternion::DataType t7 = q(1) * q(3);
//    quaternion::DataType t8 = -q(2) * q(2);
//    quaternion::DataType t9 = q(2) * q(3);
//    quaternion::DataType t10 = -q(3) * q(3);
//
//    return vector3(
//            u.x + 2 * ((t8 + t10) * u.x + (t6 - t4) * u.y + (t3 + t7) * u.z),
//            u.y + 2 * ((t4 + t6) * u.x + (t5 + t10) * u.y + (t9 - t2) * u.z),
//            u.z + 2 * ((t7 - t3) * u.x + (t2 + t9) * u.y + (t5 + t8) * u.z)
//            );
//}
////
////inline quaternion toQuat(const vector3&u) {
////    return quaternion(0, u.x, u.y, u.z);
////}
//
//inline vector3 toVec(const quaternion&q) {
//    return vector3(q(1), q(2), q(3));
//}

//inline mtnVector operator*(const )
class mtnMatrix;
//void gaussj3(mtnMatrix & a);

class mtnMatrix {
    friend class mtnVector;
public:
    typedef double Scalar;
    Scalar data[9];


    mtnMatrix() {
    }

    mtnMatrix(const mtnMatrix & m) {
        for (int i = 0; i < 9; i++) {
            data[i] = m.data[i];
        }
    }

    mtnMatrix(const mtnVector & a, const mtnVector & b, const mtnVector & c) {
        for (int i = 0; i < 3; i++) {
            data[3 * i] = a(i);
            data[3 * i + 1] = b(i);
            data[3 * i + 2] = c(i);
        }
    }

    mtnMatrix(const Scalar a, const Scalar b, const Scalar c,
            const Scalar d, const Scalar e, const Scalar f,
            const Scalar g, const Scalar h, const Scalar i) {

        data[0] = a;
        data[1] = b;
        data[2] = c;
        data[3] = d;
        data[4] = e;
        data[5] = f;
        data[6] = g;
        data[7] = h;
        data[8] = i;

    }
    inline static void badInverse(mtnMatrix & a) {
        mtnMatrix temp;
        temp = a.com();
        temp = temp.transpose();
        temp *= (1 / a.det());
        a = temp;
    }

    inline Scalar det() {
        return data[0] * data[4] * data[8]
                + data[3] * data[7] * data[2]
                + data[6] * data[1] * data[5]
                - data[2] * data[4] * data[6]
                - data[5] * data[7] * data[0]
                - data[8] * data[1] * data[3];
    }

    inline bool operator==(const mtnMatrix & b) {
        bool result = true;
        for (int i = 0; i < 9; i++) result = result && (data[i] == b.data[i]);
        return result;
    }

    inline bool operator!=(const mtnMatrix & b) {
        bool result = false;
        for (int i = 0; i < 9; i++) result = result || (data[i] != b.data[i]);
        return result;
    }

    inline mtnMatrix operator*=(const Scalar a) {
        for (int i = 0; i < 9; i++) data[i] *= a;
        return *this;
    }

    inline mtnMatrix operator+=(const mtnMatrix & a) {
        data[1] += a.data[1];
        data[2] += a.data[2];
        data[3] += a.data[3];
        data[4] += a.data[4];
        data[5] += a.data[5];
        data[6] += a.data[6];
        data[7] += a.data[7];
        data[8] += a.data[8];
        data[0] += a.data[0];

        return *this;
    }

    inline mtnMatrix com() {
        return mtnMatrix(
                data[4] * data[8] - data[7] * data[5], data[5] * data[6] - data[3] * data[8], data[3] * data[7] - data[6] * data[4],
                data[7] * data[2] - data[1] * data[8], data[8] * data[0] - data[2] * data[6], data[6] * data[1] - data[0] * data[7],
                data[1] * data[5] - data[4] * data[2], data[2] * data[3] - data[5] * data[0], data[0] * data[4] - data[3] * data[1]
                );
    }

    inline mtnMatrix transpose() {
        mtnMatrix temp;
        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) temp(i, j) = data[3 * j + i];
        return temp;
    }

    inline Scalar trace() {
        return data[0] + data[4] + data[8];
    }

    inline static const mtnMatrix Identity() {
        return mtnMatrix(1, 0, 0, 0, 1, 0, 0, 0, 1);
    }

    inline static const mtnMatrix Zeros() {
        return mtnMatrix(0, 0, 0, 0, 0, 0, 0, 0, 0);
    }

    inline Scalar& operator() (const int i, const int j) {
        return data[i * 3 + j];
    }

    inline const Scalar& operator() (const int i, const int j) const {
        return data[i * 3 + j];
    }

    inline std::string voigt(const bool eps = false) {
        std::ostringstream str;
        str.precision(6);
        str.width(12);
        if (eps) str << std::scientific << data[0] << "\t" << data[4] << "\t" << data[8] << "\t" << 2 * data[5] << "\t" << 2 * data[2] << "\t" << 2 * data[1];
        else str << std::scientific << data[0] << "\t" << data[4] << "\t" << data[8] << "\t" << data[5] << "\t" << data[2] << "\t" << data[1];
        return str.str();
    }
};

template <class OutputStreamT>
inline OutputStreamT& operator<<(OutputStreamT& s, const mtnMatrix&b) {
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) s << b(i, j) << ' ';
    return s;
}

inline mtnVector operator*(const mtnMatrix & H, const mtnVector& v) {
    return mtnVector(
            H(0, 0) * v(0) + H(0, 1) * v(1) + H(0, 2) * v(2),
            H(1, 0) * v(0) + H(1, 1) * v(1) + H(1, 2) * v(2),
            H(2, 0) * v(0) + H(2, 1) * v(1) + H(2, 2) * v(2));
}

inline mtnMatrix operator*(const mtnMatrix & D, const mtnMatrix& P) {

    mtnMatrix temp = mtnMatrix::Zeros();
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                temp(i, j) += D(i, k) * P(k, j);
            }
        }
    }
    return temp;
}

inline mtnMatrix operator*(const mtnMatrix::Scalar s, const mtnMatrix& P) {

    mtnMatrix temp = mtnMatrix::Zeros();
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            temp(i, j) = P(i, j) * s;
        }
    }
    return temp;
}

inline mtnMatrix operator+(const mtnMatrix & D, const mtnMatrix& P) {

    mtnMatrix temp = mtnMatrix::Zeros();
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            temp(i, j) = D(i, j) + P(i, j);
        }
    }
    return temp;
}

inline mtnMatrix operator-(const mtnMatrix & D, const mtnMatrix& P) {

    mtnMatrix temp = mtnMatrix::Zeros();
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            temp(i, j) = D(i, j) - P(i, j);
        }
    }
    return temp;
}














typedef mtnMatrix mMatrix;

inline mQuaternion unit_Quaternion(const mVector & u, const double angle) {
    double C1 = cos(angle / 2.0);
    double S1 = sin(angle / 2.0);
    return mQuaternion(C1, S1 * u(0), S1 * u(1), S1 * u(2));
}

inline mQuaternion quatrot(mVector & w1, double factor) {
    double NormRot = w1.norm();
    //   std::clog << w1 << "\t" ; 
    if (NormRot > 0) {
        return unit_Quaternion(w1 / NormRot, NormRot / factor);
    } else {
        //     std::clog << mQuaternion(1, 0, 0, 0) << "\n" ; 
        return mQuaternion(1, 0, 0, 0);
    }
}

inline mVector rotate(const mVector&u, const quaternion&q) {
    quaternion::Scalar t2 = q(0) * q(1);
    quaternion::Scalar t3 = q(0) * q(2);
    quaternion::Scalar t4 = q(0) * q(3);
    quaternion::Scalar t5 = -q(1) * q(1);
    quaternion::Scalar t6 = q(1) * q(2);
    quaternion::Scalar t7 = q(1) * q(3);
    quaternion::Scalar t8 = -q(2) * q(2);
    quaternion::Scalar t9 = q(2) * q(3);
    quaternion::Scalar t10 = -q(3) * q(3);

    return mVector(
            u(0) + 2 * ((t8 + t10) * u(0) + (t6 - t4) * u(1) + (t3 + t7) * u(2)),
            u(1) + 2 * ((t4 + t6) * u(0) + (t5 + t10) * u(1) + (t9 - t2) * u(2)),
            u(2) + 2 * ((t7 - t3) * u(0) + (t2 + t9) * u(1) + (t5 + t8) * u(2))
            );
}

template <typename T>
class gMatrix {
public:

    gMatrix(int m, int n = 1, int o = 1) : _m(m), _n(n), _o(o) {
        _data = new T[m * n * o];
    }

    ~gMatrix() {
        delete _data;
    }

    inline T& operator() (const int i, const int j) {
        return _data[_m * i + j];
    }

    inline const T& operator() (const int i, const int j) const {
        return _data[_m * i + j];
    }

    inline T& operator() (const int i) {
        return _data[i];
    }

    inline const T& operator() (const int i) const {
        return _data[i];
    }

    int _m;
    int _n;
    int _o;
protected:
    T * _data;
};

template <typename T>
inline gMatrix<T> operator*(const gMatrix<T> & a, const gMatrix<T> & b) {
    gMatrix<T> temp(a._m, b._n);
    for (int i = 0; i < a._m; i++) {
        for (int j = 0; j < b._n; j++) {
            temp(i, j) = 0;
            for (int k = 0; k < a._n; k++) {
                temp(i, j) += a(i, k) * b(k, j);
            }
        }
    }
    return temp;
}

template <class OutputStreamT, typename T>
inline OutputStreamT& operator<<(OutputStreamT& s, const gMatrix<T>&b) {
    for (int i = 0; i < b._m; i++) {
        s << "";
        for (int j = 0; j < b._n; j++) s << b(i, j) << "";
        s << "\n";
    }
    return s;
}

inline mtnMatrix outter_product(mtnVector & u, mtnVector & v) {
    mtnMatrix uv;
    uv(0, 0) = u(0) * v(0);
    uv(1, 0) = u(1) * v(0);
    uv(2, 0) = u(2) * v(0);
    uv(0, 1) = u(0) * v(1);
    uv(1, 1) = u(1) * v(1);
    uv(2, 1) = u(2) * v(1);
    uv(0, 2) = u(0) * v(2);
    uv(1, 2) = u(1) * v(2);
    uv(2, 2) = u(2) * v(2);
    return uv;
}

inline void matrixtoarray(mtnMatrix & m, double * a) {
    ;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            a[i * 3 + j] = m(i, j);
        }
    }
}

inline void arraytomatrix(double * a, mtnMatrix & m) {
    ;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            m(i, j) = a[i * 3 + j];
        }
    }
}

inline void vectortoarray(mtnVector & m, double * a) {
    ;
    for (int i = 0; i < 3; i++) {
        a[i] = m(i);
    }
}

inline void arraytovector(double * a, mtnVector & m) {
    ;
    for (int i = 0; i < 3; i++) {
        m(i) = a[i];
    }
}

inline void arraytoquaternion(double * a, mQuaternion & m) {
    ;
    for (int i = 0; i < 4; i++) {
        m(i) = a[i];
    }
}

inline void quaterniontoarray(mQuaternion & m, double * a) {
    ;
    for (int i = 0; i < 4; i++) {
        a[i] = m(i);
    }
}

#endif 