#ifndef TENSOR_HPP
#define TENSOR_HPP

#include <cstdlib>
#include <cstring>

template <typename K> class TensorBase1 {
protected:
    K * m;
    int n[1];
    int b[1];
    int nTotal;
    TensorBase1() : m(0), nTotal(0) {
        memset(this->n, 0, sizeof(this->n));
        memset(this->b, 0, sizeof(this->b));
    }
public:
    TensorBase1(K * m, int n1) : m(m) {
        this->n[0] = n1;
        this->b[0] = 0;
        this->nTotal = n1;
    }
    TensorBase1(K * m, int m1, int n1) : m(m) {
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->nTotal = this->n[0];
    }
    virtual ~TensorBase1() { }
    int size() const {
        return nTotal;
    }
    int size(int i) const {
        return n[i];
    }
    int baseIndex(int i) const {
        return b[i];
    }
    K & operator()(int i1) {
        int k = (i1-b[0]);
        return m[k];
    }
    const K & operator()(int i1) const {
        int k = (i1-b[0]);
        return m[k];
    }
    K * getMemory() {
        return m;
    }
};

template <typename K> class TensorBase2 {
protected:
    K * m;
    int n[2];
    int b[2];
    int nTotal;
    TensorBase2() : m(0), nTotal(0) {
        memset(this->n, 0, sizeof(this->n));
        memset(this->b, 0, sizeof(this->b));
    }
public:
    TensorBase2(K * m, int n1, int n2) : m(m) {
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->nTotal = n1 * n2;
    }
    TensorBase2(K * m, int m1, int n1, int m2, int n2) : m(m) {
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->nTotal = this->n[0] * this->n[1];
    }
    virtual ~TensorBase2() { }
    int size() const {
        return nTotal;
    }
    int size(int i) const {
        return n[i];
    }
    int baseIndex(int i) const {
        return b[i];
    }
    K & operator()(int i1, int i2) {
        int k = (i2-b[1]) +
                   (i1-b[0]) * n[1];
        return m[k];
    }
    const K & operator()(int i1, int i2) const {
        int k = (i2-b[1]) +
                   (i1-b[0]) * n[1];
        return m[k];
    }
    TensorBase1<K> operator()(int i1) {
        int k = (i1-b[0]) * n[1];
        return TensorBase1<K>(&m[k], b[1], b[1]+n[1]-1);
    }
    K * getMemory() {
        return m;
    }
};

template <typename K> class TensorBase3 {
protected:
    K * m;
    int n[3];
    int b[3];
    int nTotal;
    TensorBase3() : m(0), nTotal(0) {
        memset(this->n, 0, sizeof(this->n));
        memset(this->b, 0, sizeof(this->b));
    }
public:
    TensorBase3(K * m, int n1, int n2, int n3) : m(m) {
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->nTotal = n1 * n2 * n3;
    }
    TensorBase3(K * m, int m1, int n1, int m2, int n2, int m3, int n3) : m(m) {
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->nTotal = this->n[0] * this->n[1] * this->n[2];
    }
    virtual ~TensorBase3() { }
    int size() const {
        return nTotal;
    }
    int size(int i) const {
        return n[i];
    }
    int baseIndex(int i) const {
        return b[i];
    }
    K & operator()(int i1, int i2, int i3) {
        int k = (i3-b[2]) +
                   (i2-b[1]) * n[2] +
                   (i1-b[0]) * n[2] * n[1];
        return m[k];
    }
    const K & operator()(int i1, int i2, int i3) const {
        int k = (i3-b[2]) +
                   (i2-b[1]) * n[2] +
                   (i1-b[0]) * n[2] * n[1];
        return m[k];
    }
    TensorBase2<K> operator()(int i1) {
        int k = (i1-b[0]) * n[2] * n[1];
        return TensorBase2<K>(&m[k], b[1], b[1]+n[1]-1, b[2], b[2]+n[2]-1);
    }
    TensorBase1<K> operator()(int i1, int i2) {
        int k = (i2-b[1]) * n[2] +
                   (i1-b[0]) * n[2] * n[1];
        return TensorBase1<K>(&m[k], b[2], b[2]+n[2]-1);
    }
    K * getMemory() {
        return m;
    }
};

template <typename K> class TensorBase4 {
protected:
    K * m;
    int n[4];
    int b[4];
    int nTotal;
    TensorBase4() : m(0), nTotal(0) {
        memset(this->n, 0, sizeof(this->n));
        memset(this->b, 0, sizeof(this->b));
    }
public:
    TensorBase4(K * m, int n1, int n2, int n3, int n4) : m(m) {
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->nTotal = n1 * n2 * n3 * n4;
    }
    TensorBase4(K * m, int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4) : m(m) {
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3];
    }
    virtual ~TensorBase4() { }
    int size() const {
        return nTotal;
    }
    int size(int i) const {
        return n[i];
    }
    int baseIndex(int i) const {
        return b[i];
    }
    K & operator()(int i1, int i2, int i3, int i4) {
        int k = (i4-b[3]) +
                   (i3-b[2]) * n[3] +
                   (i2-b[1]) * n[3] * n[2] +
                   (i1-b[0]) * n[3] * n[2] * n[1];
        return m[k];
    }
    const K & operator()(int i1, int i2, int i3, int i4) const {
        int k = (i4-b[3]) +
                   (i3-b[2]) * n[3] +
                   (i2-b[1]) * n[3] * n[2] +
                   (i1-b[0]) * n[3] * n[2] * n[1];
        return m[k];
    }
    TensorBase3<K> operator()(int i1) {
        int k = (i1-b[0]) * n[3] * n[2] * n[1];
        return TensorBase3<K>(&m[k], b[1], b[1]+n[1]-1, b[2], b[2]+n[2]-1, b[3], b[3]+n[3]-1);
    }
    TensorBase2<K> operator()(int i1, int i2) {
        int k = (i2-b[1]) * n[3] * n[2] +
                   (i1-b[0]) * n[3] * n[2] * n[1];
        return TensorBase2<K>(&m[k], b[2], b[2]+n[2]-1, b[3], b[3]+n[3]-1);
    }
    TensorBase1<K> operator()(int i1, int i2, int i3) {
        int k = (i3-b[2]) * n[3] +
                   (i2-b[1]) * n[3] * n[2] +
                   (i1-b[0]) * n[3] * n[2] * n[1];
        return TensorBase1<K>(&m[k], b[3], b[3]+n[3]-1);
    }
    K * getMemory() {
        return m;
    }
};

template <typename K> class TensorBase5 {
protected:
    K * m;
    int n[5];
    int b[5];
    int nTotal;
    TensorBase5() : m(0), nTotal(0) {
        memset(this->n, 0, sizeof(this->n));
        memset(this->b, 0, sizeof(this->b));
    }
public:
    TensorBase5(K * m, int n1, int n2, int n3, int n4, int n5) : m(m) {
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->n[4] = n5;
        this->b[4] = 0;
        this->nTotal = n1 * n2 * n3 * n4 * n5;
    }
    TensorBase5(K * m, int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4, int m5, int n5) : m(m) {
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->n[4] = n5 - m5 + 1;
        this->b[4] = m5;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3] * this->n[4];
    }
    virtual ~TensorBase5() { }
    int size() const {
        return nTotal;
    }
    int size(int i) const {
        return n[i];
    }
    int baseIndex(int i) const {
        return b[i];
    }
    K & operator()(int i1, int i2, int i3, int i4, int i5) {
        int k = (i5-b[4]) +
                   (i4-b[3]) * n[4] +
                   (i3-b[2]) * n[4] * n[3] +
                   (i2-b[1]) * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[4] * n[3] * n[2] * n[1];
        return m[k];
    }
    const K & operator()(int i1, int i2, int i3, int i4, int i5) const {
        int k = (i5-b[4]) +
                   (i4-b[3]) * n[4] +
                   (i3-b[2]) * n[4] * n[3] +
                   (i2-b[1]) * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[4] * n[3] * n[2] * n[1];
        return m[k];
    }
    TensorBase4<K> operator()(int i1) {
        int k = (i1-b[0]) * n[4] * n[3] * n[2] * n[1];
        return TensorBase4<K>(&m[k], b[1], b[1]+n[1]-1, b[2], b[2]+n[2]-1, b[3], b[3]+n[3]-1, b[4], b[4]+n[4]-1);
    }
    TensorBase3<K> operator()(int i1, int i2) {
        int k = (i2-b[1]) * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[4] * n[3] * n[2] * n[1];
        return TensorBase3<K>(&m[k], b[2], b[2]+n[2]-1, b[3], b[3]+n[3]-1, b[4], b[4]+n[4]-1);
    }
    TensorBase2<K> operator()(int i1, int i2, int i3) {
        int k = (i3-b[2]) * n[4] * n[3] +
                   (i2-b[1]) * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[4] * n[3] * n[2] * n[1];
        return TensorBase2<K>(&m[k], b[3], b[3]+n[3]-1, b[4], b[4]+n[4]-1);
    }
    TensorBase1<K> operator()(int i1, int i2, int i3, int i4) {
        int k = (i4-b[3]) * n[4] +
                   (i3-b[2]) * n[4] * n[3] +
                   (i2-b[1]) * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[4] * n[3] * n[2] * n[1];
        return TensorBase1<K>(&m[k], b[4], b[4]+n[4]-1);
    }
    K * getMemory() {
        return m;
    }
};

template <typename K> class TensorBase6 {
protected:
    K * m;
    int n[6];
    int b[6];
    int nTotal;
    TensorBase6() : m(0), nTotal(0) {
        memset(this->n, 0, sizeof(this->n));
        memset(this->b, 0, sizeof(this->b));
    }
public:
    TensorBase6(K * m, int n1, int n2, int n3, int n4, int n5, int n6) : m(m) {
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->n[4] = n5;
        this->b[4] = 0;
        this->n[5] = n6;
        this->b[5] = 0;
        this->nTotal = n1 * n2 * n3 * n4 * n5 * n6;
    }
    TensorBase6(K * m, int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4, int m5, int n5, int m6, int n6) : m(m) {
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->n[4] = n5 - m5 + 1;
        this->b[4] = m5;
        this->n[5] = n6 - m6 + 1;
        this->b[5] = m6;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3] * this->n[4] * this->n[5];
    }
    virtual ~TensorBase6() { }
    int size() const {
        return nTotal;
    }
    int size(int i) const {
        return n[i];
    }
    int baseIndex(int i) const {
        return b[i];
    }
    K & operator()(int i1, int i2, int i3, int i4, int i5, int i6) {
        int k = (i6-b[5]) +
                   (i5-b[4]) * n[5] +
                   (i4-b[3]) * n[5] * n[4] +
                   (i3-b[2]) * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[5] * n[4] * n[3] * n[2] * n[1];
        return m[k];
    }
    const K & operator()(int i1, int i2, int i3, int i4, int i5, int i6) const {
        int k = (i6-b[5]) +
                   (i5-b[4]) * n[5] +
                   (i4-b[3]) * n[5] * n[4] +
                   (i3-b[2]) * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[5] * n[4] * n[3] * n[2] * n[1];
        return m[k];
    }
    TensorBase5<K> operator()(int i1) {
        int k = (i1-b[0]) * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase5<K>(&m[k], b[1], b[1]+n[1]-1, b[2], b[2]+n[2]-1, b[3], b[3]+n[3]-1, b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1);
    }
    TensorBase4<K> operator()(int i1, int i2) {
        int k = (i2-b[1]) * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase4<K>(&m[k], b[2], b[2]+n[2]-1, b[3], b[3]+n[3]-1, b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1);
    }
    TensorBase3<K> operator()(int i1, int i2, int i3) {
        int k = (i3-b[2]) * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase3<K>(&m[k], b[3], b[3]+n[3]-1, b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1);
    }
    TensorBase2<K> operator()(int i1, int i2, int i3, int i4) {
        int k = (i4-b[3]) * n[5] * n[4] +
                   (i3-b[2]) * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase2<K>(&m[k], b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1);
    }
    TensorBase1<K> operator()(int i1, int i2, int i3, int i4, int i5) {
        int k = (i5-b[4]) * n[5] +
                   (i4-b[3]) * n[5] * n[4] +
                   (i3-b[2]) * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase1<K>(&m[k], b[5], b[5]+n[5]-1);
    }
    K * getMemory() {
        return m;
    }
};

template <typename K> class TensorBase7 {
protected:
    K * m;
    int n[7];
    int b[7];
    int nTotal;
    TensorBase7() : m(0), nTotal(0) {
        memset(this->n, 0, sizeof(this->n));
        memset(this->b, 0, sizeof(this->b));
    }
public:
    TensorBase7(K * m, int n1, int n2, int n3, int n4, int n5, int n6, int n7) : m(m) {
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->n[4] = n5;
        this->b[4] = 0;
        this->n[5] = n6;
        this->b[5] = 0;
        this->n[6] = n7;
        this->b[6] = 0;
        this->nTotal = n1 * n2 * n3 * n4 * n5 * n6 * n7;
    }
    TensorBase7(K * m, int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4, int m5, int n5, int m6, int n6, int m7, int n7) : m(m) {
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->n[4] = n5 - m5 + 1;
        this->b[4] = m5;
        this->n[5] = n6 - m6 + 1;
        this->b[5] = m6;
        this->n[6] = n7 - m7 + 1;
        this->b[6] = m7;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3] * this->n[4] * this->n[5] * this->n[6];
    }
    virtual ~TensorBase7() { }
    int size() const {
        return nTotal;
    }
    int size(int i) const {
        return n[i];
    }
    int baseIndex(int i) const {
        return b[i];
    }
    K & operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7) {
        int k = (i7-b[6]) +
                   (i6-b[5]) * n[6] +
                   (i5-b[4]) * n[6] * n[5] +
                   (i4-b[3]) * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return m[k];
    }
    const K & operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7) const {
        int k = (i7-b[6]) +
                   (i6-b[5]) * n[6] +
                   (i5-b[4]) * n[6] * n[5] +
                   (i4-b[3]) * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return m[k];
    }
    TensorBase6<K> operator()(int i1) {
        int k = (i1-b[0]) * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase6<K>(&m[k], b[1], b[1]+n[1]-1, b[2], b[2]+n[2]-1, b[3], b[3]+n[3]-1, b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1);
    }
    TensorBase5<K> operator()(int i1, int i2) {
        int k = (i2-b[1]) * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase5<K>(&m[k], b[2], b[2]+n[2]-1, b[3], b[3]+n[3]-1, b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1);
    }
    TensorBase4<K> operator()(int i1, int i2, int i3) {
        int k = (i3-b[2]) * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase4<K>(&m[k], b[3], b[3]+n[3]-1, b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1);
    }
    TensorBase3<K> operator()(int i1, int i2, int i3, int i4) {
        int k = (i4-b[3]) * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase3<K>(&m[k], b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1);
    }
    TensorBase2<K> operator()(int i1, int i2, int i3, int i4, int i5) {
        int k = (i5-b[4]) * n[6] * n[5] +
                   (i4-b[3]) * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase2<K>(&m[k], b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1);
    }
    TensorBase1<K> operator()(int i1, int i2, int i3, int i4, int i5, int i6) {
        int k = (i6-b[5]) * n[6] +
                   (i5-b[4]) * n[6] * n[5] +
                   (i4-b[3]) * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase1<K>(&m[k], b[6], b[6]+n[6]-1);
    }
    K * getMemory() {
        return m;
    }
};

template <typename K> class TensorBase8 {
protected:
    K * m;
    int n[8];
    int b[8];
    int nTotal;
    TensorBase8() : m(0), nTotal(0) {
        memset(this->n, 0, sizeof(this->n));
        memset(this->b, 0, sizeof(this->b));
    }
public:
    TensorBase8(K * m, int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8) : m(m) {
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->n[4] = n5;
        this->b[4] = 0;
        this->n[5] = n6;
        this->b[5] = 0;
        this->n[6] = n7;
        this->b[6] = 0;
        this->n[7] = n8;
        this->b[7] = 0;
        this->nTotal = n1 * n2 * n3 * n4 * n5 * n6 * n7 * n8;
    }
    TensorBase8(K * m, int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4, int m5, int n5, int m6, int n6, int m7, int n7, int m8, int n8) : m(m) {
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->n[4] = n5 - m5 + 1;
        this->b[4] = m5;
        this->n[5] = n6 - m6 + 1;
        this->b[5] = m6;
        this->n[6] = n7 - m7 + 1;
        this->b[6] = m7;
        this->n[7] = n8 - m8 + 1;
        this->b[7] = m8;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3] * this->n[4] * this->n[5] * this->n[6] * this->n[7];
    }
    virtual ~TensorBase8() { }
    int size() const {
        return nTotal;
    }
    int size(int i) const {
        return n[i];
    }
    int baseIndex(int i) const {
        return b[i];
    }
    K & operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8) {
        int k = (i8-b[7]) +
                   (i7-b[6]) * n[7] +
                   (i6-b[5]) * n[7] * n[6] +
                   (i5-b[4]) * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return m[k];
    }
    const K & operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8) const {
        int k = (i8-b[7]) +
                   (i7-b[6]) * n[7] +
                   (i6-b[5]) * n[7] * n[6] +
                   (i5-b[4]) * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return m[k];
    }
    TensorBase7<K> operator()(int i1) {
        int k = (i1-b[0]) * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase7<K>(&m[k], b[1], b[1]+n[1]-1, b[2], b[2]+n[2]-1, b[3], b[3]+n[3]-1, b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1);
    }
    TensorBase6<K> operator()(int i1, int i2) {
        int k = (i2-b[1]) * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase6<K>(&m[k], b[2], b[2]+n[2]-1, b[3], b[3]+n[3]-1, b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1);
    }
    TensorBase5<K> operator()(int i1, int i2, int i3) {
        int k = (i3-b[2]) * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase5<K>(&m[k], b[3], b[3]+n[3]-1, b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1);
    }
    TensorBase4<K> operator()(int i1, int i2, int i3, int i4) {
        int k = (i4-b[3]) * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase4<K>(&m[k], b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1);
    }
    TensorBase3<K> operator()(int i1, int i2, int i3, int i4, int i5) {
        int k = (i5-b[4]) * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase3<K>(&m[k], b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1);
    }
    TensorBase2<K> operator()(int i1, int i2, int i3, int i4, int i5, int i6) {
        int k = (i6-b[5]) * n[7] * n[6] +
                   (i5-b[4]) * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase2<K>(&m[k], b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1);
    }
    TensorBase1<K> operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7) {
        int k = (i7-b[6]) * n[7] +
                   (i6-b[5]) * n[7] * n[6] +
                   (i5-b[4]) * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase1<K>(&m[k], b[7], b[7]+n[7]-1);
    }
    K * getMemory() {
        return m;
    }
};

template <typename K> class TensorBase9 {
protected:
    K * m;
    int n[9];
    int b[9];
    int nTotal;
    TensorBase9() : m(0), nTotal(0) {
        memset(this->n, 0, sizeof(this->n));
        memset(this->b, 0, sizeof(this->b));
    }
public:
    TensorBase9(K * m, int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, int n9) : m(m) {
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->n[4] = n5;
        this->b[4] = 0;
        this->n[5] = n6;
        this->b[5] = 0;
        this->n[6] = n7;
        this->b[6] = 0;
        this->n[7] = n8;
        this->b[7] = 0;
        this->n[8] = n9;
        this->b[8] = 0;
        this->nTotal = n1 * n2 * n3 * n4 * n5 * n6 * n7 * n8 * n9;
    }
    TensorBase9(K * m, int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4, int m5, int n5, int m6, int n6, int m7, int n7, int m8, int n8, int m9, int n9) : m(m) {
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->n[4] = n5 - m5 + 1;
        this->b[4] = m5;
        this->n[5] = n6 - m6 + 1;
        this->b[5] = m6;
        this->n[6] = n7 - m7 + 1;
        this->b[6] = m7;
        this->n[7] = n8 - m8 + 1;
        this->b[7] = m8;
        this->n[8] = n9 - m9 + 1;
        this->b[8] = m9;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3] * this->n[4] * this->n[5] * this->n[6] * this->n[7] * this->n[8];
    }
    virtual ~TensorBase9() { }
    int size() const {
        return nTotal;
    }
    int size(int i) const {
        return n[i];
    }
    int baseIndex(int i) const {
        return b[i];
    }
    K & operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9) {
        int k = (i9-b[8]) +
                   (i8-b[7]) * n[8] +
                   (i7-b[6]) * n[8] * n[7] +
                   (i6-b[5]) * n[8] * n[7] * n[6] +
                   (i5-b[4]) * n[8] * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return m[k];
    }
    const K & operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9) const {
        int k = (i9-b[8]) +
                   (i8-b[7]) * n[8] +
                   (i7-b[6]) * n[8] * n[7] +
                   (i6-b[5]) * n[8] * n[7] * n[6] +
                   (i5-b[4]) * n[8] * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return m[k];
    }
    TensorBase8<K> operator()(int i1) {
        int k = (i1-b[0]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase8<K>(&m[k], b[1], b[1]+n[1]-1, b[2], b[2]+n[2]-1, b[3], b[3]+n[3]-1, b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1, b[8], b[8]+n[8]-1);
    }
    TensorBase7<K> operator()(int i1, int i2) {
        int k = (i2-b[1]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase7<K>(&m[k], b[2], b[2]+n[2]-1, b[3], b[3]+n[3]-1, b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1, b[8], b[8]+n[8]-1);
    }
    TensorBase6<K> operator()(int i1, int i2, int i3) {
        int k = (i3-b[2]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase6<K>(&m[k], b[3], b[3]+n[3]-1, b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1, b[8], b[8]+n[8]-1);
    }
    TensorBase5<K> operator()(int i1, int i2, int i3, int i4) {
        int k = (i4-b[3]) * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase5<K>(&m[k], b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1, b[8], b[8]+n[8]-1);
    }
    TensorBase4<K> operator()(int i1, int i2, int i3, int i4, int i5) {
        int k = (i5-b[4]) * n[8] * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase4<K>(&m[k], b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1, b[8], b[8]+n[8]-1);
    }
    TensorBase3<K> operator()(int i1, int i2, int i3, int i4, int i5, int i6) {
        int k = (i6-b[5]) * n[8] * n[7] * n[6] +
                   (i5-b[4]) * n[8] * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase3<K>(&m[k], b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1, b[8], b[8]+n[8]-1);
    }
    TensorBase2<K> operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7) {
        int k = (i7-b[6]) * n[8] * n[7] +
                   (i6-b[5]) * n[8] * n[7] * n[6] +
                   (i5-b[4]) * n[8] * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase2<K>(&m[k], b[7], b[7]+n[7]-1, b[8], b[8]+n[8]-1);
    }
    TensorBase1<K> operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8) {
        int k = (i8-b[7]) * n[8] +
                   (i7-b[6]) * n[8] * n[7] +
                   (i6-b[5]) * n[8] * n[7] * n[6] +
                   (i5-b[4]) * n[8] * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase1<K>(&m[k], b[8], b[8]+n[8]-1);
    }
    K * getMemory() {
        return m;
    }
};

template <typename K> class TensorBase10 {
protected:
    K * m;
    int n[10];
    int b[10];
    int nTotal;
    TensorBase10() : m(0), nTotal(0) {
        memset(this->n, 0, sizeof(this->n));
        memset(this->b, 0, sizeof(this->b));
    }
public:
    TensorBase10(K * m, int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, int n9, int n10) : m(m) {
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->n[4] = n5;
        this->b[4] = 0;
        this->n[5] = n6;
        this->b[5] = 0;
        this->n[6] = n7;
        this->b[6] = 0;
        this->n[7] = n8;
        this->b[7] = 0;
        this->n[8] = n9;
        this->b[8] = 0;
        this->n[9] = n10;
        this->b[9] = 0;
        this->nTotal = n1 * n2 * n3 * n4 * n5 * n6 * n7 * n8 * n9 * n10;
    }
    TensorBase10(K * m, int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4, int m5, int n5, int m6, int n6, int m7, int n7, int m8, int n8, int m9, int n9, int m10, int n10) : m(m) {
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->n[4] = n5 - m5 + 1;
        this->b[4] = m5;
        this->n[5] = n6 - m6 + 1;
        this->b[5] = m6;
        this->n[6] = n7 - m7 + 1;
        this->b[6] = m7;
        this->n[7] = n8 - m8 + 1;
        this->b[7] = m8;
        this->n[8] = n9 - m9 + 1;
        this->b[8] = m9;
        this->n[9] = n10 - m10 + 1;
        this->b[9] = m10;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3] * this->n[4] * this->n[5] * this->n[6] * this->n[7] * this->n[8] * this->n[9];
    }
    virtual ~TensorBase10() { }
    int size() const {
        return nTotal;
    }
    int size(int i) const {
        return n[i];
    }
    int baseIndex(int i) const {
        return b[i];
    }
    K & operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, int i10) {
        int k = (i10-b[9]) +
                   (i9-b[8]) * n[9] +
                   (i8-b[7]) * n[9] * n[8] +
                   (i7-b[6]) * n[9] * n[8] * n[7] +
                   (i6-b[5]) * n[9] * n[8] * n[7] * n[6] +
                   (i5-b[4]) * n[9] * n[8] * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return m[k];
    }
    const K & operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, int i10) const {
        int k = (i10-b[9]) +
                   (i9-b[8]) * n[9] +
                   (i8-b[7]) * n[9] * n[8] +
                   (i7-b[6]) * n[9] * n[8] * n[7] +
                   (i6-b[5]) * n[9] * n[8] * n[7] * n[6] +
                   (i5-b[4]) * n[9] * n[8] * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return m[k];
    }
    TensorBase9<K> operator()(int i1) {
        int k = (i1-b[0]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase9<K>(&m[k], b[1], b[1]+n[1]-1, b[2], b[2]+n[2]-1, b[3], b[3]+n[3]-1, b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1, b[8], b[8]+n[8]-1, b[9], b[9]+n[9]-1);
    }
    TensorBase8<K> operator()(int i1, int i2) {
        int k = (i2-b[1]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase8<K>(&m[k], b[2], b[2]+n[2]-1, b[3], b[3]+n[3]-1, b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1, b[8], b[8]+n[8]-1, b[9], b[9]+n[9]-1);
    }
    TensorBase7<K> operator()(int i1, int i2, int i3) {
        int k = (i3-b[2]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase7<K>(&m[k], b[3], b[3]+n[3]-1, b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1, b[8], b[8]+n[8]-1, b[9], b[9]+n[9]-1);
    }
    TensorBase6<K> operator()(int i1, int i2, int i3, int i4) {
        int k = (i4-b[3]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase6<K>(&m[k], b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1, b[8], b[8]+n[8]-1, b[9], b[9]+n[9]-1);
    }
    TensorBase5<K> operator()(int i1, int i2, int i3, int i4, int i5) {
        int k = (i5-b[4]) * n[9] * n[8] * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase5<K>(&m[k], b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1, b[8], b[8]+n[8]-1, b[9], b[9]+n[9]-1);
    }
    TensorBase4<K> operator()(int i1, int i2, int i3, int i4, int i5, int i6) {
        int k = (i6-b[5]) * n[9] * n[8] * n[7] * n[6] +
                   (i5-b[4]) * n[9] * n[8] * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase4<K>(&m[k], b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1, b[8], b[8]+n[8]-1, b[9], b[9]+n[9]-1);
    }
    TensorBase3<K> operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7) {
        int k = (i7-b[6]) * n[9] * n[8] * n[7] +
                   (i6-b[5]) * n[9] * n[8] * n[7] * n[6] +
                   (i5-b[4]) * n[9] * n[8] * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase3<K>(&m[k], b[7], b[7]+n[7]-1, b[8], b[8]+n[8]-1, b[9], b[9]+n[9]-1);
    }
    TensorBase2<K> operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8) {
        int k = (i8-b[7]) * n[9] * n[8] +
                   (i7-b[6]) * n[9] * n[8] * n[7] +
                   (i6-b[5]) * n[9] * n[8] * n[7] * n[6] +
                   (i5-b[4]) * n[9] * n[8] * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase2<K>(&m[k], b[8], b[8]+n[8]-1, b[9], b[9]+n[9]-1);
    }
    TensorBase1<K> operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9) {
        int k = (i9-b[8]) * n[9] +
                   (i8-b[7]) * n[9] * n[8] +
                   (i7-b[6]) * n[9] * n[8] * n[7] +
                   (i6-b[5]) * n[9] * n[8] * n[7] * n[6] +
                   (i5-b[4]) * n[9] * n[8] * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase1<K>(&m[k], b[9], b[9]+n[9]-1);
    }
    K * getMemory() {
        return m;
    }
};

template <typename K> class TensorBase11 {
protected:
    K * m;
    int n[11];
    int b[11];
    int nTotal;
    TensorBase11() : m(0), nTotal(0) {
        memset(this->n, 0, sizeof(this->n));
        memset(this->b, 0, sizeof(this->b));
    }
public:
    TensorBase11(K * m, int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, int n9, int n10, int n11) : m(m) {
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->n[4] = n5;
        this->b[4] = 0;
        this->n[5] = n6;
        this->b[5] = 0;
        this->n[6] = n7;
        this->b[6] = 0;
        this->n[7] = n8;
        this->b[7] = 0;
        this->n[8] = n9;
        this->b[8] = 0;
        this->n[9] = n10;
        this->b[9] = 0;
        this->n[10] = n11;
        this->b[10] = 0;
        this->nTotal = n1 * n2 * n3 * n4 * n5 * n6 * n7 * n8 * n9 * n10 * n11;
    }
    TensorBase11(K * m, int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4, int m5, int n5, int m6, int n6, int m7, int n7, int m8, int n8, int m9, int n9, int m10, int n10, int m11, int n11) : m(m) {
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->n[4] = n5 - m5 + 1;
        this->b[4] = m5;
        this->n[5] = n6 - m6 + 1;
        this->b[5] = m6;
        this->n[6] = n7 - m7 + 1;
        this->b[6] = m7;
        this->n[7] = n8 - m8 + 1;
        this->b[7] = m8;
        this->n[8] = n9 - m9 + 1;
        this->b[8] = m9;
        this->n[9] = n10 - m10 + 1;
        this->b[9] = m10;
        this->n[10] = n11 - m11 + 1;
        this->b[10] = m11;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3] * this->n[4] * this->n[5] * this->n[6] * this->n[7] * this->n[8] * this->n[9] * this->n[10];
    }
    virtual ~TensorBase11() { }
    int size() const {
        return nTotal;
    }
    int size(int i) const {
        return n[i];
    }
    int baseIndex(int i) const {
        return b[i];
    }
    K & operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, int i10, int i11) {
        int k = (i11-b[10]) +
                   (i10-b[9]) * n[10] +
                   (i9-b[8]) * n[10] * n[9] +
                   (i8-b[7]) * n[10] * n[9] * n[8] +
                   (i7-b[6]) * n[10] * n[9] * n[8] * n[7] +
                   (i6-b[5]) * n[10] * n[9] * n[8] * n[7] * n[6] +
                   (i5-b[4]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return m[k];
    }
    const K & operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, int i10, int i11) const {
        int k = (i11-b[10]) +
                   (i10-b[9]) * n[10] +
                   (i9-b[8]) * n[10] * n[9] +
                   (i8-b[7]) * n[10] * n[9] * n[8] +
                   (i7-b[6]) * n[10] * n[9] * n[8] * n[7] +
                   (i6-b[5]) * n[10] * n[9] * n[8] * n[7] * n[6] +
                   (i5-b[4]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return m[k];
    }
    TensorBase10<K> operator()(int i1) {
        int k = (i1-b[0]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase10<K>(&m[k], b[1], b[1]+n[1]-1, b[2], b[2]+n[2]-1, b[3], b[3]+n[3]-1, b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1, b[8], b[8]+n[8]-1, b[9], b[9]+n[9]-1, b[10], b[10]+n[10]-1);
    }
    TensorBase9<K> operator()(int i1, int i2) {
        int k = (i2-b[1]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase9<K>(&m[k], b[2], b[2]+n[2]-1, b[3], b[3]+n[3]-1, b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1, b[8], b[8]+n[8]-1, b[9], b[9]+n[9]-1, b[10], b[10]+n[10]-1);
    }
    TensorBase8<K> operator()(int i1, int i2, int i3) {
        int k = (i3-b[2]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase8<K>(&m[k], b[3], b[3]+n[3]-1, b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1, b[8], b[8]+n[8]-1, b[9], b[9]+n[9]-1, b[10], b[10]+n[10]-1);
    }
    TensorBase7<K> operator()(int i1, int i2, int i3, int i4) {
        int k = (i4-b[3]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase7<K>(&m[k], b[4], b[4]+n[4]-1, b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1, b[8], b[8]+n[8]-1, b[9], b[9]+n[9]-1, b[10], b[10]+n[10]-1);
    }
    TensorBase6<K> operator()(int i1, int i2, int i3, int i4, int i5) {
        int k = (i5-b[4]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase6<K>(&m[k], b[5], b[5]+n[5]-1, b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1, b[8], b[8]+n[8]-1, b[9], b[9]+n[9]-1, b[10], b[10]+n[10]-1);
    }
    TensorBase5<K> operator()(int i1, int i2, int i3, int i4, int i5, int i6) {
        int k = (i6-b[5]) * n[10] * n[9] * n[8] * n[7] * n[6] +
                   (i5-b[4]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase5<K>(&m[k], b[6], b[6]+n[6]-1, b[7], b[7]+n[7]-1, b[8], b[8]+n[8]-1, b[9], b[9]+n[9]-1, b[10], b[10]+n[10]-1);
    }
    TensorBase4<K> operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7) {
        int k = (i7-b[6]) * n[10] * n[9] * n[8] * n[7] +
                   (i6-b[5]) * n[10] * n[9] * n[8] * n[7] * n[6] +
                   (i5-b[4]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase4<K>(&m[k], b[7], b[7]+n[7]-1, b[8], b[8]+n[8]-1, b[9], b[9]+n[9]-1, b[10], b[10]+n[10]-1);
    }
    TensorBase3<K> operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8) {
        int k = (i8-b[7]) * n[10] * n[9] * n[8] +
                   (i7-b[6]) * n[10] * n[9] * n[8] * n[7] +
                   (i6-b[5]) * n[10] * n[9] * n[8] * n[7] * n[6] +
                   (i5-b[4]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase3<K>(&m[k], b[8], b[8]+n[8]-1, b[9], b[9]+n[9]-1, b[10], b[10]+n[10]-1);
    }
    TensorBase2<K> operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9) {
        int k = (i9-b[8]) * n[10] * n[9] +
                   (i8-b[7]) * n[10] * n[9] * n[8] +
                   (i7-b[6]) * n[10] * n[9] * n[8] * n[7] +
                   (i6-b[5]) * n[10] * n[9] * n[8] * n[7] * n[6] +
                   (i5-b[4]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase2<K>(&m[k], b[9], b[9]+n[9]-1, b[10], b[10]+n[10]-1);
    }
    TensorBase1<K> operator()(int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, int i10) {
        int k = (i10-b[9]) * n[10] +
                   (i9-b[8]) * n[10] * n[9] +
                   (i8-b[7]) * n[10] * n[9] * n[8] +
                   (i7-b[6]) * n[10] * n[9] * n[8] * n[7] +
                   (i6-b[5]) * n[10] * n[9] * n[8] * n[7] * n[6] +
                   (i5-b[4]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] +
                   (i4-b[3]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] +
                   (i3-b[2]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] +
                   (i2-b[1]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] +
                   (i1-b[0]) * n[10] * n[9] * n[8] * n[7] * n[6] * n[5] * n[4] * n[3] * n[2] * n[1];
        return TensorBase1<K>(&m[k], b[10], b[10]+n[10]-1);
    }
    K * getMemory() {
        return m;
    }
};

template <typename K> class Tensor1 : public TensorBase1<K> {
public:
    Tensor1() {
        memset(this->n, 0, sizeof(this->n));
        memset(this->b, 0, sizeof(this->b));
    }
    Tensor1(const TensorBase1<K> & t) {
        memcpy(this->n, t.n, sizeof(this->n));
        memcpy(this->b, t.b, sizeof(this->b));
        this->nTotal = t.nTotal;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memcpy(this->m, t.m, this->nTotal * sizeof(K));
        }
    }
    Tensor1(const Tensor1<K> & t) {
        memcpy(this->n, t.n, sizeof(this->n));
        memcpy(this->b, t.b, sizeof(this->b));
        this->nTotal = t.nTotal;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memcpy(this->m, t.m, this->nTotal * sizeof(K));
        }
    }
    Tensor1(int n1) {
        this->n[0] = n1;
        this->b[0] = 0;
        this->nTotal = n1;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    Tensor1(int m1, int n1) {
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->nTotal = this->n[0];
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    ~Tensor1() {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
    }
    void resize(int n1) {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
        this->n[0] = n1;
        this->b[0] = 0;
        this->nTotal = n1;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    void resize(int m1, int n1) {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->nTotal = this->n[0];
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
};

template <typename K> class Tensor2 : public TensorBase2<K> {
public:
    Tensor2() {
        memset(this->n, 0, sizeof(this->n));
        memset(this->b, 0, sizeof(this->b));
    }
    Tensor2(const TensorBase2<K> & t) {
        memcpy(this->n, t.n, sizeof(this->n));
        memcpy(this->b, t.b, sizeof(this->b));
        this->nTotal = t.nTotal;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memcpy(this->m, t.m, this->nTotal * sizeof(K));
        }
    }
    Tensor2(const Tensor2<K> & t) {
        memcpy(this->n, t.n, sizeof(this->n));
        memcpy(this->b, t.b, sizeof(this->b));
        this->nTotal = t.nTotal;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memcpy(this->m, t.m, this->nTotal * sizeof(K));
        }
    }
    Tensor2(int n1, int n2) {
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->nTotal = n1 * n2;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    Tensor2(int m1, int n1, int m2, int n2) {
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->nTotal = this->n[0] * this->n[1];
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    ~Tensor2() {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
    }
    void resize(int n1, int n2) {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->nTotal = n1 * n2;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    void resize(int m1, int n1, int m2, int n2) {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->nTotal = this->n[0] * this->n[1];
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
};

template <typename K> class Tensor3 : public TensorBase3<K> {
public:
    Tensor3() {
        memset(this->n, 0, sizeof(this->n));
        memset(this->b, 0, sizeof(this->b));
    }
    Tensor3(const TensorBase3<K> & t) {
        memcpy(this->n, t.n, sizeof(this->n));
        memcpy(this->b, t.b, sizeof(this->b));
        this->nTotal = t.nTotal;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memcpy(this->m, t.m, this->nTotal * sizeof(K));
        }
    }
    Tensor3(const Tensor3<K> & t) {
        memcpy(this->n, t.n, sizeof(this->n));
        memcpy(this->b, t.b, sizeof(this->b));
        this->nTotal = t.nTotal;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memcpy(this->m, t.m, this->nTotal * sizeof(K));
        }
    }
    Tensor3(int n1, int n2, int n3) {
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->nTotal = n1 * n2 * n3;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    Tensor3(int m1, int n1, int m2, int n2, int m3, int n3) {
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->nTotal = this->n[0] * this->n[1] * this->n[2];
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    ~Tensor3() {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
    }
    void resize(int n1, int n2, int n3) {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->nTotal = n1 * n2 * n3;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    void resize(int m1, int n1, int m2, int n2, int m3, int n3) {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->nTotal = this->n[0] * this->n[1] * this->n[2];
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
};

template <typename K> class Tensor4 : public TensorBase4<K> {
public:
    Tensor4() {
        memset(this->n, 0, sizeof(this->n));
        memset(this->b, 0, sizeof(this->b));
    }
    Tensor4(const TensorBase4<K> & t) {
        memcpy(this->n, t.n, sizeof(this->n));
        memcpy(this->b, t.b, sizeof(this->b));
        this->nTotal = t.nTotal;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memcpy(this->m, t.m, this->nTotal * sizeof(K));
        }
    }
    Tensor4(const Tensor4<K> & t) {
        memcpy(this->n, t.n, sizeof(this->n));
        memcpy(this->b, t.b, sizeof(this->b));
        this->nTotal = t.nTotal;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memcpy(this->m, t.m, this->nTotal * sizeof(K));
        }
    }
    Tensor4(int n1, int n2, int n3, int n4) {
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->nTotal = n1 * n2 * n3 * n4;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    Tensor4(int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4) {
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3];
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    ~Tensor4() {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
    }
    void resize(int n1, int n2, int n3, int n4) {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->nTotal = n1 * n2 * n3 * n4;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    void resize(int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4) {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3];
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
};

template <typename K> class Tensor5 : public TensorBase5<K> {
public:
    Tensor5() {
        memset(this->n, 0, sizeof(this->n));
        memset(this->b, 0, sizeof(this->b));
    }
    Tensor5(const TensorBase5<K> & t) {
        memcpy(this->n, t.n, sizeof(this->n));
        memcpy(this->b, t.b, sizeof(this->b));
        this->nTotal = t.nTotal;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memcpy(this->m, t.m, this->nTotal * sizeof(K));
        }
    }
    Tensor5(const Tensor5<K> & t) {
        memcpy(this->n, t.n, sizeof(this->n));
        memcpy(this->b, t.b, sizeof(this->b));
        this->nTotal = t.nTotal;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memcpy(this->m, t.m, this->nTotal * sizeof(K));
        }
    }
    Tensor5(int n1, int n2, int n3, int n4, int n5) {
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->n[4] = n5;
        this->b[4] = 0;
        this->nTotal = n1 * n2 * n3 * n4 * n5;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    Tensor5(int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4, int m5, int n5) {
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->n[4] = n5 - m5 + 1;
        this->b[4] = m5;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3] * this->n[4];
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    ~Tensor5() {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
    }
    void resize(int n1, int n2, int n3, int n4, int n5) {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->n[4] = n5;
        this->b[4] = 0;
        this->nTotal = n1 * n2 * n3 * n4 * n5;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    void resize(int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4, int m5, int n5) {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->n[4] = n5 - m5 + 1;
        this->b[4] = m5;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3] * this->n[4];
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
};

template <typename K> class Tensor6 : public TensorBase6<K> {
public:
    Tensor6() {
        memset(this->n, 0, sizeof(this->n));
        memset(this->b, 0, sizeof(this->b));
    }
    Tensor6(const TensorBase6<K> & t) {
        memcpy(this->n, t.n, sizeof(this->n));
        memcpy(this->b, t.b, sizeof(this->b));
        this->nTotal = t.nTotal;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memcpy(this->m, t.m, this->nTotal * sizeof(K));
        }
    }
    Tensor6(const Tensor6<K> & t) {
        memcpy(this->n, t.n, sizeof(this->n));
        memcpy(this->b, t.b, sizeof(this->b));
        this->nTotal = t.nTotal;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memcpy(this->m, t.m, this->nTotal * sizeof(K));
        }
    }
    Tensor6(int n1, int n2, int n3, int n4, int n5, int n6) {
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->n[4] = n5;
        this->b[4] = 0;
        this->n[5] = n6;
        this->b[5] = 0;
        this->nTotal = n1 * n2 * n3 * n4 * n5 * n6;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    Tensor6(int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4, int m5, int n5, int m6, int n6) {
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->n[4] = n5 - m5 + 1;
        this->b[4] = m5;
        this->n[5] = n6 - m6 + 1;
        this->b[5] = m6;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3] * this->n[4] * this->n[5];
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    ~Tensor6() {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
    }
    void resize(int n1, int n2, int n3, int n4, int n5, int n6) {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->n[4] = n5;
        this->b[4] = 0;
        this->n[5] = n6;
        this->b[5] = 0;
        this->nTotal = n1 * n2 * n3 * n4 * n5 * n6;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    void resize(int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4, int m5, int n5, int m6, int n6) {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->n[4] = n5 - m5 + 1;
        this->b[4] = m5;
        this->n[5] = n6 - m6 + 1;
        this->b[5] = m6;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3] * this->n[4] * this->n[5];
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
};

template <typename K> class Tensor7 : public TensorBase7<K> {
public:
    Tensor7() {
        memset(this->n, 0, sizeof(this->n));
        memset(this->b, 0, sizeof(this->b));
    }
    Tensor7(const TensorBase7<K> & t) {
        memcpy(this->n, t.n, sizeof(this->n));
        memcpy(this->b, t.b, sizeof(this->b));
        this->nTotal = t.nTotal;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memcpy(this->m, t.m, this->nTotal * sizeof(K));
        }
    }
    Tensor7(const Tensor7<K> & t) {
        memcpy(this->n, t.n, sizeof(this->n));
        memcpy(this->b, t.b, sizeof(this->b));
        this->nTotal = t.nTotal;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memcpy(this->m, t.m, this->nTotal * sizeof(K));
        }
    }
    Tensor7(int n1, int n2, int n3, int n4, int n5, int n6, int n7) {
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->n[4] = n5;
        this->b[4] = 0;
        this->n[5] = n6;
        this->b[5] = 0;
        this->n[6] = n7;
        this->b[6] = 0;
        this->nTotal = n1 * n2 * n3 * n4 * n5 * n6 * n7;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    Tensor7(int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4, int m5, int n5, int m6, int n6, int m7, int n7) {
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->n[4] = n5 - m5 + 1;
        this->b[4] = m5;
        this->n[5] = n6 - m6 + 1;
        this->b[5] = m6;
        this->n[6] = n7 - m7 + 1;
        this->b[6] = m7;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3] * this->n[4] * this->n[5] * this->n[6];
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    ~Tensor7() {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
    }
    void resize(int n1, int n2, int n3, int n4, int n5, int n6, int n7) {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->n[4] = n5;
        this->b[4] = 0;
        this->n[5] = n6;
        this->b[5] = 0;
        this->n[6] = n7;
        this->b[6] = 0;
        this->nTotal = n1 * n2 * n3 * n4 * n5 * n6 * n7;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    void resize(int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4, int m5, int n5, int m6, int n6, int m7, int n7) {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->n[4] = n5 - m5 + 1;
        this->b[4] = m5;
        this->n[5] = n6 - m6 + 1;
        this->b[5] = m6;
        this->n[6] = n7 - m7 + 1;
        this->b[6] = m7;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3] * this->n[4] * this->n[5] * this->n[6];
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
};

template <typename K> class Tensor8 : public TensorBase8<K> {
public:
    Tensor8() {
        memset(this->n, 0, sizeof(this->n));
        memset(this->b, 0, sizeof(this->b));
    }
    Tensor8(const TensorBase8<K> & t) {
        memcpy(this->n, t.n, sizeof(this->n));
        memcpy(this->b, t.b, sizeof(this->b));
        this->nTotal = t.nTotal;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memcpy(this->m, t.m, this->nTotal * sizeof(K));
        }
    }
    Tensor8(const Tensor8<K> & t) {
        memcpy(this->n, t.n, sizeof(this->n));
        memcpy(this->b, t.b, sizeof(this->b));
        this->nTotal = t.nTotal;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memcpy(this->m, t.m, this->nTotal * sizeof(K));
        }
    }
    Tensor8(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8) {
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->n[4] = n5;
        this->b[4] = 0;
        this->n[5] = n6;
        this->b[5] = 0;
        this->n[6] = n7;
        this->b[6] = 0;
        this->n[7] = n8;
        this->b[7] = 0;
        this->nTotal = n1 * n2 * n3 * n4 * n5 * n6 * n7 * n8;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    Tensor8(int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4, int m5, int n5, int m6, int n6, int m7, int n7, int m8, int n8) {
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->n[4] = n5 - m5 + 1;
        this->b[4] = m5;
        this->n[5] = n6 - m6 + 1;
        this->b[5] = m6;
        this->n[6] = n7 - m7 + 1;
        this->b[6] = m7;
        this->n[7] = n8 - m8 + 1;
        this->b[7] = m8;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3] * this->n[4] * this->n[5] * this->n[6] * this->n[7];
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    ~Tensor8() {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
    }
    void resize(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8) {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->n[4] = n5;
        this->b[4] = 0;
        this->n[5] = n6;
        this->b[5] = 0;
        this->n[6] = n7;
        this->b[6] = 0;
        this->n[7] = n8;
        this->b[7] = 0;
        this->nTotal = n1 * n2 * n3 * n4 * n5 * n6 * n7 * n8;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    void resize(int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4, int m5, int n5, int m6, int n6, int m7, int n7, int m8, int n8) {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->n[4] = n5 - m5 + 1;
        this->b[4] = m5;
        this->n[5] = n6 - m6 + 1;
        this->b[5] = m6;
        this->n[6] = n7 - m7 + 1;
        this->b[6] = m7;
        this->n[7] = n8 - m8 + 1;
        this->b[7] = m8;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3] * this->n[4] * this->n[5] * this->n[6] * this->n[7];
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
};

template <typename K> class Tensor9 : public TensorBase9<K> {
public:
    Tensor9() {
        memset(this->n, 0, sizeof(this->n));
        memset(this->b, 0, sizeof(this->b));
    }
    Tensor9(const TensorBase9<K> & t) {
        memcpy(this->n, t.n, sizeof(this->n));
        memcpy(this->b, t.b, sizeof(this->b));
        this->nTotal = t.nTotal;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memcpy(this->m, t.m, this->nTotal * sizeof(K));
        }
    }
    Tensor9(const Tensor9<K> & t) {
        memcpy(this->n, t.n, sizeof(this->n));
        memcpy(this->b, t.b, sizeof(this->b));
        this->nTotal = t.nTotal;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memcpy(this->m, t.m, this->nTotal * sizeof(K));
        }
    }
    Tensor9(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, int n9) {
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->n[4] = n5;
        this->b[4] = 0;
        this->n[5] = n6;
        this->b[5] = 0;
        this->n[6] = n7;
        this->b[6] = 0;
        this->n[7] = n8;
        this->b[7] = 0;
        this->n[8] = n9;
        this->b[8] = 0;
        this->nTotal = n1 * n2 * n3 * n4 * n5 * n6 * n7 * n8 * n9;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    Tensor9(int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4, int m5, int n5, int m6, int n6, int m7, int n7, int m8, int n8, int m9, int n9) {
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->n[4] = n5 - m5 + 1;
        this->b[4] = m5;
        this->n[5] = n6 - m6 + 1;
        this->b[5] = m6;
        this->n[6] = n7 - m7 + 1;
        this->b[6] = m7;
        this->n[7] = n8 - m8 + 1;
        this->b[7] = m8;
        this->n[8] = n9 - m9 + 1;
        this->b[8] = m9;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3] * this->n[4] * this->n[5] * this->n[6] * this->n[7] * this->n[8];
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    ~Tensor9() {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
    }
    void resize(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, int n9) {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->n[4] = n5;
        this->b[4] = 0;
        this->n[5] = n6;
        this->b[5] = 0;
        this->n[6] = n7;
        this->b[6] = 0;
        this->n[7] = n8;
        this->b[7] = 0;
        this->n[8] = n9;
        this->b[8] = 0;
        this->nTotal = n1 * n2 * n3 * n4 * n5 * n6 * n7 * n8 * n9;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    void resize(int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4, int m5, int n5, int m6, int n6, int m7, int n7, int m8, int n8, int m9, int n9) {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->n[4] = n5 - m5 + 1;
        this->b[4] = m5;
        this->n[5] = n6 - m6 + 1;
        this->b[5] = m6;
        this->n[6] = n7 - m7 + 1;
        this->b[6] = m7;
        this->n[7] = n8 - m8 + 1;
        this->b[7] = m8;
        this->n[8] = n9 - m9 + 1;
        this->b[8] = m9;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3] * this->n[4] * this->n[5] * this->n[6] * this->n[7] * this->n[8];
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
};

template <typename K> class Tensor10 : public TensorBase10<K> {
public:
    Tensor10() {
        memset(this->n, 0, sizeof(this->n));
        memset(this->b, 0, sizeof(this->b));
    }
    Tensor10(const TensorBase10<K> & t) {
        memcpy(this->n, t.n, sizeof(this->n));
        memcpy(this->b, t.b, sizeof(this->b));
        this->nTotal = t.nTotal;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memcpy(this->m, t.m, this->nTotal * sizeof(K));
        }
    }
    Tensor10(const Tensor10<K> & t) {
        memcpy(this->n, t.n, sizeof(this->n));
        memcpy(this->b, t.b, sizeof(this->b));
        this->nTotal = t.nTotal;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memcpy(this->m, t.m, this->nTotal * sizeof(K));
        }
    }
    Tensor10(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, int n9, int n10) {
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->n[4] = n5;
        this->b[4] = 0;
        this->n[5] = n6;
        this->b[5] = 0;
        this->n[6] = n7;
        this->b[6] = 0;
        this->n[7] = n8;
        this->b[7] = 0;
        this->n[8] = n9;
        this->b[8] = 0;
        this->n[9] = n10;
        this->b[9] = 0;
        this->nTotal = n1 * n2 * n3 * n4 * n5 * n6 * n7 * n8 * n9 * n10;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    Tensor10(int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4, int m5, int n5, int m6, int n6, int m7, int n7, int m8, int n8, int m9, int n9, int m10, int n10) {
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->n[4] = n5 - m5 + 1;
        this->b[4] = m5;
        this->n[5] = n6 - m6 + 1;
        this->b[5] = m6;
        this->n[6] = n7 - m7 + 1;
        this->b[6] = m7;
        this->n[7] = n8 - m8 + 1;
        this->b[7] = m8;
        this->n[8] = n9 - m9 + 1;
        this->b[8] = m9;
        this->n[9] = n10 - m10 + 1;
        this->b[9] = m10;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3] * this->n[4] * this->n[5] * this->n[6] * this->n[7] * this->n[8] * this->n[9];
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    ~Tensor10() {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
    }
    void resize(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, int n9, int n10) {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->n[4] = n5;
        this->b[4] = 0;
        this->n[5] = n6;
        this->b[5] = 0;
        this->n[6] = n7;
        this->b[6] = 0;
        this->n[7] = n8;
        this->b[7] = 0;
        this->n[8] = n9;
        this->b[8] = 0;
        this->n[9] = n10;
        this->b[9] = 0;
        this->nTotal = n1 * n2 * n3 * n4 * n5 * n6 * n7 * n8 * n9 * n10;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    void resize(int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4, int m5, int n5, int m6, int n6, int m7, int n7, int m8, int n8, int m9, int n9, int m10, int n10) {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->n[4] = n5 - m5 + 1;
        this->b[4] = m5;
        this->n[5] = n6 - m6 + 1;
        this->b[5] = m6;
        this->n[6] = n7 - m7 + 1;
        this->b[6] = m7;
        this->n[7] = n8 - m8 + 1;
        this->b[7] = m8;
        this->n[8] = n9 - m9 + 1;
        this->b[8] = m9;
        this->n[9] = n10 - m10 + 1;
        this->b[9] = m10;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3] * this->n[4] * this->n[5] * this->n[6] * this->n[7] * this->n[8] * this->n[9];
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
};

template <typename K> class Tensor11 : public TensorBase11<K> {
public:
    Tensor11() {
        memset(this->n, 0, sizeof(this->n));
        memset(this->b, 0, sizeof(this->b));
    }
    Tensor11(const TensorBase11<K> & t) {
        memcpy(this->n, t.n, sizeof(this->n));
        memcpy(this->b, t.b, sizeof(this->b));
        this->nTotal = t.nTotal;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memcpy(this->m, t.m, this->nTotal * sizeof(K));
        }
    }
    Tensor11(const Tensor11<K> & t) {
        memcpy(this->n, t.n, sizeof(this->n));
        memcpy(this->b, t.b, sizeof(this->b));
        this->nTotal = t.nTotal;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memcpy(this->m, t.m, this->nTotal * sizeof(K));
        }
    }
    Tensor11(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, int n9, int n10, int n11) {
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->n[4] = n5;
        this->b[4] = 0;
        this->n[5] = n6;
        this->b[5] = 0;
        this->n[6] = n7;
        this->b[6] = 0;
        this->n[7] = n8;
        this->b[7] = 0;
        this->n[8] = n9;
        this->b[8] = 0;
        this->n[9] = n10;
        this->b[9] = 0;
        this->n[10] = n11;
        this->b[10] = 0;
        this->nTotal = n1 * n2 * n3 * n4 * n5 * n6 * n7 * n8 * n9 * n10 * n11;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    Tensor11(int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4, int m5, int n5, int m6, int n6, int m7, int n7, int m8, int n8, int m9, int n9, int m10, int n10, int m11, int n11) {
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->n[4] = n5 - m5 + 1;
        this->b[4] = m5;
        this->n[5] = n6 - m6 + 1;
        this->b[5] = m6;
        this->n[6] = n7 - m7 + 1;
        this->b[6] = m7;
        this->n[7] = n8 - m8 + 1;
        this->b[7] = m8;
        this->n[8] = n9 - m9 + 1;
        this->b[8] = m9;
        this->n[9] = n10 - m10 + 1;
        this->b[9] = m10;
        this->n[10] = n11 - m11 + 1;
        this->b[10] = m11;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3] * this->n[4] * this->n[5] * this->n[6] * this->n[7] * this->n[8] * this->n[9] * this->n[10];
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    ~Tensor11() {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
    }
    void resize(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, int n9, int n10, int n11) {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
        this->n[0] = n1;
        this->b[0] = 0;
        this->n[1] = n2;
        this->b[1] = 0;
        this->n[2] = n3;
        this->b[2] = 0;
        this->n[3] = n4;
        this->b[3] = 0;
        this->n[4] = n5;
        this->b[4] = 0;
        this->n[5] = n6;
        this->b[5] = 0;
        this->n[6] = n7;
        this->b[6] = 0;
        this->n[7] = n8;
        this->b[7] = 0;
        this->n[8] = n9;
        this->b[8] = 0;
        this->n[9] = n10;
        this->b[9] = 0;
        this->n[10] = n11;
        this->b[10] = 0;
        this->nTotal = n1 * n2 * n3 * n4 * n5 * n6 * n7 * n8 * n9 * n10 * n11;
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
    void resize(int m1, int n1, int m2, int n2, int m3, int n3, int m4, int n4, int m5, int n5, int m6, int n6, int m7, int n7, int m8, int n8, int m9, int n9, int m10, int n10, int m11, int n11) {
        if (this->m != 0) {
            delete[] this->m;
            this->m = 0;
        }
        this->n[0] = n1 - m1 + 1;
        this->b[0] = m1;
        this->n[1] = n2 - m2 + 1;
        this->b[1] = m2;
        this->n[2] = n3 - m3 + 1;
        this->b[2] = m3;
        this->n[3] = n4 - m4 + 1;
        this->b[3] = m4;
        this->n[4] = n5 - m5 + 1;
        this->b[4] = m5;
        this->n[5] = n6 - m6 + 1;
        this->b[5] = m6;
        this->n[6] = n7 - m7 + 1;
        this->b[6] = m7;
        this->n[7] = n8 - m8 + 1;
        this->b[7] = m8;
        this->n[8] = n9 - m9 + 1;
        this->b[8] = m9;
        this->n[9] = n10 - m10 + 1;
        this->b[9] = m10;
        this->n[10] = n11 - m11 + 1;
        this->b[10] = m11;
        this->nTotal = this->n[0] * this->n[1] * this->n[2] * this->n[3] * this->n[4] * this->n[5] * this->n[6] * this->n[7] * this->n[8] * this->n[9] * this->n[10];
        this->m = new K[this->nTotal];
        if (this->m == 0) {
            memset(this->n, 0, sizeof(this->n));
            memset(this->b, 0, sizeof(this->b));
            this->nTotal = 0;
        } else {
            memset(this->m, 0, this->nTotal * sizeof(K));
        }
    }
};

#endif

