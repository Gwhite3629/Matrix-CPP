#include "matrix.h"
#include "vector.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>

template <class T>
Vector<T>::Vector(unsigned int len, bool row_vec)
{
    this->is_row_vec = row_vec;
    this->length = len;
    this->data = (T*)malloc(sizeof(T) * len);
    memset(this->data, 0, sizeof(T)*len);
}

template <class T>
Vector<T>::Vector(unsigned int len, bool row_vec, T start, T end)
{
    time_t t;
    srand((unsigned int) time(&t));
    this->length = len;
    this->is_row_vec = row_vec;
    this->data = (T*)malloc(sizeof(T) * len);

    for (unsigned int i = 0; i < len; i++) {
        this->data[i] = ((T)rand()/(T)(RAND_MAX)) * (end - start) + start;
    }
}

template <class T>
Vector<T>::Vector(unsigned int len, bool row_vec, T v)
{
    this->is_row_vec = row_vec;
    this->length = len;
    this->data = (T*)malloc(sizeof(T) * len);
    memset(this->data, v, sizeof(T) * len);
}

template <class T>
Vector<T>::~Vector(void)
{
    free(this->data);
}

template <class T>
void Vector<T>::clean(void)
{
    memset(this->data, 0, this->length*sizeof(T));
}

template <class T>
Vector<T> Vector<T>::operator+(const Vector& v)
{
    Vector<T> out(this->length, this->is_row_vec);
    for (unsigned int i = 0; i < this->length; i++) {
        out.data[i] = this->data[i] + v.data[i];
    }

    return out;
}

template <class T>
void Vector<T>::operator+=(const Vector& v)
{
    for (unsigned int i = 0; i < this->length; i++)  {
        this->data[i] += v.data[i];
    }
}

template <class T>
Vector<T> Vector<T>::operator+(T v)
{
    Vector<T> out(this->length, this->is_row_vec);
    for (unsigned int i = 0; i < this->length; i++) {
        out.data[i] = this->data[i] + v;
    }

    return out;
}

template <class T>
void Vector<T>::operator+=(T v)
{
    for (unsigned int i = 0; i < this->length; i++) {
        this->data[i] += v;
    }
}

template <class T>
Vector<T> Vector<T>::operator-(const Vector& v)
{
    Vector<T> out(this->length, this->is_row_vec);
    for (unsigned int i = 0; i < this->length; i++)
    {
        out.data[i] = this->data[i] - v.data[i];
    }

    return out;
}

template <class T>
void Vector<T>::operator-=(const Vector& v)
{
    for (unsigned int i = 0; i < this->length; i++)
    {
        this->data[i] -= v.data[i];
    }
}

template <class T>
Vector<T> Vector<T>::operator-(const T v)
{
    Vector<T> out(this->length, this->is_row_vec);
    for (unsigned int i = 0; i < this->length; i++)
    {
        out.data[i] = this->data[i] - v.data[i];
    }

    return out;
}

template <class T>
void Vector<T>::operator-=(const T v)
{
    for (unsigned int i = 0; i < this->length; i++)
    {
        this->data[i] -= v;
    }
}

template <class T>
Vector<T> Vector<T>::operator*(const Vector& v)
{
    Vector<T> out(this->length, this->is_row_vec);
    out = this->copy();
    for (unsigned int i = 0; i < this->length; i++) {
        out.set(i, this->data[i]*v.get(i));
    }
    return out;
}

template <class T>
void Vector<T>::operator*=(const Vector& v)
{
    for (unsigned int i = 0; i < this->length; i++) {
        this->data[i] *= v.get(i);
    }
}

template <class T>
Vector<T> Vector<T>::operator*(const T v)
{
    Vector<T> out(this->length, this->is_row_vec);
    out = this->copy();
    for (unsigned int i = 0; i < this->length; i++) {
        out.set(i, this->data[i]*v);
    }
    return out;
}

template <class T>
void Vector<T>::operator*=(const T v)
{
    for (unsigned int i = 0; i < this->length; i++) {
        this->data[i] *= v;
    }
}

template <class T>
Vector<T> Vector<T>::operator/(const Vector& v)
{
    Vector<T> out(this->length, this->is_row_vec);
    out = this->copy();
    for (unsigned int i = 0; i < this->length; i++) {
        out.set(i, this->data[i]/v.get(i));
    }
    return out;
}

template <class T>
void Vector<T>::operator/=(const Vector& v)
{
    for (unsigned int i = 0; i < this->length; i++) {
        this->data[i] /= v.get(i);
    }
}

template <class T>
Vector<T> Vector<T>::operator/(const T v)
{
    Vector<T> out(this->length, this->is_row_vec);
    out = this->copy();
    for (unsigned int i = 0; i < this->length; i++) {
        out.set(i, this->data[i]/v);
    }
    return out;
}

template <class T>
void Vector<T>::operator/=(const T v)
{
    for (unsigned int i = 0; i < this->length; i++) {
        this->data[i] *= v;
    }
}

template <class T>
void Vector<T>::operator=(const Vector& v)
{
    for (unsigned int i = 0; i < this->length; i++) {
        this->data[i] = v.data[i];
    }
}

template <class T>
void Vector<T>::operator=(const T v)
{
    for (unsigned int i = 0; i < this->length; i++) {
        this->data[i] = v;
    }
}

template <class T>
T Vector<T>::sum(void)
{
    T output = 0;
    for (unsigned int i = 0; i < this->length; i++) {
        output += this->data[i];
    }

    return output;
}

template <class T>
Vector<T> Vector<T>::copy(void)
{
    Vector<T> out(this->length, this->is_row_vec);
    for (unsigned int i = 0; i < this->length; i++) {
        out.set(i, this->get(i));
    }

    return out;
}

template <class T>
void Vector<T>::transpose(void)
{
    this->is_row_vec = !(this->is_row_vec);
}

template <class T>
Vector<T> Vector<T>::copy_transpose(void)
{
    Vector<T> out(this->length, !(this->is_row_vec));
    for (unsigned int i = 0; i < this->length; i++) {
        out.set(i, this->get(i));
    }

    return out;
}

template <class T>
void Vector<T>::linspace(T start, T stop)
{
    T dif = (stop - start + 1)/(this->length);

    for(unsigned int i = 0; i < this->length; i++) {
        this->data[i] = start + dif*i;
    }
}

template <class T>
T Vector<T>::inner(const Vector& v)
{
    T p;
    for (unsigned int i = 0; i < this->length; i++) {
        p += this->data[i] * v.data[i];
    }
    return p;
}

template <class T>
Matrix<T> Vector<T>::outer(const Vector& v)
{
    Matrix<T> M(this->length,v.length);
    for (unsigned int i = 0; i < this->length; i++) {
        for (unsigned int j = 0; j < v.length; j++) {
            M.set(i, j, this->data[i]*v.data[j]);
        }
    }

    return M;
}

template <class T>
T Vector<T>::get(unsigned int i) const
{
    return this->data[i];
}

template <class T>
void Vector<T>::set(unsigned int i, T v)
{
    this->data[i] = v;
}

template <class T>
void Vector<T>::fill(T numbers[])
{
    for (unsigned int i = 0; i < this->length; i++) {
        this->set(i, numbers[i]);
    }
}

template <class T>
void Vector<T>::random(void)
{
    srand(time(NULL));
    for (unsigned int i = 0; i < this->length; i++) {
        this->data.set(i, (T)rand());
    }
}

template <class T>
T Vector<T>::max(void) const
{
    T m = this->get(0);
    for (unsigned int i = 0; i <  this->length; i++) {
        if (abs(this->get(i)) > abs(m))
            m = this->get(i);
    }

    return m;
}

template <class T>
void Vector<T>::print(unsigned int w) const
{
    std::cout << "[";
    for (unsigned int i = 0; i < this->length; i++) {
        std::cout << std::setw(w) << this->data[i];
        if (i != this->length - 1) {
            std::cout << ",";
        }
    }
    std::cout << "]" << std::endl;
}