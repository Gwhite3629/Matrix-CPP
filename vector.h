#ifndef _VECTOR_H_
#define _VECTOR_H_

template <class T> class Matrix;

template <class T>
class Vector
{
public:
    unsigned int length;
    bool is_row_vec;

    Vector(unsigned int, bool);
    Vector(unsigned int, bool, T, T);
    Vector(unsigned int, bool, T);
    ~Vector(void);
    void clean(void);

    // Add two vectors
    Vector operator+(const Vector&);
    void operator+=(const Vector&);
    // Add constant to vector
    Vector operator+(const T);
    void operator+=(const T);

    // Subtract two vectors
    Vector operator-(const Vector&);
    void operator-=(const Vector&);
    // Subtract constant from vector
    Vector operator-(const T);
    void operator-=(const T);

    // Multiply two vectors elementwise
    Vector operator*(const Vector&);
    void operator*=(const Vector&);
    // Multiply vector by constant
    Vector operator*(const T);
    void operator*=(const T);

    // Divide two vectors elementwise
    Vector operator/(const Vector&);
    void operator/=(const Vector&);
    // Divide vector by constant
    Vector operator/(const T);
    void operator/=(const T);

    void operator=(const Vector&);
    void operator=(const T);

    T sum(void);

    void linspace(T, T);

    Vector copy(void);
    void transpose(void);
    Vector copy_transpose(void);
    T inner(const Vector&);
    Matrix<T> outer(const Vector&);

    void unit(void);
    T magnitude(void) const;

    T& operator()(unsigned int);
    T operator()(unsigned int) const;

    T get(unsigned int) const;
    void set(unsigned int, T);
    void fill(T[]);
    void random(void);
    T max(void) const;

    void print(unsigned int) const;

private:
    T *data;
};

#endif // _VECTOR_H_