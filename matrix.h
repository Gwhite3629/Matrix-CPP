#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <complex>

template <class T> class Vector;

template <class T>
class Matrix
{
public:
    unsigned int rows;
    unsigned int cols;
    bool is_square;

    Matrix(unsigned int, unsigned int);
    Matrix(unsigned int, unsigned int, T, T);
    Matrix(unsigned int, unsigned int, T);
    ~Matrix(void);

    void I(void);

    // Add two matrices
    Matrix operator+(const Matrix&) const;
    void operator+=(const Matrix&);
    // Add constant to matrix
    Matrix operator+(const T) const;
    void operator+=(const T);

    // Subtract two matrices
    Matrix operator-(const Matrix&) const;
    void operator-=(const Matrix&);
    // Subtract constant from matrix
    Matrix operator-(const T) const;
    void operator-=(const T);

    // Multiply matrices element wise
    Matrix operator*(const Matrix&) const;
    void operator*=(const Matrix&);
    // Multiply constant to matrix
    Matrix operator*(const T) const;
    void operator*=(const T);

    // Multiply matrices element wise
    Matrix operator/(const Matrix&) const;
    void operator/=(const Matrix&);
    // Multiply constant to matrix
    Matrix operator/(const T) const;
    void operator/=(const T);

    void operator=(const Matrix&);
    void operator=(const T);

    T sum(void) const;

    // Add row
    Matrix AddVec(const Vector<T>&, unsigned int);
    void AddVecTo(const Vector<T>&, unsigned int);

    T reduce(void);
    T echelon(void);
    void invert(void);
    Matrix augment(const Matrix&) const;
    Matrix slice(unsigned int, unsigned int, unsigned int, unsigned int) const;

    T trace(void) const;
    Matrix copy(void) const;
    void transpose(void);
    Matrix copy_transpose(void) const;
    void conjugate_transpose(void);
    Matrix copy_conjugate_transpose(void) const;
    T inner(const Matrix&) const;
    Matrix outer(const Matrix&) const;
    Matrix submatrix(unsigned int, unsigned int) const;
    T minor(unsigned int, unsigned int) const;
    T cofactor(unsigned int, unsigned int) const;
    T echelon_determinant(void) const;
    T laplace_determinant(void) const;
    Vector<std::complex<T>> eigenvalues(unsigned int) const;

    T condition(unsigned int) const;
    T onenorm(void) const;
    T twonorm(void) const;
    T infnorm(void) const;

    void QR(Matrix *, Matrix *) const;

    T& operator()(unsigned int, unsigned int);
    T operator()(unsigned int, unsigned int) const;

    Vector<T> get_row(unsigned int) const;
    Vector<T> get_col(unsigned int) const;
    void set_row(unsigned int, const Vector<T>&);
    void set_col(unsigned int, const Vector<T>&);

    void swap(unsigned int, unsigned int);

    Matrix del_row(unsigned int) const;
    Matrix del_col(unsigned int) const;

    unsigned int argmax(unsigned int, unsigned int, unsigned int, unsigned int) const;
    bool compare(const Matrix&) const;

    void random(T, T);

    void print(unsigned int) const;

private:
    T *data;

    T get(unsigned int, unsigned int) const;
    void set(unsigned int, unsigned int, T);
};

#endif // _MATRIX_H_