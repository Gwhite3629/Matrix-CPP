#include "matrix.h"
#include "vector.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <cstring>
#include <complex>

template <class T>
Matrix<T>::Matrix(unsigned int R, unsigned int C)
{
    this->rows = R;
    this->cols = C;
    if (R == C) {
        this->is_square = 1;
    } else {
        this->is_square = 0;
    }
    this->data = (T*)malloc(sizeof(T) * R * C);
    memset(this->data, 0, sizeof(T) * R * C);
}

template <class T>
Matrix<T>::Matrix(unsigned int R, unsigned int C, T start, T end)
{
    time_t t;
    srand((unsigned int) time(&t));
    this->rows = R;
    this->cols = C;
    if (R == C) {
        this->is_square = 1;
    } else {
        this->is_square = 0;
    }
    this->data = (T*)malloc(sizeof(T) * R * C);
    
    for (unsigned int i = 0; i < R*C; i++) {
        this->data[i] = ((T)rand()/(T)(RAND_MAX)) * (end - start) + start;
    }
}

template <class T>
Matrix<T>::Matrix(unsigned int R, unsigned int C, T v)
{
    this->rows = R;
    this->cols = C;
    if (R == C) {
        this->is_square = 1;
    } else {
        this->is_square = 0;
    }
    this->data = (T*)malloc(sizeof(T) * R * C);
    
    for (unsigned int i = 0; i < R*C; i++) {
        this->data[i] = v;
    }
}

template <class T>
Matrix<T>::~Matrix(void)
{
    if (this->data)
        free(this->data);
    this->data = NULL;
}

template <class T>
void Matrix<T>::I(void)
{
    assert(is_square);
    memset(this->data, 0, this->cols*this->rows*sizeof(T));
    for (unsigned int i = 0; i < this->rows; i++) {
        this->set(i, i, 1);
    }
}

template <class T>
Matrix<T> Matrix<T>::operator+(const Matrix& M) const
{
    Matrix<T> out(this->rows, this->cols);
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            out.set(i, j, (this->get(i, j) + M.get(i, j)));
        }
    }

    return out;
}

template <class T>
void Matrix<T>::operator+=(const Matrix& M)
{
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            this->set(i, j, (this->get(i, j) + M.get(i, j)));
        }
    }
}

template <class T>
Matrix<T> Matrix<T>::operator+(const T v) const
{
    Matrix<T> out(this->rows, this->cols);
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            out.set(i, j, (this->get(i, j) + v));
        }
    }

    return out;
}

template <class T>
void Matrix<T>::operator+=(const T v)
{
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            this->set(i, j, (this->get(i, j) + v));
        }
    }
}

template <class T>
Matrix<T> Matrix<T>::operator-(const Matrix& v) const
{
    Matrix<T> out(this->rows, this->cols);
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            out.set(i, j, (this->get(i, j) - v.get(i, j)));
        }
    }

    return out;
}

template <class T>
void Matrix<T>::operator-=(const Matrix& v)
{
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            this->set(i, j, (this->get(i, j) - v.get(i, j)));
        }
    }
}

template <class T>
Matrix<T> Matrix<T>::operator-(const T v) const
{
    Matrix<T> out(this->rows, this->cols);
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            out.set(i, j, (this->get(i, j) - v));
        }
    }

    return out;
}

template <class T>
void Matrix<T>::operator-=(const T v)
{
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            this->set(i, j, (this->get(i, j) - v));
        }
    }
}

template <class T>
Matrix<T> Matrix<T>::operator*(const Matrix& v) const
{
    Matrix<T> out(this->rows, this->cols);
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            out.set(i, j, (this->get(i, j) * v.get(i, j)));
        }
    }

    return out;
}

template <class T>
void Matrix<T>::operator*=(const Matrix& v)
{
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            this->set(i, j, (this->get(i, j) * v.get(i, j)));
        }
    }
}

template <class T>
Matrix<T> Matrix<T>::operator*(const T v) const
{
    Matrix<T> out(this->rows, this->cols);
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            out.set(i, j, (this->get(i, j) * v));
        }
    }

    return out;
}

template <class T>
void Matrix<T>::operator*=(const T v)
{
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            this->set(i, j, (this->get(i, j) * v));
        }
    }
}

template <class T>
Matrix<T> Matrix<T>::operator/(const Matrix& v) const
{
    Matrix<T> out(this->rows, this->cols);
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            out.set(i, j, (this->get(i, j) / v.get(i, j)));
        }
    }

    return out;
}

template <class T>
void Matrix<T>::operator/=(const Matrix& v)
{
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            this->set(i, j, (this->get(i, j) / v.get(i, j)));
        }
    }
}

template <class T>
Matrix<T> Matrix<T>::operator/(const T v) const
{
    Matrix<T> out(this->rows, this->cols);
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            out.set(i, j, (this->get(i, j) / v));
        }
    }

    return out;
}

template <class T>
void Matrix<T>::operator/=(const T v)
{
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            this->set(i, j, (this->get(i, j) / v));
        }
    }
}

template <class T>
void Matrix<T>::operator=(const Matrix& M)
{
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            T v = M.get(i, j);
            this->set(i, j, v);
        }
    }
}

template <class T>
void Matrix<T>::operator=(const T v)
{
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            this->set(i, j, v);
        }
    }
}

template <class T>
T Matrix<T>::sum(void) const
{
    T out = 0;

    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            out += this->get(i, j);
        }
    }

    return out;
}

template <class T>
Matrix<T> Matrix<T>::AddVec(const Vector<T>& v, unsigned int idx)
{
    Matrix<T> out(this->rows, this->cols);
    out = this->copy();
    if (v.is_row_vec) {
        assert(v.length == this->get_row(idx).length);
            out.set_row(idx, out.get_row(idx) + v);
    } else {
        assert(v.length == this->get_col(idx).length);
        for (unsigned int i = 0; i < v.length; i++) {
            out.set_col(idx, out.get_col(idx) + v);
        }
    }

    return out;
}

template <class T>
void Matrix<T>::AddVecTo(const Vector<T>& v, unsigned int idx)
{
    if (v.is_row_vec) {
        assert(v.length == this->get_row(idx).length);
        this->set_row(idx, this->get_row(idx) + v);
    } else {
        assert(v.length == this->get_col(idx).length);
        this->set_col(idx, this->get_col(idx) + v);
    }
}

template <class T>
Matrix<T> Matrix<T>::copy(void) const
{
    Matrix<T> out(this->rows, this->cols);
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            out.set(i, j, this->get(i, j));
        }
    }

    return out;
}

template <class T>
T Matrix<T>::reduce(void)
{
    Vector<T> * row = new Vector<T> (this->cols, 1);
    unsigned int lead = 0;
    T d = 1;
    T lv;

    for (unsigned int r = 0; r < this->rows; r++) {
        if (lead >= this->cols) {
            return d;
        }
        unsigned int i = r;
        while (this->get(i, lead) == 0) {
            i += 1;
            if (i == this->rows) {
                i = r;
                lead += 1;
                if (this->cols == lead) {
                    return d;
                }
            }
        }
        *row = this->get_row(r);
        this->set_row(r, this->get_row(i));
        this->set_row(i, *row);
        d *= -1;
        lv = this->get(r, lead);
        *row = this->get_row(r);
        this->set_row(r, *row*(1/lv));
        d *= 1/lv;
        for (unsigned int j = 0; j < this->rows; j++) {
            if (j != r) {
                lv = this->get(j, lead);
                this->set_row(j, this->get_row(j)-(this->get_row(r)*lv));
            }
        }
        lead +=1;
    }

    delete row;

    return d;
}

template <class T>
T Matrix<T>::echelon(void)
{
    T h = 0;
    T k = 0;
    T d = 1;

    while ((h < this->rows) & (k < this->cols)) {
        unsigned int imax = this->argmax(0, 1, h, k);
        if (this->get(imax, k) == 0) {
            k += 1;
        } else {
            this->swap(h, imax);
            d *= -1;
            for (unsigned int i = h + 1; i < this->rows; i++) {
                T f = this->get(i, k)/this->get(h, k);
                this->set(i, k, 0);
                for (unsigned int j = k + 1; j < this->cols; j++) {
                    this->set(i, j, (this->get(i, j) - (this->get(h, j) * f)));
                }
            }

            h += 1;
            k += 1;
        }
    }

    return d;
}

template <class T>
void Matrix<T>::invert(void)
{
    assert(is_square);

    Matrix<T> Id(this->rows, this->cols);
    Id.I();
    Matrix<T> c(this->rows, this->cols);
    c = this->copy();
    Matrix<T> a(this->rows, this->cols*2);
    a = c.augment(Id);
    a.reduce();
    Matrix<T> test(a.rows, a.rows);
    test = a.slice(0,a.rows-1,0,a.rows-1);
    assert(test.compare(Id));
    if (test.compare(Id)) {
        Matrix<T> inv(a.rows, a.rows);
        inv = a.slice(0,a.rows-1,a.rows,2*a.rows-1);
        for (unsigned int i = 0; i < this->rows; i++) {
            this->set_row(i, inv.get_row(i));
        }
    } else {
        for (unsigned int i = 0; i < this->rows; i++) {
            this->set_row(i, Id.get_row(i));
        }
    }
}

template <class T>
Matrix<T> Matrix<T>::augment(const Matrix& M)const
{
    assert(this->rows == M.rows);
    Matrix<T> out(this->rows, this->cols+M.cols);
    for (unsigned int i = 0; i < this->cols; i++) {
        out.set_col(i,this->get_col(i));
    }
    for (unsigned int i = this->cols; i < out.cols; i++) {
        out.set_col(i, M.get_col(i-this->cols));
    }

    return out;
}

template <class T>
Matrix<T> Matrix<T>::slice(unsigned int r1, unsigned int r2, unsigned int c1, unsigned int c2) const
{
    Matrix<T> out(r2-r1+1,c2-c1+1);
    for (unsigned int i = 0; i < out.rows; i++) {
        for (unsigned int j = 0; j < out.cols; j++) {
            out.set(i, j, this->get(r1+i,c1+j));
        }
    }

    return out;
}

template <class T>
T Matrix<T>::trace(void) const
{
    T out = 0;

    for (unsigned int i = 0; i < this->rows; i++) {
        out += this->get(i, i);
    }

    return out;
}

template <class T>
void Matrix<T>::transpose(void)
{
    Matrix<T> t(this->cols, this->rows);
    t = this->copy_transpose();
    unsigned int temp = this->rows;
    this->rows = this->cols;
    this->cols = temp;
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            this->set(i, j, t.get(i, j));
        }
    }
}

template <class T>
Matrix<T> Matrix<T>::copy_transpose(void) const
{
    Matrix<T> out(this->cols, this->rows);
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            out.set(j, i, this->get(i, j));
        }
    }

    return out;
}

template <class T>
void Matrix<T>::conjugate_transpose(void)
{
    Matrix<T> t(this->cols, this->rows);
    t = this->copy_conjugate_transpose();
    unsigned int temp = this->rows;
    this->cols = temp;
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            this->set(i, j, t.get(i, j));
        }
    }
}

template <class T>
Matrix<T> Matrix<T>::copy_conjugate_transpose(void) const
{
    Matrix<std::complex<T> > out(this->cols, this->rows);
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            out.set(j, i, std::conj(this->get(i, j)));
        }
    }

    return out;
}

template <class T>
T Matrix<T>::inner(const Matrix& M) const
{
    T out = 0;
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            out += this->get(i, j)*M.get(i, j);
        }
    }

    return out;
}

template <class T>
Matrix<T> Matrix<T>::outer(const Matrix& M) const
{
    assert(this->cols == M.rows);
    Matrix<T> out(this->rows, M.cols);
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int k = 0; k < this->cols; k++) {
            for (unsigned int j = 0; j < M.cols; j++) {
                out(i, j) += this->get(i, k)*M.get(k, j);
            }
        }
    }

    return out;
}

template <class T>
Matrix<T> Matrix<T>::submatrix(unsigned int row, unsigned int col) const
{
    Matrix<T> R(this->rows-1, this->cols);
    Matrix<T> out(this->rows-1, this->cols-1);
    R = this->del_row(row);
    out = R.del_col(col);
    return out;
}

template <class T>
T Matrix<T>::minor(unsigned int row, unsigned int col) const
{
    Matrix<T> sub(this->rows-1, this->cols-1);
    sub = this->submatrix(row, col);
    T m = sub.echelon_determinant();
    return m;
}

template <class T>
T Matrix<T>::cofactor(unsigned int row, unsigned int col) const
{
    T m = this->minor(row, col);
    T c = m*pow(-1, row+col);
    return c;
}

template <class T>
T Matrix<T>::echelon_determinant(void) const
{
    Matrix<T> U(this->rows, this->cols);
    Matrix<T> B(this->rows, this->cols);
    Matrix<T> L(this->rows, this->cols);
    U = this->copy();
    B = this->copy();
    B.invert();
    L = this->copy();
    T d = U.reduce();
    T D = 1/d;
    for (unsigned int i = 0; i < B.rows; i++) {
        D *= B.get(i, i);
    }
    d = L.echelon();
    D = D*d;
    for (unsigned int i = 0; i < B.rows; i++) {
        assert(B.get(i, i) != 0);
        D *= 1/(B.get(i, i));
    }

    return D;
}

template <class T>
T Matrix<T>::laplace_determinant(void) const
{
    T D = 0;
    for (unsigned int j = 0; j < this->cols; j++) {
        D += this->get(0, j)*this->cofactor(0, j);
    }

    return D;
}

template <class T>
Vector<std::complex<T> > Matrix<T>::eigenvalues(unsigned int iterations) const
{
    Matrix<T> A(this->rows, this->cols);
    A = this->copy();
    Matrix<T> Q(this->rows, this->cols);
    Matrix<T> R(this->rows, this->cols);

    A.QR(&Q, &R);
    Matrix<T> q(this->rows, this->cols);
    q = Q.copy();
    q.transpose();
    A = (Q.outer(*this)).outer(q);
    for (unsigned int i = 0; i < iterations; i++) {
        A.QR(&Q, &R);
        A = R.outer(Q);
    }
    Vector<std::complex<T> > E(this->rows, 1);
    bool flag = 0;
    for (unsigned int j = 0; j < this->rows; j++) {
        if ((j <= this->rows-2) | flag) {
            if (flag) {
                flag = 0;
            } else if (abs(A.get(j+1, j)) > pow(10,-5)) {
                std::complex<T> z(((A.get(j, j)+A.get(j+1, j+1))/2),std::sqrt(abs(A.get(j, j+1))*abs(A.get(j+1, j))));
                E.set(j, z);
                E.set(j+1, std::conj(E.get(j)));
                flag = 1;
            } else {
                E.set(j, A.get(j, j));
            }
        } else {
            E.set(j, A.get(j, j));
        }
    }

    return E;
}

template <class T>
T Matrix<T>::condition(unsigned int n) const
{
    Matrix<T> A(this->rows, this->cols);
    Matrix<T> AI(this->rows, this->cols);
    A = this->copy();
    AI = this->copy();
    AI.invert();

    T ret = 0;

    switch (n) {
        case 0:
            ret = AI.infnorm()*A.infnorm();
            break;
        case 1:
            ret = AI.onenorm()*A.onenorm();
            break;
        case 2:
            ret = AI.twonorm()*A.twonorm();
            std::cout << A.twonorm() << std::endl;
            break;
    }
    return ret;
}

template <class T>
T Matrix<T>::onenorm(void) const
{
    T max = 0;
    T sum = 0;
    for (unsigned int j = 0; j < this->cols; j++) {
        for (unsigned int i = 0; i < this->rows; i++) {
            sum += abs(this->get(i, j));
        }
        if (sum > max)
            max = sum;
        sum = 0;
    }
    return max;
}

template <class T>
T Matrix<T>::twonorm(void) const
{
    Matrix<T> A(this->rows, this->cols);
    Matrix<T> AC(this->cols, this->rows);
    Matrix<T> P(this->cols, this->cols);
    A = this->copy();
    AC = this->copy_transpose();
    P = AC.outer(A);
    Vector<std::complex<T> > E(this->cols, 1);
    E = P.eigenvalues(10000);
    std::complex<T> m = E.max();
    return sqrt(abs(m));
}

template <class T>
T Matrix<T>::infnorm(void) const
{
    T max = 0;
    T sum = 0;
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            sum += abs(this->get(i, j));
        }
        if (sum > max)
            max = sum;
        sum = 0;
    }
    return max;
}

template <class T>
void Matrix<T>::QR(Matrix *Q, Matrix *R) const
{
    T m = this->rows;
    T n = this->cols;
    Matrix<T> A(m, n);
    A = this->copy();
    Vector<T> d(n, 0);
    for (unsigned int j = 0; j < n; j++) {
        T s = 0;
        for (unsigned int i = j; i < m; i++) {
            s += pow(A.get(i, j), 2);
        }
        s = sqrt(s);
        if (A.get(j, j) > 0) {
            d.set(j, -s);
        } else {
            d.set(j, s);
        }
        T fak = sqrt(s*(s+abs(A.get(j, j))));
        A.set(j, j, A.get(j, j) - d.get(j));
        for (unsigned int k = j; k < m; k++) {
            A.set(k, j, A.get(k, j)/fak);
        }
        for (unsigned int i = j + 1; i < m; i++) {
            s = 0;
            for (unsigned int k = j; k < m; k++) {
                s += A.get(k, j)*A.get(k, i);
            }
            for (unsigned int k = j; k < m; k++) {
                A.set(k, i, A.get(k, i)-(A.get(k, j)*s));
            }
        }
    }

    for (unsigned int i = 0; i < A.rows; i++) {
        for (unsigned int j = 0; j < A.cols; j++) {
            if (abs(A.get(i, j)) < pow(10,-7.5))
                A.set(i, j, 0.0f);
        }
    }

    // Populate R
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            if (j > i) {
                R->set(i, j, A.get(i, j));
            } else if (i == j) {
                R->set(i, j, d.get(j));
            }
        }
    }

    // Populate Q

    Q->I();
    Matrix<T> I(this->rows, this->cols);
    I.I();
    Vector<T> * w = new Vector<T>(this->cols, 0);
    Vector<T> * W = new Vector<T>(this->cols, 1);
    Matrix<T> * P = new Matrix<T>(this->cols, this->cols);
    for (unsigned int i = 0; i < this->cols; i++) {
        w->clean();
        for (unsigned int j = 0; j < this->cols; j++) {
            if (j >= i) {
                w->set(j, A.get(j, i));
            }
        }
        *W = w->copy();
        *P = I - W->outer(*w);
        *Q = Q->outer(*P);
    }

    delete w;
    delete W;
    delete P;
}

template <class T>
void Matrix<T>::QR_fast(Matrix *Q, Matrix *R) const
{
    unsigned int s = this->rows;
    unsigned int n = this->rows;
    unsigned int z = 0;
    Matrix<T> A(this->rows, this->cols);
    A = this->copy();
    Matrix<T> M(this->rows, this->cols);
    M = this->copy();
    Vector<T> v(this->rows, 0);
    Vector<T> u(this->rows, 0);

    for (unsigned int c = 0; c < n; c++) {
        if (c == (z + s)) {
            Matrix<T> BQ(n, (c - 1) - z);
            Matrix<T> BM(n, n - c);
            BQ = Q->slice(0, n, z, c - 1);
            BM = M.slice(0, n, c, n);
            Matrix<T> C((c - 1) - z, n - c);
            C = BQ.copy_transpose().outer(BM);
            Matrix<T> S(n, n - c);
            S = BQ.outer(C);
            for (unsigned int i = 0; i < n; i++) {
                for (unsigned int j = c; j < n; j++) {
                    M(i, j) -= S(i, j - c + 1);
                }
            }
            z = c;
        }
        v = M.get_col(c);
        for (unsigned int j = z; j < c; j++) {
            u = Q->get_col(j);
            v -= (u * ((u.copy_transpose()).outer(v))(0, 0));
        }
        v.unit();
        for (unsigned int i = 0; i < n; i++) {
            (*Q)(i, c) = v(i);
        }
    }
    (*R) = (Q->copy_transpose()).outer((*this));
}

template <class T>
Matrix<T> Matrix<T>::Cholesky_fast(void) const
{
    unsigned int s = this->rows;
    unsigned int n = this->rows;
    unsigned int z = 0;
    Matrix<T> L(this->rows, this->cols);
    Matrix<T> A(this->rows, this->cols);
    A = this->copy();

    for (unsigned int c = 0; c < n; c++) {
        /*if (c == (z + s)) {
            Matrix<T> R(n - c, (c - 1) - z);
            R = L.slice(c, n, z, c - 1);
            Matrix<T> S(n - c, n - c);
            S = R.outer(R.copy_transpose());
            for (unsigned int i = c; i < n; i++) {
                for (unsigned int j = c; j < n; j++) {
                    A(i, j) -= S(i - c, j - c + 1);
                }
            }
            z = c;
        }*/
        L(c, c) = A(c, c);
        for(unsigned int k = z; k < c; k++) {
            L(c, c) -= (L(c, k) * L(c, k));
        }
        L(c, c) = sqrt(L(c, c));
        for (unsigned int i = (c + 1); i < n; i++) {
            L(i, c) = A(i, c);
            for (unsigned int k = z; k < c; k++) {
                L(i, c) -= (L(i, k)*L(c, k));
            }
            L(i, c) /= L(c, c);
        }
    }

    return L;
}

template <class T>
void Matrix<T>::LU_fast(Matrix *L, Matrix *U) const
{
    unsigned int s = this->rows;
    unsigned int n = this->rows;
    unsigned int z = 0;
    Matrix<T> A(this->rows, this->cols);
    A = this->copy();

    for (unsigned int c = 0; c < n; c++) {
        if (c == (z + s)) {
            Matrix<T> RL(n - c, (c - 1) - z);
            Matrix<T> RU((c - 1) - z, n - c);
            RL = L->slice(c, n, z, c - 1);
            RU = U->slice(z, c - 1, c, n);
            Matrix<T> S(n - c, n - c);
            S = RL.outer(RU);
            for (unsigned int i = c; i < n; i++) {
                for (unsigned int j = c; j < n; j++) {
                    A(i, j) -= S(i - c + 1, j - c + 1);
                }
            }
        }
        for (unsigned int i = c; i < n; i++) {
            (*L)(i, c) = A(i, c);
            for (unsigned int k = z; k < c; k++) {
                (*L)(i, c) -= ((*L)(i, k)*(*U)(k, c));
            }
        }
        for (unsigned int i = c; i < n; i++) {
            (*U)(c, i) = A(c, i);
            for (unsigned int k = z; k < c; k++) {
                (*U)(c, i) -= ((*L)(c, k)*(*U)(k, i));
            }
            (*U)(c, i) /= (*L)(c, c);
        }
    }
}

template <class T>
void Matrix<T>::SVD(Matrix *E, Matrix *U, Matrix *V) const
{

}

template <class T>
T& Matrix<T>::operator()(unsigned int row, unsigned int col)
{
    return this->data[row*this->cols + col];
}

template <class T>
T Matrix<T>::operator()(unsigned int row, unsigned int col) const
{
    return this->data[row*this->cols + col];
}

template <class T>
T Matrix<T>::get(unsigned int row, unsigned int col) const
{
    return this->data[row*this->cols + col];
}

template <class T>
void Matrix<T>::set(unsigned int row, unsigned int col, T val)
{
    this->data[row*this->cols + col] = val;
}

template <class T>
Vector<T> Matrix<T>::get_row(unsigned int r) const
{
    Vector<T> out(this->cols, 1);
    for (unsigned int j = 0; j < this->cols; j++) {
        out.set(j, this->get(r, j));
    }

    return out;
}

template <class T>
Vector<T> Matrix<T>::get_col(unsigned int c) const
{
    Vector<T> out(this->rows, 0);
    for (unsigned int i = 0; i < this->rows; i++) {
        out.set(i, this->get(i, c));
    }

    return out;
}

template <class T>
void Matrix<T>::set_row(unsigned int r, const Vector<T>& row)
{
    for (unsigned int j = 0; j < this->cols; j++) {
        this->set(r, j, row.get(j));
    }
}

template <class T>
void Matrix<T>::set_col(unsigned int c, const Vector<T>& col)
{
    for (unsigned int i = 0; i < this->rows; i++) {
        this->set(i, c, col.get(i));
    }
}

template <class T>
void Matrix<T>::swap(unsigned int row1, unsigned int row2)
{
    Vector<T> r(this->cols, 1);
    r = this->get_row(row1);
    this->set_row(row1, this->get_row(row2));
    this->set_row(row2, r);
}

template <class T>
Matrix<T> Matrix<T>::del_row(unsigned int row) const
{
    Matrix<T> out(this->rows-1, this->cols);
    bool flag = 0;
    for (unsigned int i = 0; i < out.rows; i++) {
        for (unsigned int j = 0; j < out.cols; j++) {
            if ((i == (row)) | flag) {
                out.set(i, j, this->get(i+1, j));
                flag = 1;
            } else {
                out.set(i, j, this->get(i, j));
            }
        }
    }

    return out;
}

template <class T>
Matrix<T> Matrix<T>::del_col(unsigned int col) const
{
    Matrix<T> out(this->rows, this->cols-1);
    bool flag = 0;
    for (unsigned int j = 0; j < out.cols; j++) {
        for (unsigned int i = 0; i < out.rows; i++) {
            if ((j == (col)) | flag) {
                out.set(i, j, this->get(i, j+1));
                flag = 1;
            } else {
                out.set(i, j, this->get(i, j));
            }
        }
    }

    return out;
}

template <class T>
unsigned int Matrix<T>::argmax(unsigned int row, unsigned int col, unsigned int i0, unsigned int i1) const
{
    T max = 0;
    unsigned int ind = 0;

    if ((row == 0) & (col != 0)) {
        for (unsigned int i = i0; i < this->rows; i++) {
            if (abs(this->get(i, i1)) > abs(max)) {
                max = this->get(i, i1);
                ind = i;
            }
        }
    } else if ((row != 0) & (col == 0)) {
        for (unsigned int j = i1; j < this->cols; j++) {
            if (abs(this->get(i0, j)) > abs(max)) {
                max = this->get(i0, j);
                ind = j;
            }
        }
    }

    return ind;
}

template <class T>
bool Matrix<T>::compare(const Matrix<T>& M) const
{
    assert((this->rows == M.rows) & (this->cols == M.cols));
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            if (this->get(i, j) - M.get(i, j) > pow(10,-10)) {
                return 0;
            }
        }
    }

    return 1;
}

template <class T>
void Matrix<T>::random(T min, T max)
{
    srand(time(NULL));
    for (unsigned int i = 0; i < this->rows; i++) {
        for (unsigned int j = 0; j < this->cols; j++) {
            this->set(i, j, ((T)rand()/(T)(RAND_MAX)) * (max - min) + min);
        }
    }
}

template <class T>
void Matrix<T>::print(unsigned int w) const
{
    for (unsigned int i = 0; i < this->rows; i++) {
        std::cout << "|";
        for (unsigned int j = 0; j < this->cols; j++) {
            std::cout << std::setw(w) << std::setprecision(w) << std::fixed << std::setfill(' ') << this->get(i, j);
            if (j != cols - 1) {
                std::cout << ",";
            }
        }
        std::cout << "|" << std::endl;
    }
}

template <class T>
void Matrix<T>::print(unsigned int w, unsigned int p) const
{
    for (unsigned int i = 0; i < this->rows; i++) {
        std::cout << "|";
        for (unsigned int j = 0; j < this->cols; j++) {
            std::cout << std::setw(w) << std::setprecision(p) << std::fixed << std::setfill(' ') << this->get(i, j);
            if (j != cols - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "|" << std::endl;
    }
}

template <class T>
void Matrix<T>::verify_square(void)
{
    this->square = ((this->rows == this->cols) ? 1 : 0);
}

template <class T>
void Matrix<T>::verify_symmetric(void)
{
    this->is_symmetric = 1;
    if (this->is_square == 0) {
        this->is_symmetric = 0;
        return;
    }
    this->is_symmetric = this->compare(this->copy_transpose());
}

template <class T>
void Matrix<T>::verify_definite(void)
{

}

template <class T>
void Matrix<T>::verify_hermitian(void)
{
    this->is_hermitian = 1;
    if (this->is_square == 0) {
        this->is_hermitian = 0;
        return;
    }
    this->is_hermitian = this->compare(this->copy_conjugate_transpose());
}
