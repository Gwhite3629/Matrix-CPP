#include "solver.h"
#include "vector.h"
#include "matrix.h"

template <class T>
Vector<T> RK4(Vector<T> initial, T v_t0, T tf, unsigned int res, T (*f)(T, T))
{
    T dt = tf/res;
    Vector<T> out(initial.length);
    out = initial.copy();
    for (unsigned int i = 0; i < res; i++) {
        for (unsigned int j = 0; j < initial.length; j++) {
            T k1 = f(i*dt,out.get(j));
            T k2 = f(i*dt+dt/2,out.get(j)+dt*k1/2);
            T k3 = f(i*dt+dt/2,out.get(j)+dt*k2/2);
            T k4 = f(i*dt+dt,out.get(j)+dt*k3);
            out.set(j, out.get(j) + (1/6)*(k1+2*k2+2*k3+k4)*dt);
        }
    }

    return out;
}