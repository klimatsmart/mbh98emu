#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Implementation of Businger and Golub's singular value decomposition
algorithm for real matrices.
"""

import numpy as np

ETA = 1.2e-7
TOL = 2.4e-32


def svd(a):
    """
    Singular value decomposition.
    
    Parameters
    ----------
    a : array
        A real 2D array.
    
    Returns
    -------
    u : array
        Left singular vectors.
    s : array
        Singular values.
    v : array
        Right singular vectors.
    """
    a = a.copy()
    if a.shape[0] < a.shape[1]:
        transpose = True
        a = a.T
    else:
        transpose = False
    
    a, b, c = _bidiagonal_reduction(a)
    u, s, v = _diagonal_reduction(a, b, c)
    u, s, v = _sorted_reduction(u, s, v)
    u = _u_back_transformation(u, a, b)
    v = _v_back_transformation(v, a, c)
    
    if transpose:
        u, v = v, u
    return u, s, v


def _bidiagonal_reduction(a):
    # Reduce a to bidiagonal form by Householder reflections.
    m, n = a.shape
    b = np.empty(n)
    c = np.empty(n)
    c[0] = 0
    for k in range(n):
        # Eliminate column k of a.
        z = np.float32(0)
        for i in range(k, m):
            z += np.square(a[i, k], dtype=np.float32)
        b[k] = 0
        if z > TOL:
            z = np.sqrt(z)
            b[k] = z
            w = np.float32(np.abs(a[k, k]))
            q = 1
            if w != 0:
                q = a[k, k] / w
            a[k, k] = q * (z+w)
            if k < n-1:
                p = np.ascontiguousarray(a[k:, k, np.newaxis] * a[k:, k+1:])
                if k < n-2:
                    q = p.sum(0)
                else:
                    q = 0
                    for x in p:
                        q += x
                q /= (z*(z+w))
                a[k:, k+1:] -= q * a[k:, k, np.newaxis]
                # Phase transformation.
                q = -np.sign(a[k, k])
                a[k, k+1:] *= q
        # Eliminate row k of a.
        if k < n-1:
            z = np.float32(0)
            for j in range(k+1, n):
                z += np.square(a[k, j], dtype=np.float32)
            c[k+1] = 0
            if z > TOL:
                z = np.sqrt(z)
                c[k+1] = z
                w = np.float32(np.abs(a[k, k+1]))
                q = 1
                if w != 0:
                    q = a[k, k+1] / w
                a[k, k+1] = q * (z+w)
                p = np.asfortranarray(a[k, k+1:] * a[k+1:, k+1:])
                q = p.sum(1)
                q /= (z*(z+w))
                a[k+1:, k+1:] -= q[:, np.newaxis] * a[k, k+1:]
                # Phase transformation.
                q = -np.sign(a[k, k+1])
                a[k+1:, k+1] *= q
    return a, b, c


def _diagonal_reduction(a, b, c):
    # Initialize u and v.
    m, n = a.shape
    u = np.eye(m, n)
    v = np.eye(n, n)
    
    # Copy diagonals.
    s = b.copy()
    t = c.copy()
    
    # Tolerance for negligible elements.
    eps = np.max(s + t)
    eps *= ETA
    
    # QR diagonalization of bidiagonal matrix.
    for k in reversed(range(n)):
        convergence = False
        while not convergence:
            for l in reversed(range(k+1)):
                if np.abs(t[l]) > eps and np.abs(s[l-1]) > eps:
                    continue
                if np.abs(t[l]) > eps and np.abs(s[l-1]) <= eps:
                    # Cancellation by Givens rotations.
                    cs = np.float32(0)
                    sn = np.float32(1)
                    for i in range(l, k+1):
                        f = np.float32(sn * t[i])
                        t[i] *= cs
                        if np.abs(f) <= eps:
                            # Go to test for convergence.
                            break
                        h = np.float32(s[i])
                        w = np.sqrt(f*f + h*h, dtype=np.float32)
                        s[i] = w
                        cs = np.float32(h / w)
                        sn = -np.float32(f / w)
                        x = u[:n, l-1].astype(np.float32)
                        y = u[:n, i].astype(np.float32)
                        u[:n, l-1] = x*cs + y*sn
                        u[:n, i] = y*cs - x*sn
                
                # Test for convergence.
                w = np.float32(s[k])
                if l == k:
                    convergence = True
                    break
                
                # Origin shift.
                x = np.float32(s[l])
                y = np.float32(s[k-1])
                g = np.float32(t[k-1])
                h = np.float32(t[k])
                f = np.float32(((y-w)*(y+w) + (g-h)*(g+h)) / (2*h*y))
                g = np.sqrt(f*f + 1, dtype=np.float32)
                if f < 0:
                    g = -g
                f = np.float32(((x-w)*(x+w) + (y/(f+g)-h)*h) / x)
                
                # QR step by Givens rotations.
                cs = np.float32(1)
                sn = np.float32(1)
                for i in range(l+1, k+1):
                    g = np.float32(t[i])
                    y = np.float32(s[i])
                    h = np.float32(sn * g)
                    g = np.float32(cs * g)
                    w = np.sqrt(h*h + f*f, dtype=np.float32)
                    t[i-1] = w
                    cs = np.float32(f / w)
                    sn = np.float32(h / w)
                    f = np.float32(x*cs + g*sn)
                    g = np.float32(g*cs - x*sn)
                    h = np.float32(y * sn)
                    y = np.float32(y * cs)
                    x = v[:, i-1].astype(np.float32)
                    w = v[:, i].astype(np.float32)
                    v[:, i-1] = x*cs + w*sn
                    v[:, i] = w*cs - x*sn
                    w = np.sqrt(h*h + f*f, dtype=np.float32)
                    s[i-1] = w
                    cs = np.float32(f / w)
                    sn = np.float32(h / w)
                    f = np.float32(cs*g + sn*y)
                    x = np.float32(cs*y - sn*g)
                    y = u[:n, i-1].astype(np.float32)
                    w = u[:n, i].astype(np.float32)
                    u[:n, i-1] = y*cs + w*sn
                    u[:n, i] = w*cs - y*sn
                t[l] = 0
                t[k] = f
                s[k] = x
                
                # Restart loop over l.
                break
        
        # Convergence.
        if w < 0:
            s[k] = -w
            v[:, k] = -v[:, k]
    
    return u, s, v


def _sorted_reduction(u, s, v):
    index = np.flip(np.argsort(s))
    s = s[index]
    u = u[:, index]
    v = v[:, index]
    return u, s, v


def _u_back_transformation(u, a, b):
    n = a.shape[1]
    for k in reversed(range(n)):
        if b[k] != 0:
            q = -np.sign(a[k, k])
            u[k, :] *= q
            p = np.ascontiguousarray(a[k:, k, np.newaxis] * u[k:, :])
            if n == 1:
                q = 0
                for x in p:
                    q += x
            else:
                q = p.sum(0)
            q /= (np.abs(a[k, k]) * b[k])
            u[k:, :] -= q * a[k:, k, np.newaxis]
    return u


def _v_back_transformation(v, a, c):
    n = a.shape[1]
    for k in reversed(range(n-1)):
        if c[k+1] != 0:
            q = -np.sign(a[k, k+1])
            v[k+1, :] *= q
            p = np.ascontiguousarray(a[k, k+1:, np.newaxis] * v[k+1:, :])
            q = p.sum(0)
            q /= (np.abs(a[k, k+1]) * c[k+1])
            v[k+1:, :] -= q * a[k, k+1:, np.newaxis]
    return v
