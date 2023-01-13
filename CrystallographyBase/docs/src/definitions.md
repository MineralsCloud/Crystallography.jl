# Definitions and conventions

```@contents
Pages = ["definitions.md"]
Depth = 3
```

## Basis vectors

In `CrystallographyBase`, basis vectors are represented by three-column vectors:

```math
\mathbf{a} = \begin{bmatrix}
    a_x \\
    a_y \\
    a_z
\end{bmatrix},
\quad
\mathbf{b} = \begin{bmatrix}
    b_x \\
    b_y \\
    b_z
\end{bmatrix},
\quad
\mathbf{c} = \begin{bmatrix}
    c_x \\
    c_y \\
    c_z
\end{bmatrix},
```

in Cartesian coordinates. Depending on the situation,
``\begin{bmatrix} \mathbf{a}_1 & \mathbf{a}_2 & \mathbf{a}_3 \end{bmatrix}``
is used instead of
``\begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix}``.

Therefore, a lattice is represented as

```math
\mathrm{A} =
\begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix} =
\begin{bmatrix}
    a_x & b_x & c_x \\
    a_y & b_y & c_y \\
    a_z & b_z & c_z
\end{bmatrix}.
```

A reciprocal lattice is its inverse, represented as three row vectors:

```math
\mathrm{B} =
\mathrm{A}^{-1} =
\begin{bmatrix}
    \mathbf{b}_1 \\
    \mathbf{b}_2 \\
    \mathbf{b}_3
\end{bmatrix},
```

so that

```math
\mathrm{A} \mathrm{B} = \mathrm{B} \mathrm{A} = \mathrm{I},
```

where ``\mathrm{I}`` is the ``3 \times 3`` identity matrix.

## Crystal coordinates

Coordinates of an atomic point ``\mathbf{x}`` are represented
as three fractional values relative to basis vectors as follows,

```math
\mathbf{x} = \begin{bmatrix}
    x_1 \\
    x_2 \\
    x_3
\end{bmatrix},
```

where ``0 \le x_i < 1``. A position vector ``\mathbf{r}`` in
Cartesian coordinates is obtained by

```math
\mathbf{r} = \mathrm{A} \mathbf{x} = \begin{bmatrix} \mathbf{a}_1 & \mathbf{a}_2 & \mathbf{a}_3 \end{bmatrix} \mathbf{x},
```

or

```math
\mathbf{r} = \sum_i x_i \mathbf{a}_i.
```

## Transformation to the primitive cell

In the standardized unit cells, there are five different centring
types available, base centrings of A and C, rhombohedral (R), body-centred (I),
and face-centred (F). The transformation is applied to the
standardized unit cell by

```math
\begin{bmatrix} \mathbf{a}_p & \mathbf{b}_p & \mathbf{c}_p \end{bmatrix} =
\begin{bmatrix} \mathbf{a}_s & \mathbf{b}_s & \mathbf{c}_s \end{bmatrix}
\mathrm{P}
```

where ``\mathbf{a}_p``, ``\mathbf{b}_p``, and ``\mathbf{c}_p``
are the basis vectors of the primitive cell and ``\mathrm{P}`` is the
transformation matrix from the standardized unit cell to the primitive
cell. Matrices ``\mathrm{P}`` for different centring types are given as follows:

```math
\mathrm{P}_\text{A} = \begin{bmatrix}
    1 & 0 & 0 \\
    0 & \dfrac{1}{2} & \dfrac{-1}{2} \\
    0 & \dfrac{1}{2} & \dfrac{1}{2}
\end{bmatrix},
\quad
\mathrm{P}_\text{C} = \begin{bmatrix}
    \dfrac{1}{2} & \dfrac{1}{2} & 0 \\
    \dfrac{-1}{2} & \dfrac{1}{2} & 0 \\
    0 & 0 & 1
\end{bmatrix},
\quad
\mathrm{P}_\text{R} = \begin{bmatrix}
    \dfrac{2}{3} & \dfrac{-1}{3} & \dfrac{-1}{3} \\
    \dfrac{1}{3} & \dfrac{1}{3} & \dfrac{\bar{2}}{3} \\
    \dfrac{1}{3} & \dfrac{1}{3} & \dfrac{1}{3}
\end{bmatrix},
\quad
\mathrm{P}_\text{I} = \begin{bmatrix}
    \dfrac{-1}{2} & \dfrac{1}{2} & \dfrac{1}{2} \\
    \dfrac{1}{2} & \dfrac{-1}{2} & \dfrac{1}{2} \\
    \dfrac{1}{2} & \dfrac{1}{2} & \dfrac{-1}{2}
\end{bmatrix},
\quad
\mathrm{P}_\text{F} = \begin{bmatrix}
    0 & \dfrac{1}{2} & \dfrac{1}{2} \\
    \dfrac{1}{2} & 0 & \dfrac{1}{2} \\
    \dfrac{1}{2} & \dfrac{1}{2} & 0
\end{bmatrix}.
```

The choice of transformation matrix depends on the purpose.

For rhombohedral lattice systems with the H setting (hexagonal lattice),
``\mathrm{P}_\text{R}`` is applied to obtain
primitive basis vectors. However, with the R setting (rhombohedral lattice),
no transformation matrix is used because it is already a primitive cell.
