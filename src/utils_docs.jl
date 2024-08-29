@doc raw"""
    separation_matrix(ϕ::NTuple{3,<:Real}, μ::NTuple{3,<:Real}=(1, 1, 1))
    separation_matrix(ϕ_0::Real,...) # Assume equidistant phases

Construct a separation matrix `M` with `eltype(M) = Complex{T}` for given phase shifts `ϕ` and modulations `μ`.

The matrix `M` is such that for acquisitions ``D_1(\tilde{k})``, ``D_2(\tilde{k})`` and ``D_3(\tilde{k})`` the
components are separated as
```math
    \begin{pmatrix}
    C_{0}(\tilde{k}) \\
    C_{-1}(\tilde{k}) \\
    C_{+1}(\tilde{k})
    \end{pmatrix} =
    \begin{pmatrix}
    H(\tilde{k})S(\tilde{k}) \\
    H(\tilde{k})S(\tilde{k} - \tilde{\nu}) \\
    H(\tilde{k})S(\tilde{k} + \tilde{\nu})
    \end{pmatrix} =
    \begin{pmatrix}
    1 & \tfrac{\mu_1}{2} e^{-i\phi_1} & \tfrac{\mu_1}{2} e^{i \phi_1} \\
    1 & \tfrac{\mu_2}{2} e^{-i\phi_2} & \tfrac{\mu_2}{2} e^{i \phi_2} \\
    1 & \tfrac{\mu_3}{2} e^{-i\phi_3} & \tfrac{\mu_3}{2} e^{i \phi_3}
    \end{pmatrix}^{-1}
    \begin{pmatrix}
    D_{1}(\tilde{k}) \\
    D_{2}(\tilde{k}) \\
    D_{3}(\tilde{k})
    \end{pmatrix}
```
"""
separation_matrix
