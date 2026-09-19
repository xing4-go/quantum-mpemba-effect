# QuTiP 数值研究 Workflow：从 Long-Range QMPE 到 QMPE Reopening

## 0. 项目总目标

本项目研究 Quantum Mpemba Effect（QMPE）对以下因素的依赖：

1. 相互作用空间衰减形式；
2. 相互作用范围；
3. disorder / localization；
4. engineered local potential。

整个研究按照四个 Stage 展开：

| Stage | Hamiltonian | Scan | 核心问题 |
|---|---|---|---|
| I | Power-law XY | \(\alpha,\theta,N\) | 复现论文并确定 QMPE 存在区域 |
| II | Exponential XY | \(\xi,\theta,N\) | 缩短 interaction range 能否关闭 QMPE |
| III | XY + disorder | \(W,\theta,N\) | localization 如何抑制 QMPE |
| IV | Exponential XY + local potential | \(\xi,h_i,\theta,N\) | engineered potential 能否重新打开 QMPE |

研究逻辑为

\[
\boxed{
\text{QMPE baseline}
\rightarrow
\text{interaction engineering}
\rightarrow
\text{QMPE suppression}
\rightarrow
\text{QMPE reopening?}
}
\]

其中 Stage I 和 Stage III 有原论文作为 benchmark；Stage II 和 Stage IV 是扩展研究。

特别注意：

> Stage IV 的 “reopen QMPE” 是需要数值验证的研究假设，而不是预设结论。

---

# 1. 统一理论框架

四个 Stage 尽量使用同一套数值框架，只修改 Hamiltonian。

## 1.1 总 Hamiltonian

统一写成

\[
H
=
\sum_{i<j}
\frac{J_{ij}}{2}
\left(
\sigma_i^x\sigma_j^x+
\sigma_i^y\sigma_j^y
\right)
+
\sum_i h_i\sigma_i^z.
\]

利用

\[
\sigma_i^x\sigma_j^x+\sigma_i^y\sigma_j^y
=
2
\left(
\sigma_i^+\sigma_j^-
+
\sigma_i^-\sigma_j^+
\right),
\]

也可以写成

\[
H
=
\sum_{i<j}
J_{ij}
\left(
\sigma_i^+\sigma_j^-
+
\sigma_i^-\sigma_j^+
\right)
+
\sum_i h_i\sigma_i^z.
\]

这里：

- \(J_{ij}\)：控制 interaction spatial profile；
- \(\alpha\)：power-law exponent；
- \(\xi\)：exponential interaction length；
- \(h_i\)：local site energy / disorder potential；
- \(W\)：disorder strength。

---

# 2. 必须始终保持的 \(U(1)\) symmetry

定义 conserved charge

\[
Q
=
\frac12\sum_{i=1}^{N}\sigma_i^z.
\]

所有主要模型都应满足

\[
\boxed{
[H,Q]=0.
}
\]

这是整个研究的基础。

因为我们研究的是：

\[
\text{initially broken }U(1)
\rightarrow
\text{local }U(1)\text{ restoration}.
\]

因此每次构造新 Hamiltonian 后，程序第一件事都应该检查

```python
comm = H * Q - Q * H
print(comm.norm())