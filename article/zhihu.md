# Quantum Mpemba Effect in Closed Quantum Systems

为了更加直观地研究量子 Mpemba 效应（Quantum Mpemba Effect, QME），我们关注封闭量子系统中的非平衡动力学。

传统的 Mpemba 效应通常讨论与外部环境耦合的开放系统，其中具有不同初始温度的状态最终弛豫到相同的热平衡态。其反常之处在于：一个初始时刻距离最终平衡态更远的状态，反而可能更早到达平衡态。

对于封闭量子系统，情况有所不同：系统与外界之间不存在能量或信息交换，其动力学由幺正演化描述：

$$
|\psi(t)\rangle=e^{-iHt}|\psi(0)\rangle .
$$

因此，在封闭量子系统中，Mpemba 效应不能简单地通过“系统耗散并达到最终平衡态所需要的时间”来定义。相反，我们考察的是**局域对称性的恢复（local symmetry restoration）**。

首先，我们制备一系列具有不同倾斜角的**倾斜铁磁态（tilted ferromagnetic states）**。不同的倾斜角代表对哈密顿量所具有的 $U(1)$ 自旋旋转对称性的不同破缺程度。

在量子淬火（quantum quench）之后，尽管整个封闭系统始终进行幺正演化，但其子系统仍然可以发生局域弛豫，并逐渐恢复这种对称性。

在这一背景下，QME 表现为：

**初始对称性破缺程度更大的状态，反而比初始对称性破缺程度较小的状态更快地实现局域对称性恢复。**

我们主要研究带有 $z$ 方向外场的 XY 模型：

$$
H =
\sum_{i<j}\frac{J_{ij}}{2}
\left(
\sigma_i^x\sigma_j^x+\sigma_i^y\sigma_j^y
\right)
+
\sum_i h_i\sigma_i^z .
$$


## 1. Tilted Ferromagnetic State

在本文中，**倾斜铁磁态（tilted ferromagnetic state）**可以明确地写成一个张量积态。

从完全极化态出发：

$$
|\uparrow\uparrow\cdots\uparrow\rangle
=
|\uparrow\rangle^{\otimes N}.
$$

将每一个自旋旋转相同的极角 $\theta$ 和方位角 $\phi$，定义

$$
|TF(\theta,\phi)\rangle
=
e^{-i\phi \hat M_z}
e^{-i\theta \hat M_y}
|\uparrow\rangle^{\otimes N}.
$$

因此，

$$
|TF(\theta,\phi)\rangle
=
\bigotimes_{j=1}^{N}
\left(
\cos\frac{\theta}{2}|\uparrow\rangle_j
+
e^{i\phi}\sin\frac{\theta}{2}|\downarrow\rangle_j
\right).
$$

对于 $\phi=0$ 的情况，

$$
|TF(\theta,0)\rangle
=
\left(
\cos\frac{\theta}{2}|\uparrow\rangle
+
\sin\frac{\theta}{2}|\downarrow\rangle
\right)^{\otimes N}.
$$


![Initial state and tilted state](E:/physics%20notes/condensed%20matter%20physics/quantum%20mpemba%20effect/Weixin%20Image_20260911133309_320_7.jpg)

初始态与倾斜态的简单示意图。对于倾斜态，其磁化方向在 Bloch 球上发生了旋转。


## 2. Entanglement Asymmetry

对于一个子系统 $A$，我们定义纠缠不对称度（entanglement asymmetry, EA）：

$$
\Delta S_A
=
\ln\left[\mathrm{Tr}_A(\rho_A^2)\right]
-
\ln\left[\mathrm{Tr}_A(\rho_{\tilde A}^2)\right].
$$

其中 $\rho_{\tilde A}$ 定义为

$$
\rho_{\tilde A}
=
\sum_m
\Pi_m \rho_A \Pi_m ,
$$

其中 $\Pi_m$ 是 $M_z^A$ 对应于本征值 $m$ 的本征子空间上的投影算符。

EA 描述了倾斜态偏离 $U(1)$ 对称性的程度。EA 只有在

$$
[\rho_A,M_z^A]=0
$$

时才会消失。

因此，我们可以通过研究

$$
\Delta S_A(t)
$$

的时间演化来判断局域 $U(1)$ 对称性的恢复过程。

如果一个初始 EA 更大的状态，其 EA 随时间衰减得反而更快，并最终与初始 EA 更小的状态发生交叉，那么这种反常的对称性恢复动力学就是我们所研究的 QME。


# 3. Macroscopic Simulations

## 3.1 Power-Law Coupling

首先考虑幂律衰减的自旋耦合。不同位置的自旋之间的耦合强度随距离按照幂律发生变化。


![EA dynamics](E:/physics%20notes/condensed%20matter%20physics/quantum%20mpemba%20effect/screenshot007.png)

EA 的动力学。由于系统的时间演化是幺正的，不同曲线会发生多次交叉并重新靠近。因此，我们主要关注早期的 $0<J_0t<3$ 动力学阶段。


![Finite-size effect](E:/physics%20notes/condensed%20matter%20physics/quantum%20mpemba%20effect/screenshot016.png)

通过改变系统尺寸 $N$ 对 EA 动力学进行比较，以排除有限尺寸效应（finite-size effect）。


![Angular distribution](E:/physics%20notes/condensed%20matter%20physics/quantum%20mpemba%20effect/screenshot012.png)

对于 $N=12$ 的自旋链，随着对称性逐渐恢复，角分布之间的差异逐渐被抹平。


## 3.2 Husimi $Q$ Function

为了更加直观地观察局域 $U(1)$ 对称性的恢复，可以进一步研究 Husimi $Q$ function。

对于子系统密度矩阵 $\rho_A(t)$，定义

$$
Q(\theta,\phi;t)
=
\langle\theta,\phi|
\rho_A(t)
|\theta,\phi\rangle,
$$

其中 $|\theta,\phi\rangle$ 表示相应的自旋相干态。

Husimi $Q$ function 描述了子系统量子态在自旋相干态相空间中的角分布。


![Husimi Q distribution](E:/physics%20notes/condensed%20matter%20physics/quantum%20mpemba%20effect/screenshot014.png)

$Q(\theta,\phi;t)$ 给出了演化后的子系统密度矩阵在自旋相干态 $|\theta,\phi\rangle$ 上的投影。


![Ideal Q distribution](E:/physics%20notes/condensed%20matter%20physics/quantum%20mpemba%20effect/screenshot015.png)

理想情况下，随着 $U(1)$ 对称性逐渐恢复，$Q$ function 在方位角 $\phi$ 方向上的不均匀性逐渐消失，最终趋向更加均匀的角分布。

因此，从 Bloch 球上的几何图像来看，对称性恢复可以理解为：

$$
\text{azimuthally nonuniform distribution}
\quad\longrightarrow\quad
\text{azimuthally uniform distribution}.
$$

也就是说，随着时间演化，初始态在 $\phi$ 方向上的结构逐渐被“抹平”。


## 3.3 Dependence on Power-Law Exponent

接下来改变幂律耦合中的指数 $\alpha$，研究相互作用范围对于 QME 的影响。


![EA with different alpha](E:/physics%20notes/condensed%20matter%20physics/quantum%20mpemba%20effect/screenshot006.png)

不同幂律指数 $\alpha$ 下的 EA$(t)$。改变 $\alpha$ 会改变自旋之间有效相互作用的空间范围，从而改变局域对称性恢复的动力学。


# 4. Exponentially Decaying Coupling

除了幂律相互作用，我们进一步考虑指数衰减形式的耦合。

其中，特征长度 $\xi$ 控制相互作用能够延伸的空间范围。


![xi=0.5](E:/physics%20notes/condensed%20matter%20physics/quantum%20mpemba%20effect/screenshot008.png)

对于 $\xi=0.5$ 的情况，在有限的时间范围内我们没有观察到明显的 QME。


![xi=2](E:/physics%20notes/condensed%20matter%20physics/quantum%20mpemba%20effect/screenshot009.png)

对于 $\xi=2$ 的情况，QME 重新出现。


![QME versus characteristic length](E:/physics%20notes/condensed%20matter%20physics/quantum%20mpemba%20effect/screenshot010.png)

随着特征长度 $\xi$ 增大，QME 重新出现，同时曲线的交叉时间提前。

这一结果说明，相互作用的空间范围可能是影响 QME 的一个重要参数。

对于短程相互作用系统中的 QME，目前存在一些不同的物理解释，例如可以从电荷输运（charge transport）的图像理解对称性恢复过程。


# 5. Disorder

最后，我们进一步加入 $z$ 方向的无序场，研究无序对于 QME 的影响。


![Disorder](E:/physics%20notes/condensed%20matter%20physics/quantum%20mpemba%20effect/screenshot011.png)

在当前模拟参数下，$z$ 方向的无序破坏了原本观察到的 QME。


# 6. Coupling Range, Entanglement and QME

综合上述数值模拟可以看到，QME 与相互作用的空间结构密切相关。

改变幂律指数 $\alpha$ 或指数衰减特征长度 $\xi$，本质上都是在改变体系中自旋之间的有效相互作用范围。这种变化进一步影响量子信息、关联以及纠缠在系统中的传播，并最终改变子系统局域对称性恢复的动力学。

因此，可以将目前观察到的物理关系概括为

$$
\text{coupling range}
\longrightarrow
\text{correlation / entanglement dynamics}
\longrightarrow
\text{local symmetry restoration}
\longrightarrow
\text{QME}.
$$

需要注意的是，这并不意味着简单的

$$
\text{stronger coupling}
\Rightarrow
\text{stronger entanglement}
\Rightarrow
\text{stronger QME}.
$$

更重要的是**相互作用的空间结构以及不同初态在这种相互作用下产生的动力学差异**。

在我们的指数衰减耦合模拟中，当 $\xi$ 很小时，相互作用高度局域，在有限时间窗口内没有观察到明显的 EA crossover；随着 $\xi$ 增大，相互作用能够覆盖更远的自旋，QME 重新出现，而且 crossover time 提前。

因此，目前的数值结果表明：增加有效相互作用范围能够促进在有限时间窗口内观察到 QME。

但仅仅观察 EA crossover 还不足以确定这一现象的微观机制。为了进一步理解这一结果，需要研究纠缠增长、关联传播以及守恒荷输运与 EA 恢复速度之间的关系。


# 7. Summary

在封闭量子系统中，由于整体动力学始终保持幺正演化，因此 QME 并不是通过系统到达某个最终热平衡态的时间来定义，而是通过**局域对称性恢复的反常速度**来定义。

对于不同倾斜角 $\theta$ 的 tilted ferromagnetic states，初始倾斜角决定了 $U(1)$ 对称性的破缺程度。通过研究子系统的 entanglement asymmetry

$$
\Delta S_A(t),
$$

可以直接追踪局域对称性的恢复过程。

QME 的核心特征是：

$$
\text{greater initial symmetry breaking}
\quad\longrightarrow\quad
\text{faster symmetry restoration}.
$$

从数值模拟可以看到，相互作用的空间范围对这一过程具有明显影响。

对于幂律耦合，改变 $\alpha$ 会改变 QME 的动力学；对于指数衰减耦合，增加特征长度 $\xi$ 可以使有限时间内原本消失的 QME 重新出现，并使 EA crossover 提前。

与此同时，Husimi $Q$ function 为这一过程提供了更加直观的几何图像：随着局域 $U(1)$ 对称性恢复，量子态在方位角 $\phi$ 方向上的不均匀分布逐渐被抹平。

因此，目前的模拟结果揭示了一条值得进一步研究的关系：

$$
\text{interaction range}
\rightarrow
\text{information / entanglement propagation}
\rightarrow
\text{symmetry restoration}
\rightarrow
\text{Quantum Mpemba Effect}.
$$