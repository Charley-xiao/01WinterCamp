# 计算流体力学概述

## 1 有限差分法

**有限差分法**：基于待求解偏微分方程的微分形式，通过有限差分来近似导数，从而寻求微分方程的近似解.

差分格式构建方法：**泰勒展开法**、**多项式逼近法**.

### 1.1 泰勒展开

使用泰勒展开法，对于物理量 $\phi(x)$，有：
$$
\phi(x)=\phi(x_i)+\sum\limits_{i=1}^n\dfrac{(x-x_i)^n}{i!}\left(\dfrac{\part^i\phi}{\part x^i}\right)_i+H
$$
在 $x_{i+1}$ 处，有：
$$
\left(\dfrac{\part\phi}{\part x}\right)_i=\\\dfrac{\phi_{i+1}-\phi_i}{x_{i+1}-x_i}-\dfrac{x_{i+1}-x_i}2\left(\dfrac{\part^2\phi}{\part x^2}\right)_i-\dfrac{(x_{i+1}-x_i)^2}6\left(\dfrac{\part^3\phi}{\part x^3}\right)_i+H\tag{1}
$$
**前向差分**：
$$
\left(\dfrac{\part\phi}{\part x}\right)_i\approx\dfrac{\phi_{i+1}-\phi_i}{x_{i+1}-x_i}
$$
在 $x_{i-1}$ 处，有：
$$
\left(\dfrac{\part\phi}{\part x}\right)_i=\\\dfrac{\phi_{i}-\phi_{i-1}}{x_{i}-x_{i-1}}+\dfrac{x_{i}-x_{i-1}}2\left(\dfrac{\part^2\phi}{\part x^2}\right)_i-\dfrac{(x_{i}-x_{i-1})^2}6\left(\dfrac{\part^3\phi}{\part x^3}\right)_i+H\tag{2}
$$
**后向差分**：
$$
\left(\dfrac{\part\phi}{\part x}\right)_i\approx\dfrac{\phi_{i}-\phi_{i-1}}{x_{i}-x_{i-1}}
$$
结合 $(1),(2)$ 可得**中心差分**：
$$
\left(\dfrac{\part\phi}{\part x}\right)_i=\\
\dfrac{\phi_{i+1}-\phi_{i-1}}{x_{i+1}-x_{i-1}}-\dfrac{(x_{i+1}-x_i)^2-(x_i-x_{i-1})^2}{2(x_{i+1}-x_{i-1})}\left(\dfrac{\part^2\phi}{\part x^2}\right)_i\\-\dfrac{(x_{i+1}-x_i)^3-(x_{i}-x_{i-1})^3}{6(x_{i+1}-x_{i-1})}\left(\dfrac{\part^3\phi}{\part x^3}\right)_i+H\\\approx\\
\dfrac{\phi_{i+1}-\phi_{i-1}}{x_{i+1}-x_{i-1}}
$$

### 1.2 多项式逼近

使用多项式逼近法，有以下步骤：

1. 确定差分网格节点；
2. 选择多项式函数分布；
3. 计算多项式函数；
4. 对多项式求导.

假如我们使用抛物线拟合，需要三个点 $x_{i+1},x_i,x_{i-1}$.
$$

$$
