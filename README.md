# A Computational Study: Primal-Dual Hybrid Gradient (PDHG) Algorithm for Linear Programming

The code `my_pdhg.m` implemented a **first-order method** for solving the following standard LP pair:
min c^T x, s.t. Ax = b, x >= 0
max b^T y, s.t. A^T y <= c


## Main Algorithm

### Inilitalization

- Initialize x\textbf{x} and y\textbf{y} arbitarily as long as x≥0\textbf{x}\ge \textbf{0}.
- In my code, I initialized both to be a zero vector.

### Update Rule

- $\textbf{x}^{k+1} = proj_{\mathbb{R}^n_+}(\textbf{x}^k-s(c-A^T\textbf{y}^k))$;
- $\textbf{y}^{k+1}=\textbf{y}^k+s(b-A(2\textbf{x}^{k+1}-\textbf{x}^k))$.

## Algorithm Enhancement

### Scaling

- Set the scale parameter α=‖b‖‖c‖\alpha=\frac{\|\textbf{b}\|}{\|\textbf{c}\|}.
- Do the scaling:
  - At the beginning of the algorithm: c←α⋅c\textbf{c} \leftarrow \alpha \cdot \textbf{c};
  - At the end of the algorithm: y←y/α\textbf{y} \leftarrow \textbf{y}/\alpha.

### Preconditioning

- D1AD2(D−12x)=D1bD_1AD_2(D_2^{-1}\textbf{x}) = D_1\textbf{b}, where D1D1 and D2D2 are diagonal matrices. Correspondingly, $D_2A^TD_1(D_1^{-1}\textbf{y})\le D_2\textbf{c}$.
- We scale the matrix AA to make it well conditioned: A←D1AD2A\leftarrow D_1AD_2.
- Correspondingly, we need to scale b\textbf{b} and cc: b←D1b\textbf{b} \leftarrow D_1\textbf{b} and c←D1c\textbf{c}\leftarrow D_1 \textbf{c}.
- Thus, at the end of the algorithm, we need to scale x\textbf{x} and y\textbf{y} back by x=D2x\textbf{x}=D_2\textbf{x} and y=D1y\textbf{y}=D_1\textbf{y}.
- Refer to the paper proposed by Pock and Chambolle, we set $D_1=diag\{ \frac{1}{\sum_{i=1}^{m}|A_{i,j}|^{2-\beta}}, i=1, \dots, m \}^{-\frac{1}{2}},and, and D_2=diag\{ \frac{1}{\sum_{j=1}^{n}|A_{i,j}|^{\beta}}, j=1, \dots, n \}^{-\frac{1}{2}}.Inthecode,Iset. In the code, I set \beta=1.2$.

### Adaptive step

- Set the initial stepsize s=γ‖A‖s=\frac{\gamma}{\|A\|} after scaling AA in the precoditioning part, where γ∈(0,1)\gamma\in(0, 1). I set γ=0.99\gamma=0.99 in my code to be greedy.
- After every iteration, we adapt the stepsize $\tilde{s}^{k+1} = max\{s, t\cdot\frac{\|d\textbf{x}^{k+1}\|}{\|Ad\textbf{x}^{k+1}\|}\},where, where d\textbf{x}^{k+1}=\textbf{x}^{k+1}-\textbf{x}^{k}and and tistheshrinkparameterin is the shrink parameter in (0, 1)$.

### Restart

- Set N=min{3⋅n,maxit2}N=min\{3\cdot n, \frac{maxit}{2}\} and do the restart for every NN iterations from the average x\textbf{x} and y\textbf{y} of the past 1000 iterations.

## Implementation Details

The code is highly optimized by the following tricks:

- **Only calculate necessary errors among the three.** If the primal relative error is not less than tolerance, and the error information is not needed to printe, then other two errors are not calculated. If the primal relative error is less than the tolerance, then we calculate the relative gap error and see whether we need to calculate the relative dual error.
- **Only calculate AxAx and ATyA^Ty once per iteration**. These are the bottleneck for the algorithm efficiency in every iteration and nothing can be done to avoid these two multiplications between large-scale matrix and vector.

## Observations and Code Performance

- My code achieves the efficiency of the instructor's code for all but one cases (i.e. for the case m=400, n=4000, mine is slower than his by 0.09s. for the others, mine is much faster). The code maybe further fastered by choosing better hyperparameters.
- My obversation for the code behavior under different enhancement for different problems are as follows:
  - The **scaling and adaptive step** enhancement are crucial for all the problems for good performance.
  - For the **random problem of small size with very low tolerance** 1e−141e-14, the **restarting** enhancement is crucial to reach high accuracy. 
  - For the **random problem of midium size with midium tolerance 1e-6**, the algorithm costs a lot of time since the matrices are much larger.
  - For the **image impainting problem**, I turn off the **restarting** part since it normally leads to longer elapsed time. Restarting is not necessary for this problem because the tolerance 1e−21e-2 is so large that the algorithm can easily converges.
  - For the **Israel problem**, the **preconditioning** enhancement is very crucial since the matrix AA are special (maybe ill-conditioned or high-unbalanced elements).
