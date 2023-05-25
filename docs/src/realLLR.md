# Likelihood ratio tests


Consider correlated genes $A$, $B$, and a third variable $E$ upstream of $A$ and $B$, such as a significant eQTL of $A$. The eQTLs can be obtained either *de novo* using eQTL identification tools such as matrix-eQTL [@Shabalin:2012] or kruX [@Qi:2014], or from published analyses. Throughout this article, we assume that $E$ is a significant eQTL of $A$, whereas extension to other data types is straightforward. We use $A_i$ and $B_i$ for the expression levels of gene $A$ and $B$ respectively, which are assumed to have gone through the supernormalization, and optionally the genotypes of the best eQTL of $A$ as $E_i$, where $i=1,\dots,n$ across samples. Genotypes are assumed to have a total of $n_a$ alleles, so $E_i\in\{0,\dots,n_a\}$. We define the null and alternative hypotheses for a total of six tests, as shown in **Table [1](#tab-tests){reference-type="ref" reference="tab-tests"}**. LLRs of every test are calculated separately as follows:


## Correlation test

Define the null hypothesis as $A$ and $B$ are independent, and the alternative hypothesis as they are correlated:

$${\mathcal H}_{\mathrm{null}}^{\mathrm{(0)}}=A\qquad B,\hspace{4em}{\mathcal H}_{\mathrm{alt}}^{\mathrm{(0)}}=A \to B.$$

The superscript $(0)$ is the numbering of the test. Both hypotheses are modeled with gene expression levels following bivariate normal distributions, as

$$\begin{pmatrix}
    A_i\\\\
    B_i
\end{pmatrix} \sim 
N\left(
    \begin{pmatrix}
        0\\\\
        0
    \end{pmatrix},
    \begin{pmatrix}
        \sigma_{A0}^2 & \rho\,\sigma_{A0}\sigma_{B0}\\\\
        \rho\,\sigma_{A0}\sigma_{B0} & \sigma_{B0}^2
    \end{pmatrix}
\right)$$

for $i=1,\dots,n$. The null hypothesis corresponds to $\rho=0$.

Maximum likelihood estimators (MLE) for the model parameters $\rho$, $\sigma_{A0}$, and $\sigma_{B0}$ are

$$\hat\rho=\frac{1}{n}\sum_{i=1}^nA_iB_i,\hspace{2em}\hat\sigma_{A0}=\hat\sigma_{B0}=1,\label{eq-cor-rho}$$

and the LLR is simply

$$\mathrm{LLR}^{\mathrm{(0)}}=-\frac{n}{2}\ln(1-\hat\rho^2).\label{eq-cor-llr}$$

In the absence of genotype information, we use nonzero correlation between $A$ and $B$ as the indicator for $A\rightarrow B$ regulation, giving the posterior probability

$$P(A \to B)=P({\mathcal H}_{\mathrm{alt}}^{\mathrm{(0)}}\mid\mathrm{LLR}^{\mathrm{(0)}}).\nonumber$$

false negative

1.  **Primary (linkage) test:** Verify that $E$ regulates $A$ from ${\mathcal H}_{\mathrm{alt}}^{\mathrm{(1)}}\equiv E\rightarrow A$ and ${\mathcal H}_{\mathrm{null}}^{\mathrm{(1)}}\equiv E\hspace{1.6em}A$. For ${\mathcal H}_{\mathrm{alt}}^{\mathrm{(1)}}$, we model $E\rightarrow A$ as $A$ follows a normal distribution whose mean is determined by $E$ categorically, i.e.

$$A_i\mid E_i\sim N(\mu_{E_i},\sigma_A^2).\label{eq-ba-hm-d1}$$ 

From the total likelihood $p(A\mid E)=\prod_{i=1}^np(A_i\mid E_i)$, we find MLEs for model parameters $\mu_j,j=0,1,\dots,n_a$, and $\sigma_A$, as


$$\hat\mu_j=\frac{1}{n_j}\sum_{i=1}^nA_i\delta_{E_ij},\hspace{2em}\hat\sigma_A^2=1-\sum_{j=0}^{n_a}\frac{n_j}{n}\hat\mu_j^2,\nonumber$$

where $n_j$ is the sample count by genotype category,

$$n_j\equiv\sum_{i=1}^n\delta_{E_ij}.\nonumber$$ 

The Kronecker delta function is defined as $\delta_{xy}=1$ for $x=y$, and 0 otherwise. When summing over all genotype values ($j=0,\dots,n_a$), we only pick those that exist ($n_j>0$) throughout this article. Since the null hypothesis is simply that $A_i$ is sampled from a genotype-independent normal distribution, with MLEs of mean zero and standard deviation one due to the supernormalization (**Section [\[ssec-algorithm\]](#ssec-algorithm){reference-type="ref"    reference="ssec-algorithm"}**), the LLR for test 1 becomes

$$\mathrm{LLR}^{\mathrm{(1)}}=-\frac{n}{2}\ln\hat\sigma_A^2.\label{eq-ba-hm-llr1}$$

By favoring a large $\mathrm{LLR}^{\mathrm{(1)}}$, we select  ${\mathcal H}_{\mathrm{alt}}^{\mathrm{(1)}}$ and verify that $E$ regulates $A$, with
    
$$P(E\rightarrow A)=P({\mathcal H}_{\mathrm{alt}}^{\mathrm{(1)}}\mid\mathrm{LLR}^{\mathrm{(1)}}).\nonumber$$
