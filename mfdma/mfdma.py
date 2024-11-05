emimport numpy as np

def width(alpha):
    return float(np.max(alpha) - np.min(alpha))

def assy(f,alpha):
    i = np.argmax(f)
    a_0 = alpha[i]
    return float((np.max(alpha) - a_0)/ (a_0 - np.min(alpha)))   

def MFDMA_1D(x, n_min, n_max, N, theta, q):

    """
    Performs Multifractal Detrended Moving Average (MF-DMA) analysis on a 1D time series.

    Parameters:
    -----------
    x : array-like, shape (M,)
        Input time series data. If the input is not 1D, it will be reshaped.
    n_min : int
        Minimum scale window length for the log-spaced window sizes.
    n_max : int
        Maximum scale window length for the log-spaced window sizes.
    N : int
        Number of scales to consider between n_min and n_max.
    theta : float
        Offset parameter controlling the detrending method. Typically 0 <= theta <= 1.
    q : array-like
        Array of moments for the multifractal analysis, including both positive and negative values.
        For q = 0, a logarithmic averaging is used.

    Returns:
    --------
    n : array-like
        The window sizes used for the analysis.
    Fq : array-like
        The fluctuation function values for each window size and each q.
    tau : array-like
        The multifractal scaling exponent for each q.
    alpha : array-like
        The singularity strength values.
    f : array-like
        The multifractal spectrum f(α), indicating the distribution of singularities.

    Description:
    ------------
    The function applies the MF-DMA method by computing the cumulative sum of the input series, dividing it into segments, and calculating the root-mean-square fluctuation for each segment. It then estimates the scaling exponent tau(q) from the fluctuation function across different window sizes, and derives the singularity strength alpha(q) and multifractal spectrum f(α).
    """
    # Verifica se x é uma série temporal unidimensional
    if x.ndim == 1:
        x = x.reshape(-1, 1)
    
    M = len(x)
    MIN = np.log10(n_min)
    MAX = np.log10(n_max)
    n = np.unique(np.round(np.logspace(MIN, MAX, N))).astype(int)

    # Calcula a soma cumulativa y
    y = np.cumsum(x)
    
    F = []
    for i in range(len(n)):
        lgth = n[i]
        y1 = np.array([np.mean(y[j:j + lgth]) for j in range(M - lgth + 1)])

        # Calcula o residual e, assegurando que os tamanhos sejam compatíveis
        start_idx = max(0, int(np.floor(lgth * (1 - theta))))
        e = y[start_idx:start_idx + len(y1)] - y1[:len(y[start_idx:start_idx + len(y1)])]

        # Estima a função root-mean-square F
        F_i = [np.sqrt(np.mean(e[k * lgth:(k + 1) * lgth] ** 2)) for k in range(len(e) // lgth)]
        F.append(F_i)

    Fq = np.zeros((len(F), len(q)))
    for i in range(len(q)):
        for j in range(len(F)):
            f = np.array(F[j])
            if q[i] == 0:
                Fq[j, i] = np.exp(0.5 * np.mean(np.log(f ** 2)))
            else:
                Fq[j, i] = (np.mean(f ** q[i])) ** (1 / q[i])

    # Calcula o expoente de escalamento multifractal tau(q)
    h = np.zeros(len(q))
    for i in range(len(q)):
        f_q = Fq[:, i]
        k = np.polyfit(np.log(n), np.log(f_q), 1)[0]
        h[i] = k
    tau = h * q - 1

    # Calcula a força de singularidade alpha(q) e o espectro multifractal f(alpha)
    dx = 7
    dx = (dx - 1) // 2
    valid_length = len(tau) - 2 * dx

    # Verifica se valid_length é positivo
    if valid_length <= 0:
        raise ValueError("O comprimento válido para o cálculo de alpha é negativo ou zero. "
                     "Verifique os parâmetros de entrada.")

    alpha = np.zeros(valid_length)
    for i in range(dx, len(tau) - dx):
        xx = q[i - dx:i + dx + 1]
        yy = tau[i - dx:i + dx + 1]
        alpha[i - dx] = np.polyfit(xx, yy, 1)[0]

    alpha = alpha[:valid_length]
    q_cut = q[dx:dx + valid_length]
    f = q_cut * alpha - tau[dx:dx + valid_length]


    return n, Fq, tau, h, alpha, f


def MFDMA_2D(X, n_min, n_max, N, theta, q):
    """
    Perform Multifractal Detrended Moving Average Analysis (MFDMA) on 2D data.

    Parameters
    ----------
    X : numpy.ndarray
        2D array representing the time series data to be analyzed.
    n_min : int
        Minimum window size for the moving average.
    n_max : int
        Maximum window size for the moving average.
    N : int
        Number of window sizes between `n_min` and `n_max` to consider.
    theta : float
        Parameter controlling the overlap between windows.
    q : list or numpy.ndarray
        List of q-order moments used to compute the multifractal spectrum.

    Returns
    -------
    n : numpy.ndarray
        Array of window sizes used in the analysis.
    Fq : numpy.ndarray
        Matrix containing the fluctuation functions for each window size and each q value.
    tau : numpy.ndarray
        Multifractal scaling exponents, related to the Hurst exponent.
    alpha : numpy.ndarray
        Singularity strengths (or Holder exponents) for each q.
    f : numpy.ndarray
        Multifractal spectrum (f(alpha)), providing information on the fractal dimensions of the structure.
    """
    
    N1, N2 = X.shape
    MIN = np.log10(n_min)
    MAX = np.log10(n_max)
    n = np.unique(np.round(np.logspace(MIN, MAX, N))).astype(int)
    
    F = []
    for lgth in n:
        Y = np.zeros((N1 - lgth + 1, N2 - lgth + 1))
        Y1 = np.zeros((N1 - lgth + 1, N2 - lgth + 1))
        
        for j in range(N1 - lgth + 1):
            for k in range(N2 - lgth + 1):
                Z = X[j:j + lgth, k:k + lgth]
                Z1 = np.cumsum(np.cumsum(Z, axis=0), axis=1)
                Y[j, k] = Z1[-1, -1]
                Y1[j, k] = np.mean(Z1)
        
        # Determine the residual e
        x0 = np.arange(0, Y.shape[0] - min(int(np.floor(lgth * theta)), lgth - 1))
        y0 = np.arange(0, Y.shape[1] - min(int(np.floor(lgth * theta)), lgth - 1))
        x1 = np.arange(Y1.shape[0] - len(x0), Y1.shape[0])
        y1 = np.arange(Y1.shape[1] - len(y0), Y1.shape[1])
        e = Y[np.ix_(x0, y0)] - Y1[np.ix_(x1, y1)]
        
        # Estimate the root-mean-square function F
        F_i = np.zeros((e.shape[0] // lgth, e.shape[1] // lgth))
        for k1 in range(F_i.shape[0]):
            for k2 in range(F_i.shape[1]):
                E = e[k1 * lgth:(k1 + 1) * lgth, k2 * lgth:(k2 + 1) * lgth]
                F_i[k1, k2] = np.sqrt(np.mean(E ** 2))
        F.append(F_i)
    
    Fq = np.zeros((len(F), len(q)))
    for i in range(len(q)):
        for j in range(len(F)):
            f = F[j].ravel()
            if q[i] == 0:
                Fq[j, i] = np.exp(0.5 * np.mean(np.log(f ** 2)))
            else:
                Fq[j, i] = (np.mean(f ** q[i])) ** (1 / q[i])
    
    # Calculate the multifractal scaling exponent tau(q)
    h = np.zeros(len(q))
    for i in range(len(q)):
        fq = Fq[:, i]
        k = np.polyfit(np.log(n), np.log(fq), 1)[0]
        h[i] = k
    tau = h * q - 2
    
    # Calculate the singularity strength function alpha(q) and spectrum f(alpha)
    dx = 7
    dx = (dx - 1) // 2
    alpha = np.zeros(len(tau) - 2 * dx)
    for i in range(dx, len(tau) - dx):
        xx = q[i - dx:i + dx + 1]
        yy = tau[i - dx:i + dx + 1]
        alpha[i - dx] = np.polyfit(xx, yy, 1)[0]
    
    alpha = alpha[dx:len(alpha)-dx]
    q_cut = q[dx:len(alpha)+dx]
    f = q_cut * alpha - tau[dx:dx + len(q_cut)]
    
    return n, Fq, tau, alpha, f

