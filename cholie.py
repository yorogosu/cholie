def cholie(g, normalize=False):
    """
    Calculate the Cholesky decomposition of a GRM based on a genotype matrix g

    This function simply shuffles the rows in g until no error is returned
    """
    n_snps, n_indiv = np.shape(g)
    resolved = False
    while not resolved:
        try:
            h = pd.DataFrame.sample(g, n_snps)
            grm = sp.dot(h.T, h) / n_snps
            cho = cholesky(pinv(grm))
            resolved = True
            return cho
        except:
            pass
