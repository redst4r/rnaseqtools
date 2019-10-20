import multiprocessing as mp
import logging
import numpy as np
from scipy.sparse import dok_matrix, csc_matrix
import tqdm

def kallisto_bootstrap(ec_counts, ec_dict, transcript_lengths, n_bootstraps, n_iter, n_cores=4):

    assert ec_counts.ndim==2
    assert ec_counts.shape[1] == 1

    def parallel_em(x, ec, trans_len):
        EM = KallistoEM(x, ec, trans_len)
        rho_list, lopP_list = EM.em_iterator(n_iter, rho_init=None, reportfreq = 20)
        return rho_list[-1,:]


    args = []
    for i in range(n_bootstraps):
        p = ec_counts.flatten()/ec_counts.sum()
        counts_boot = np.random.multinomial(np.sum(ec_counts), p, size=1).flatten().reshape(-1,1)
        args.append((counts_boot, ec_dict, transcript_lengths))

    with mp.Pool(n_cores) as pool:
        bootstrapped_rhos = pool.starmap(parallel_em, args)



class KallistoEM(object):

    def __init__(self, ec_counts, ec_dict, transcript_lengths):
        self.ec_counts = ec_counts
        self.ec_dict = ec_dict
        self.transcript_lengths = transcript_lengths

        self.n_transcripts = len(transcript_lengths)
        self.n_ec = len(ec_dict)
        assert len(ec_counts) == self.n_ec

        # turn it from a dict into a sparse matrix
        logging.info('Calculating EC-matrix')
        self.ec_matrix = ec_dict_2_ecmatrix_csc(self.ec_dict, self.n_ec, self.n_transcripts)

    def EMstep(self, rho):
        """
            rho: the relative abundances of the transcripts  (irrspecitve of size)
            ec_counts: counts in the equvalence classes

        """

        rho_slack = 1e-15
        rho = rho + rho_slack
        rho /= rho.sum()

        # here we distribute the reads over the eq.classes
        # since this depends on the size of the transcripts,
        # it will be biased by length (hence we call it alpha, not rho)

        tr_alphas = self.ec_matrix.multiply(rho)

        #    tr_alphas = tr_alphas/ tr_alphas.sum(1):
        tr_alphas = tr_alphas.multiply(1/tr_alphas.sum(1))

        # distributing the reads over the ECs according to their weight
        count_update = tr_alphas.multiply(self.ec_counts)

        # counts per transcript, summing up contributions from all ECs
        counts_per_alpha = count_update.sum(0)


        new_alpha = counts_per_alpha/ counts_per_alpha.sum()

        # now correct for the length bias!
        new_rho = new_alpha / self.transcript_lengths
        new_rho = new_rho / new_rho.sum()

        new_rho = new_rho.A.flatten()
        new_alpha = new_alpha.A.flatten()

        return new_rho

    def loglikelihood(self, rho):
        logp = np.zeros_like(self.ec_counts).astype('float64')
        rho_over_t =  rho / self.transcript_lengths

        p_tmp = self.ec_matrix.multiply(rho_over_t).sum(1)

        logP_tmp = np.log(p_tmp).A
        logp = self.ec_counts * logP_tmp
        logp[self.ec_counts == 0] = 0

        return logp.sum()

    def em_iterator(self, n_iter, rho_init=None, reportfreq = 20):
        if not isinstance(rho_init, np.ndarray):
            rho = np.ones_like(self.transcript_lengths) / len(self.transcript_lengths) # uniform dist
        else:
            rho = rho_init

        logging.info("starting EM iterations")

        lopP_list = [self.loglikelihood(rho)]
        rho_list = [rho]
        for i in range(n_iter):

            rho = self.EMstep(rho)
            logP = self.loglikelihood(rho)

            if i % reportfreq == 0:
                non_zero_rho = np.sum(rho > 0)
                print(f'Iter {i}:\tLogP={logP} (rho>0: {non_zero_rho})')
                # print(np.exp(logP), rho, a)

            lopP_list.append(logP)
            rho_list.append(rho)

        rho_list = np.stack(rho_list)
        logp_list = np.array(lopP_list)

        return rho_list, lopP_list


def ec_dict_2_ecmatrix(ec_dict, n_ec, n_transcripts):
    logging.info("building dok-matrix")
    S = dok_matrix((n_ec, n_transcripts), dtype=np.int32)
    for i, vals in tqdm.tqdm(ec_dict.items()):
        for j in vals:
            S[i, j] = 1
    return S

def ec_dict_2_ecmatrix_csc(ec_dict, n_ec, n_transcripts):
    """
    converts the EC-comatibility dict directly into a csc-sparse matrix
    """
    logging.info("building csc-matrix")
    row_col = []
    for i,vals in tqdm.tqdm(ec_dict.items()):
        update = [(i, _) for _ in vals]
        row_col.extend(update)
    rows, cols = zip(*row_col)

    data = np.ones(len(rows))

    S = csc_matrix((data, (rows, cols)), shape=(n_ec, n_transcripts))
    return S


# def kallisto_EM_update(rho_vector, ec_vector, ec_2_tr, transcript_length):
#     """
#         rho_vector: the relative abundances of the transcripts  (irrspecitve of size)
#         ec_vector: counts in the equvalence classes
#
#     """
#
#     N_total = ec_vector.sum()
#     # here we distribute the reads over the eq.classes
#     # since this depends on the size of the transcripts, it will be biased by length (hence we call it alpha, not rho)
#
#     counts_per_alpha = np.zeros_like(rho_vector).astype('float64')
#     for i in range(len(ec_vector)):
#         # what transcript are associated
#
#         # if no counts are to be distributed, just skip, otherise we get 0/0 divisions
#         if ec_vector[i] == 0:
#             continue
#
#         ix_tr = ec_2_tr[i]
#
#         tr_alphas = rho_vector[ix_tr]# this determines how we distribute the counts in EC
#         tr_alphas = tr_alphas / tr_alphas.sum()
#
#         count_update = tr_alphas * ec_vector[i]
#         counts_per_alpha[ix_tr] = counts_per_alpha[ix_tr] + count_update
#
#     np.testing.assert_approx_equal(counts_per_alpha.sum(), N_total)
#
#     new_alpha = counts_per_alpha/ counts_per_alpha.sum()
#
#     # now correct for the length bias!
#     new_rho = new_alpha / transcript_length
#     new_rho = new_rho / new_rho.sum()
#
#     return new_rho, new_alpha
#
# def likelihood(rho, ec_vector, ec_2_tr, transcript_length):
#     logp = np.zeros_like(ec_vector).astype('float64')
#
#     rho_over_t =  rho / transcript_length
#
#     for i in range(len(ec_vector)):
#
#         # if ec_vector[i] is zero, and alhpa is zero, we will get a nan
#         # but since there's zero observations, this doesnt matter
#         # hence set logP[i] = 0 if ec_vector[i]
#         if ec_vector[i] == 0:
#             continue
#
#         ix_tr = ec_2_tr[i]
#         p_tmp = np.sum(rho_over_t[ix_tr])
#         # p_tmp = p_tmp ** ec_vector[i]
#         logp[i] = ec_vector[i] * np.log(p_tmp)
#
#     # logp[ec_vector==0] = 0
#     return logp.sum()

# def plot_results(lopP_list, rho_list):
#     import matplotlib.pyplot as plt
#     plt.figure()
#     plt.plot(lopP_list)
#     plt.xlabel('iter')
#     plt.ylabel('logP')
#
#     plt.figure()
#     plt.plot(rho_list)
#     plt.xlabel('iter')
#     plt.ylabel('rho')
#
#
# if __name__ == '__main__':
#     # main()
#
#     rho = np.array([1,1,1]) / 3
#     transcript_length = np.array([109, 23, 17])
#     ec_vector = np.array([1, 2, 1])
#
#     # ec_2_tr = {0: [0,1],
#     #            1: [1,2],
#     #            2: [2],
#     #            3: [0]}
#     ec_2_tr = {0: [0,1],
#                1: [0,1,2],
#                2: [0]}
#
#     print(rho)
#     lopP_list = []
#     rho_list = []
#     for i in range(200):
#         rho, a = kallisto_EM_update(rho, ec_vector, ec_2_tr, transcript_length)
#         # print(rho)
#         logP = likelihood(rho, ec_vector, ec_2_tr, transcript_length)
#         print(np.exp(logP), rho, a)
#
#         lopP_list.append(logP)
#         rho_list.append(rho)
#
#     import matplotlib.pyplot as plt
#     plt.plot(lopP_list)
#     plt.show()
#
#     plt.figure()
#     plt.plot(rho_list)
#
#
#     kallisto_em_iterator(200, ec_vector, ec_2_tr, transcript_length, rho_init=None, reportfreq = 20)
#
#
#
#     rho = np.array([1,1]) / 2
#     transcript_length = np.array([2, 1])
#     ec_vector = np.array([30, 40])
#
#     ec_2_tr = {0: [0],
#                1: [0,1]}
#     for i in range(200):
#
#         rho, a = kallisto_EM_update(rho, ec_vector, ec_2_tr, transcript_length)
#         print(likelihood(rho, ec_vector, ec_2_tr, transcript_length))
#         print(rho)
