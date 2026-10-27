# -*- coding: utf-8 -*-
# pylint: disable=C0103
"""
Estimate cost of solving LWE using primal attacks.

See :ref:`LWE Primal Attacks` for an introduction what is available.

"""
from functools import partial
from sage.all import oo, cached_function, ceil, exp, floor, log, RR, sqrt

from .conf import red_cost_model as red_cost_model_default
from .conf import red_shape_model as red_shape_model_default
from .cost import Cost
from .io import Logging
from .lwe_comb import log_comb, split_weight, sum_log
from .lwe_parameters import LWEParameters
from .lwe_primal import primal_usvp, PrimalUSVP, PrimalHybrid
from .nd import SparseTernary
from .prob import babai_gaussian, mitm_babai_probability, mitm_babai_probability_single, amplify as prob_amplify
from .reduction import cost as costf, delta as deltaf
from .simulator import normalize as simulator_normalize
from .util import local_minimum


class PrimalMeet:
    """
    Estimate cost of solving LWE via a hybrid of primal attack and Meet-LWE [HKLS22].
    """

    @staticmethod
    def cost_meet_lwe(search_space: SparseTernary, r, sigma: float, prob_adm: float):
        """
        Compute the cost of doing REP-1 depth 2 Meet-LWE for a given search space,

        :param search_space: search space to perform REP-1 over.
        :param r: simulated profile of reduced basis over which to meet.
        :param sigma: stddev of error (noise).
        :param prob_adm: admissability probability of `r`.

        """
        n, w, is_balanced = search_space.n, search_space.hamming_weight, search_space.ones is not None
        n1, n2 = split_weight(n)
        w0 = w // 2  # The secret has `w` nonzero coefficients

        def cost_epsilon(epsilon):
            logS = {}  # log(size of search space)
            logL = {}  # log(size of list of stored candidates)
            if is_balanced:
                # If balanced, the secret has `w0` coefficients equal to 1, and `w0` equal to -1.
                w1, w2 = split_weight(w0 + 2 * epsilon)  # Adds epsilon to each of w1 and w2.
                w11, w12 = split_weight(w1)
                w21, w22 = split_weight(w2)

                if n - w < 2 * epsilon or w0 < w2 - epsilon:
                    return None

                # Number of ways that we can construct s from s1 + s2:
                # s has n-w zero entries:
                # `epsilon` come from 1 - 1, `epsilon` from -1 + 1, and rest from 0 + 0.
                # Half of the w0 +1's in s come from s1, and the other half from s2.
                # Similarly for the -1's.
                log_R1 = log_comb(n - w, epsilon, epsilon) + log_comb(w0, w1 - epsilon) + log_comb(w0, w2 - epsilon)

                # Analyse level 2.
                logL[11] = logS[11] = log_comb(n1, w11, w11)
                logL[12] = logS[12] = log_comb(n2, w12, w12)
                logL[21] = logS[21] = log_comb(n1, w21, w21)
                logL[22] = logS[22] = log_comb(n2, w22, w22)

                # Analyse level 1.
                logS[1] = log_comb(n, w1, w1)
                logS[2] = log_comb(n, w2, w2)
            else:
                # If unbalanced, the secret has just `w` nonzero coefficients that could be `1` or `-1`.
                w1 = w0 + epsilon  # weight at level 1
                w11, w12 = split_weight(w1)  # weights at level 2: left and right

                if w + epsilon > n:
                    return None

                # Number of ways that we can construct s from s1 + s2:
                # s has n-w zero entries, and `epsilon` are formed by 0 = 1 + (-1) or 0 = (-1) + 1.
                log_R1 = log_comb(n - w, epsilon) + log_comb(w, w0) + epsilon * RR(log(2))

                # Analyse level 2.
                logL[11] = logS[11] = logL[21] = logS[21] = log_comb(n1, w11) + w11 * RR(log(2))
                logL[12] = logS[12] = logL[22] = logS[22] = log_comb(n2, w12) + w12 * RR(log(2))

                # Analyse level 1.
                logS[1] = logS[2] = log_comb(n, w1) + w1 * RR(log(2))

            logL[1] = logS[1] - log_R1
            logL[2] = logS[2] - log_R1

            # Analyse the runtime
            log_runtime = RR(sum_log(*logL.values()))

            # Analyse the success probability
            bet_s1 = logS[11] + logS[12] - logS[1]
            bet_s2 = logS[21] + logS[22] - logS[2]
            assert bet_s1 <= 0 and bet_s2 <= 0
            log_bet = bet_s1 + bet_s2

            # Logarithm of the number of buckets needed to have an efficient meet algorithm.
            # When using less, some buckets might get overcrowded, so listing near-collisions becomes slow.
            log_meet_leftover = RR(max(logS[1], logS[2]))
            my_prob_adm, dim_babai = prob_adm, 0
            bound_meet = RR(12 * sigma)**2
            for i in range(len(r))[::-1]:
                dim_babai += 1
                if r[i] >= 4 * bound_meet:
                    c_i = floor((r[i] / bound_meet)**.5)  # take square root, because `r` contains square norms.
                    # Compute p_adm(r') using p_adm(r), and by recomputing only the modified part.
                    my_prob_adm /= mitm_babai_probability_single(r[i], sigma)
                    my_prob_adm *= mitm_babai_probability_single(r[i] / c_i**2, sigma)
                    log_meet_leftover -= RR(log(c_i))
                    if log_meet_leftover < 0:
                        break

            # Multiply the runtime by the number of calls to Babai in dimension `d`.
            log_runtime += RR(log(PrimalHybrid.babai_cost(dim_babai)["rop"]))

            # Repeat the MEET algorithm 1 / (exp(log_bet) my_prob_adm) times.
            # All these repetitions fail with probability of 1/e.
            runtime, prob = RR(exp(log_runtime - log_bet) / my_prob_adm), RR(1 - exp(-1))
            return Cost(rop=runtime, mem=runtime, prob=prob, h_1=w1, h_2=w11, epsilon=epsilon, d_=dim_babai)

        # Find the best epsilon by exhaustive search.
        # In practice, `epsilon` is 1, 2 or 3 so keep the range limited for performance.
        best_cost = Cost(rop=oo, mem=oo, prob=1.0)
        max_epsilon = 4

        # Compute costs for all `epsilon <= max_epsilon`, and pick the best.
        for epsilon in range(1, max_epsilon + 1):
            cost = cost_epsilon(epsilon)
            if cost and cost["rop"] < best_cost["rop"]:
                best_cost = cost  # pick this one with smaller runtime

        # Try to improve the epsilon until a worse one is found
        while best_cost["rop"] == oo or best_cost["epsilon"] == max_epsilon:
            max_epsilon += 1
            cost = cost_epsilon(max_epsilon)
            if not cost:
                break  # weight is too high, stop.
            if cost["rop"] < best_cost["rop"]:
                best_cost = cost  # pick this one with smaller runtime
        return best_cost

    @staticmethod
    @cached_function
    def cost(
        params: LWEParameters,
        beta: int,
        zeta: int,
        succ_prob: float = 0.99,
        red_shape_model=red_shape_model_default,
        red_cost_model=red_cost_model_default,
        log_level=5,
    ):
        """
        Cost of the Primal Meet-LWE attack.

        :param beta: Block size.
        :param params: LWE parameters.
        :param zeta: Guessing dimension ζ ≥ 0.
        :param m: We accept the number of samples to consider from the calling function.

        .. note :: This is the lowest level function that runs no optimization. It merely reports costs.

        """
        assert (
            isinstance(params.Xs, SparseTernary) and
            (params.Xs.ones is None or 2 * params.Xs.ones == params.Xs.hamming_weight)  # [un|]balanced
        )

        delta = deltaf(beta)
        n = params.n - zeta
        d = min(max(beta, ceil(sqrt(n * log(params.q) / log(delta)))), params.m + n)
        xi = PrimalUSVP._xi_factor(params.Xs, params.Xe)

        # 1. Simulate BKZ-β
        sigma = params.Xe.stddev
        simulator = simulator_normalize(red_shape_model)
        r = simulator(d, n, params.q, beta, xi=xi, tau=None, dual=True)

        cost_bkz = RR(costf(red_cost_model, beta, d)["rop"])

        prob_np = RR(babai_gaussian(r, sigma))
        prob_adm = RR(mitm_babai_probability(r, sigma, fast=1000))

        def cost_hw(hw):
            search_space = params.Xs.split_balanced(zeta, hw)[0]
            prob_hw = params.Xs.split_probability(zeta, hw)

            # 3. Determine cost of doing Meet LWE.
            cost_meet = PrimalMeet.cost_meet_lwe(search_space, r, sigma, prob_adm)
            if cost_meet["rop"] == oo:
                return cost_meet

            probability = (
                prob_hw  # prob. secret splits with given weights. Note: p_HW ~ 1/n roughly
                * prob_np  # prob. correct guess lifts with np. Note: p_NP ~ 1.0 (p_NP > p_adm).
                # * prob_adm  # prob. s = s1 + s2 has both (s1, s2) in the same bucket. (now part of MEET-LWE)
                * cost_meet["prob"]  # prob. Meet-LWE gives the correct answer.
            )

            cost = Cost({
                "rop": cost_bkz + cost_meet["rop"], "red": cost_bkz,
                "mem": cost_meet["mem"],
                "beta": beta, "zeta": zeta, "d": d,
                "h_": hw, "h_1": cost_meet["h_1"], "h_2": cost_meet["h_2"],
                "epsilon": cost_meet["epsilon"], "d_": cost_meet["d_"],
                "|S|": search_space, "prob": probability,
            })

            # 4. Repeat whole experiment ~1/prob times
            # assert not (not probability or RR(probability).is_NaN())
            return cost.repeat(prob_amplify(succ_prob, probability))

        # First, find a HW such that prob_hw(HW) >= 2^{-10}.
        hw, max_hw = 0, min(params.Xs.hamming_weight, zeta)
        while hw < max_hw and params.Xs.split_probability(zeta, hw) < 2**-10:
            hw += 2

        best_cost, orig_hw = cost_hw(hw), hw
        # Then, decrease the HW as long as cost decreases.
        for hw_ in range(hw - 2, 0, -2):
            cost = cost_hw(hw_)
            if cost > best_cost:
                break
            best_cost = cost
            hw = hw_
        if hw == orig_hw:
            # Or, increase the HW as long as cost decreases.
            for hw_ in range(hw + 2, max_hw, 2):
                cost = cost_hw(hw_)
                if cost > best_cost:
                    break
                best_cost = cost
                hw = hw_
        return best_cost

    @classmethod
    def cost_zeta(
        cls,
        zeta: int,
        params: LWEParameters,
        succ_prob: float = 0.99,
        red_shape_model=red_shape_model_default,
        red_cost_model=red_cost_model_default,
        log_level=5,
        **kwds,
    ):
        """
        This function optimizes costs for a fixed guessing dimension ζ.
        """

        # step 0. establish baseline
        baseline_cost = primal_usvp(
            params,
            red_shape_model=red_shape_model,
            red_cost_model=red_cost_model,
            optimize_d=False,
            log_level=log_level + 1,
            **kwds,
        )
        Logging.log("bdd", log_level, f"H0: {repr(baseline_cost)}")

        f = partial(
            PrimalMeet.cost,
            params,
            zeta=zeta,
            succ_prob=succ_prob,
            red_shape_model=red_shape_model,
            red_cost_model=red_cost_model,
            **kwds,
        )

        # step 1. optimize β
        with local_minimum(
            40, baseline_cost["beta"] + 1, precision=2, log_level=log_level + 1
        ) as it:
            for beta in it:
                it.update(f(beta))
            cost = it.y

        Logging.log("bdd", log_level, f"H1: {cost!r}")

        if cost is None:
            return Cost(rop=oo)
        return cost

    @classmethod
    def __call__(
        cls,
        params: LWEParameters,
        zeta: int = None,
        succ_prob: float = 0.99,
        red_shape_model=red_shape_model_default,
        red_cost_model=red_cost_model_default,
        log_level=1,
        **kwds,
    ):
        """
        Estimate cost of solving LWE via a hybrid of primal attack and Meet-LWE [HKLS22].
        """
        params = LWEParameters.normalize(params)
        Cost.register_impermanent(h_1=False, h_2=False, d_=False, epsilon=False)

        f = partial(
            cls.cost_zeta,
            params=params,
            succ_prob=succ_prob,
            red_shape_model=red_shape_model,
            red_cost_model=red_cost_model,
            log_level=log_level + 1,
        )

        if zeta is None:
            cost = Cost(rop=oo)
            with local_minimum(0, params.n, precision=2, log_level=log_level) as it:
                for zeta_ in it:
                    it.update(f(zeta=zeta_, **kwds))
                cost = it.y
        else:
            cost = f(zeta=zeta)

        cost["tag"] = "primal-meet-hybrid"
        cost["problem"] = params
        return cost.sanity_check()

    __name__ = "primal_meet"


primal_meet = PrimalMeet()
