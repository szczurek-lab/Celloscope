import numpy as np
import pandas as pd
import scipy.special as sc
from scipy.stats import uniform
from scipy.stats import nbinom
from scipy.stats import gamma
from scipy.stats import norm
from scipy.stats import beta
from scipy.stats import binom
import time
from scipy.stats import poisson
from joblib import Parallel, delayed


import os
import json
import sys




def Celloscope(a_data, a_results, how_many_chains):
    
    def my_logpmf_nb(x, n, p):
        coeff = sc.gammaln(n+x) - sc.gammaln(x+1) - sc.gammaln(n)
        return coeff + n*np.log(p) + sc.xlog1py(x, -p)
    
    def my_trunc_norm_sampling_matrix(mu, sigma):
        n_col = mu.shape[1]
        n_row = mu.shape[0]
        
        U = np.random.mtrand._rand.uniform(size=(n_row,n_col))
        y = mu + sigma*sc.ndtri(U+sc.ndtr(-mu/sigma)*(1-U))
        return y
    
    def my_trunc_norm_sampling_vector(mu, sigma):
        n = len(mu)
        U = np.random.mtrand._rand.uniform(size=n)
        y = mu + sigma*sc.ndtri(U+sc.ndtr(-mu/sigma)*(1-U))
        return y
    
    def my_trunc_norm_sampling_0_1(mu, sigma):
        n = len(mu)
        U = np.random.mtrand._rand.uniform(size=n)
        y = mu + sigma*sc.ndtri(  U*(    sc.ndtr( (1-mu)/sigma  )  - sc.ndtr(-mu/sigma) )   +  sc.ndtr(-mu/sigma) )
        return y
    
    
    def proposal_p_g(curr_p_g, step_size_p_g):
       # nGens = len(curr_p_g)
        prop =  my_trunc_norm_sampling_0_1(curr_p_g, step_size_p_g)
        return prop
    
    def log_lik_p_g(C, h, lambdas, p_g, n_cells):
        nSpots = len(n_cells)
        pg_factor =(1-p_g)/p_g
        p_g_m = np.tile(1-p_g,(nSpots,1) )
        mu_gs = np.matmul(lambdas, np.transpose(h))
        mu_gs = mu_gs * n_cells
        r_gs = ((mu_gs.T)*pg_factor).T
        return np.sum( my_logpmf_nb(C,r_gs,p_g_m.T), axis=1) 
    
    def helper_p_g(p_g, sigma):
        return  norm.cdf( (1 - p_g)/sigma  )  - norm.cdf( -p_g/sigma )
        
    def update_p_g(C, curr_h, curr_lambdas, curr_p_g, curr_n_cells, step_size_p_g):
        nGens = len(curr_p_g)
        prop_p_g = proposal_p_g(curr_p_g ,step_size_p_g) 
        log_curr_lik = log_lik_p_g(C, curr_h, curr_lambdas, curr_p_g, curr_n_cells)
        log_prop_lik = log_lik_p_g(C, curr_h, curr_lambdas, prop_p_g, curr_n_cells) 
        
        log_bias = np.log( helper_p_g(curr_p_g, step_size_p_g) / helper_p_g(prop_p_g, step_size_p_g)  )
        
        r = log_prop_lik - log_curr_lik + log_bias
        los = uniform.rvs(size=nGens)
        decision = r> np.log(los)
        curr_p_g[decision] = prop_p_g[decision]
        
        return (curr_p_g, decision, r)   
    
    def proposal_thetas(thetas_prev, step_size_thetas):
        return(my_trunc_norm_sampling_matrix(thetas_prev, step_size_thetas) )
    
    
    def proposal_lambda_0(curr_lambda_0, step_size_lambda_0):
        return my_trunc_norm_sampling_vector(curr_lambda_0, step_size_lambda_0)
    
    def proposal_over_lambdas(over_lambdas, step_size_over_lambdas):
        return(my_trunc_norm_sampling_vector(over_lambdas, step_size_over_lambdas))
    
    def log_prior_thetas(thetas, Z, aa, b, a_0, b_0):
        prior = gamma.logpdf(thetas, a=aa, scale=b)
        prior_epsilon = gamma.logpdf(thetas, a=a_0, scale=b_0)
        prior[Z==0] = prior_epsilon[Z==0]
        return prior.sum(axis=1)
    
    def log_lik_thetas(C, thetas, lambdas, p_g, n_cells):
        nSpots = len(n_cells)
        row_sums = thetas.sum(axis=1)
        h = thetas/row_sums[:, np.newaxis]
        pg_factor =(1-p_g)/p_g
        p_g_m = np.tile(1-p_g,(nSpots,1) )
        mu_gs = np.matmul(lambdas, np.transpose(h))
        mu_gs = mu_gs * n_cells
        r_gs = ((mu_gs.T)*pg_factor).T
        return my_logpmf_nb(C,r_gs,p_g_m.T).sum(axis=0)
    
    def update_Z(current_thetas, current_pi, a, b, a_0, b_0):
        prob_0 = gamma.logpdf(current_thetas, a_0, scale=b_0) + np.log(1-current_pi)
        prob_1 = gamma.logpdf(current_thetas, a, scale=b) + np.log(current_pi)
        prob = np.exp(prob_0-prob_1)
        p = prob/(1+prob)
        return   binom.rvs(1, 1- p ) 
    
    
    def update_thetas(curr_thetas, curr_Z, curr_lambdas, curr_p_g, curr_n_cells, C, step_size_thetas, a, b, a_0, b_0):
        nSpots = curr_thetas.shape[0]
        prop_thetas = proposal_thetas(curr_thetas.T, step_size_thetas.T).T
    
        log_lik_prop = log_lik_thetas(C, prop_thetas, curr_lambdas, curr_p_g, curr_n_cells)
        log_lik_curr = log_lik_thetas(C, curr_thetas, curr_lambdas, curr_p_g, curr_n_cells)
            
        log_prior_prop = log_prior_thetas(prop_thetas, curr_Z, a, b, a_0, b_0)
        log_prior_curr = log_prior_thetas(curr_thetas, curr_Z, a, b, a_0, b_0)
        
        bias_prop = norm.logcdf(prop_thetas.T, loc=0, scale=step_size_thetas.T).sum(axis=0)
        bias_curr = norm.logcdf(curr_thetas.T, loc=0, scale=step_size_thetas.T).sum(axis=0)
        
    
        r = log_lik_prop -  log_lik_curr + bias_curr - bias_prop  + log_prior_prop - log_prior_curr 
        
        los = uniform.rvs(size=nSpots)
        decision = r> np.log(los)
        curr_thetas[decision] = prop_thetas[decision]
        
        return (curr_thetas, decision, r )
    
    
    def log_lik_lambda_0(C_gs, h, lambda_0, over_lambda, p_g, n_cells):
        nSpots = C_gs.shape[1]
        pg_factor =(1-p_g)/p_g
        lambdas = lambda_0 + over_lambda
        mu_gs = np.matmul(lambdas, np.transpose(h))
        mu_gs = mu_gs * n_cells
        p_g_m = np.tile(1-p_g,(nSpots,1) )
        r_gs = ((mu_gs.T)*pg_factor).T
        return np.sum(nbinom._logpmf(C_gs,r_gs,p_g_m.T))
    
    
    def log_prior_lambda_0(lambda_0, a, b):
        return gamma.logpdf(lambda_0, a=a, scale=b) 
    
    
    def update_lambda_0(C, curr_h, curr_lambda_0, curr_over_lambda, curr_p_g, curr_n_cells, step_size_lambda_0):
        prop_lambda_0 = proposal_lambda_0(curr_lambda_0, step_size_lambda_0)
        
        log_prop_lik = log_lik_lambda_0(C, curr_h, prop_lambda_0, curr_over_lambda, curr_p_g, curr_n_cells)
        log_curr_lik = log_lik_lambda_0(C, curr_h, curr_lambda_0, curr_over_lambda, curr_p_g, curr_n_cells)
        
        bias_prop = norm.logcdf(prop_lambda_0, loc=0, scale=step_size_lambda_0)
        bias_curr = norm.logcdf(curr_lambda_0, loc=0, scale=step_size_lambda_0)
        
        r = log_prop_lik - log_curr_lik + bias_curr - bias_prop  
        los = uniform.rvs(size=1)
    
        decision = r> np.log(los)
        if(decision):
            curr_lambda_0 = prop_lambda_0
        
        return (curr_lambda_0, decision, r) 
    
    def log_lik_n_cells(C, h, lambdas, p_g, n_cells):
        nSpots = len(n_cells)
        pg_factor =(1-p_g)/p_g
        p_g_m = np.tile(1-p_g,(nSpots,1) )
        mu_gs = np.matmul(lambdas, np.transpose(h))
        mu_gs = mu_gs * n_cells
        r_gs = ((mu_gs.T)*pg_factor).T
        return  my_logpmf_nb(C,r_gs,p_g_m.T).sum(axis=0)  
    
    def proposal_n_cells(n_cells, step_size_n_cells):
        return my_trunc_norm_sampling_vector(n_cells, step_size_n_cells)
    
    def est_n_cells_MH(curr_n_cells, C, h, p_g, lambdas, step_size_n_cells, n_cells, mode_number_of_cells, n_cells_prior_strength):
        
        nSpots = len(curr_n_cells)
        prop_n_cells = np.round(proposal_n_cells(curr_n_cells, step_size_n_cells))
        prop_n_cells[prop_n_cells==0] = 1
    
        bias_prop = norm.logcdf(prop_n_cells, loc=0, scale=step_size_n_cells)
        bias_curr = norm.logcdf(curr_n_cells, loc=0, scale=step_size_n_cells)
        
    
        if (mode_number_of_cells=="UNKNOWN"):                  
            log_prior_curr =0
            log_prior_prop = 0
    
        else:
            log_prior_curr = norm.logpdf(curr_n_cells, loc=n_cells, scale= n_cells_prior_strength)
            log_prior_prop = norm.logpdf(prop_n_cells, loc=n_cells, scale= n_cells_prior_strength)
            
        log_curr_lik = log_lik_n_cells(C, h, lambdas, p_g, curr_n_cells)
        log_prop_lik = log_lik_n_cells(C, h, lambdas, p_g, prop_n_cells)
        
        r_n_cells = log_prop_lik - log_curr_lik + bias_curr - bias_prop   + log_prior_prop - log_prior_curr
    
        los = uniform.rvs(size=nSpots)
        decision = r_n_cells> np.log(los)
    
        curr_n_cells[decision] = prop_n_cells[decision]   
        
        return (curr_n_cells, decision, r_n_cells)
    
    
    def log_lik_over_lambdas(C, curr_h, curr_over_lam, curr_lambda_0, curr_p_g, curr_n_cells):
        nSpots = C.shape[1]
        curr_lambdas =  curr_over_lam + curr_lambda_0
        pg_factor =(1-curr_p_g)/curr_p_g
        p_g_m = np.tile(1-curr_p_g,(nSpots,1) )
        mu_gs = np.matmul(curr_lambdas, np.transpose(curr_h))
        mu_gs = mu_gs * curr_n_cells
        r_gs = ((mu_gs.T)*pg_factor).T
        
        return  np.sum(my_logpmf_nb(C,r_gs,p_g_m.T), axis=1)   
    
    
    def update_over_lambdas(C, B, curr_h, curr_over_lam, curr_lambda_0, curr_p_g, curr_n_cells, step_size_over_lambdas, xx0, yy0, xx1, yy1, ind, prior_lambdas, how_many):
        nTypes = B.shape[1]
        nGens = B.shape[0]
        
        vector_curr_over_lambdas = curr_over_lam[ind]
        step_size = np.repeat(step_size_over_lambdas,  how_many)
        vector_proposal_over_lambdas = proposal_over_lambdas(vector_curr_over_lambdas, step_size)
    
        prop_over_lam = np.zeros((nGens,nTypes))
        prop_over_lam[xx1,yy1] = vector_proposal_over_lambdas
        log_curr_lik = log_lik_over_lambdas(C, curr_h, curr_over_lam, curr_lambda_0, curr_p_g, curr_n_cells)
        log_prop_lik = log_lik_over_lambdas(C, curr_h, prop_over_lam, curr_lambda_0, curr_p_g, curr_n_cells)
        
        bias_prop = np.zeros((nGens,nTypes), dtype=float)
        bias_curr = np.zeros((nGens,nTypes), dtype=float)
        
        mat_prior_lambdas = np.zeros((nGens,nTypes), dtype=float)
        mat_prior_lambdas[ind] = np.repeat(prior_lambdas, how_many)
         
        bias_prop[xx1, yy1] = norm.logcdf(vector_proposal_over_lambdas, loc=0, scale=step_size)
        bias_curr[xx1, yy1] = norm.logcdf(vector_curr_over_lambdas, loc=0, scale=step_size)  
        
        bias_prop = bias_prop.sum(axis=1)
        bias_curr = bias_curr.sum(axis=1)
        
        r = log_prop_lik - log_curr_lik   + bias_curr - bias_prop 
        los = uniform.rvs(size=nGens)
        decision = r> np.log(los)
    
        curr_over_lam[decision,] = prop_over_lam[decision,]
        vector_curr_over_lambdas = curr_over_lam[ind]
            
        return (curr_over_lam, decision, vector_curr_over_lambdas, r)
    
    
    with open(a_data + '/params.txt') as f:
        dd = f.read()
        
    params = json.loads(dd) 

    

    
    def Gibbs(chain_nr):
        
        address_results = a_results
        address_data = a_data
        
    
        np.random.seed(int(chain_nr)+1)
    
        burn_in = params['burn in']
        number_of_iterations = params['number of iterations']
        mode_number_of_cells = params['mode number of cells']
        a = params['a']
        b = params['b']
        a_0 = params['a_0']
        b_0 = params['b_0']
        alpha = params['alpha']
    

        thinning_param = params['thinning_parameter']
        n_cells_prior_strength = params['number of cells prior strength']
        how_often_step_size = params['how often update step size']
        how_often_drop = params['how often drop']
    
        
        chain_nr = str("%02d" % chain_nr)
        
        #############################################################################
        ## results addresses  #######################################################
        #############################################################################
        
        address_results = address_results + "/chain" + chain_nr + "/"
    
        
        if not os.path.exists(address_results):
            os.makedirs(address_results)
        
        result_h = pd.DataFrame()
        result_h.to_csv(address_results +"/result_h.csv", header=False, index=False)
    

        result_n_cells = pd.DataFrame()
        result_n_cells.to_csv(address_results +"/result_n_cells.csv", header=False, index=False)
    
        result_p_g = pd.DataFrame()
        result_p_g.to_csv(address_results +"/result_p_g.csv", header=False, index=False)
    
        #############################################################################
        ### loading model input #####################################################
        #############################################################################
    
        C = pd.read_csv(address_data  +   '/C_gs.csv', index_col=0).to_numpy()
        B = pd.read_csv(address_data  +   '/matB.csv', index_col=0).to_numpy().astype(int)
        xx1, yy1 = np.where(B == 1)   
        xx0, yy0 = np.where(B == 0)
        how_many = B.sum(axis=1)
        ind = np.array(B, dtype=bool)
        
        n_cells = pd.read_csv(address_data +   'n_cells.csv')
        n_cells = n_cells['cellCount'].to_numpy().flatten()
    
        prior_lambdas = np.apply_along_axis(lambda v: np.mean(v[v!=0]), 1, C)
        lambdas_start = np.repeat(prior_lambdas, how_many)
        
        nSpots = C.shape[1]
        nTypes = B.shape[1]
        nGens =  C.shape[0]
        
        sum_dec_thetas = np.zeros(nSpots)
        sum_dec_n_cells = np.zeros(nSpots)
        sum_dec_over_lambdas = np.zeros(nGens)
        sum_dec_p_g = np.zeros(nGens)
        sum_dec_lambda_0 = 0
        
        st_thetas = np.repeat(params['step size thetas'], nSpots)
        st_n_cells =  np.repeat(params['step size number of cells'], nSpots)
        st_p_g = np.repeat(params['step size p_g'], nGens)
        st_lambda_0 = params['step size lambda_0']
    
        st_over_lambdas = prior_lambdas/2
        out = np.where(st_over_lambdas<1)
        st_over_lambdas[out] = 1
    
        current_pi = uniform.rvs(loc=0, scale=1, size=(nSpots,nTypes ), random_state=None) 
        current_Z = binom.rvs(1, current_pi)
        current_Z[:,(nTypes-1)] = 0
    
        ############ thetas #####################################################################################################################
        
        current_thetas =  gamma.rvs(np.ones((nSpots,nTypes)))
    
        sum_thetas =  0
    
        current_lambda_0 = 2*uniform.rvs(size=1)  
        current_p_g = uniform.rvs(size=nGens)
        current_n_cells = poisson.rvs(size=nSpots, mu=50) 
        current_over_lambdas = np.zeros((nGens, nTypes), dtype=float)
        current_over_lambdas[xx1,yy1] = lambdas_start
        
        if (mode_number_of_cells=="KNOWN"):
            current_n_cells = n_cells    
            
        assert np.all((B==0) | (B==1)), "elements of B should be 0 or 1"
        assert len(n_cells.shape)==1, "number of cells should be one dimentional vector"
        n_cells_int = n_cells.astype(int)
        assert np.all(n_cells - n_cells_int ==0), "number of cells should be integers"
        C_int = C.astype(int)
        assert np.all(C - C_int ==0), "gene expression counts should be integers"
        assert C.shape[0] == B.shape[0], "number of rows of C should be equal to number of rows of B"
        assert C.shape[0] == B.shape[0], "number of columns of C should be equal length of number of cells vector"
        assert burn_in < number_of_iterations, "burn-in should be smaller than number of iterations"
    
    
        counter = 0
        for i in range(number_of_iterations):
            
            print(i)
    
            current_pi = beta.rvs( alpha/nTypes +current_Z, 2-current_Z)
            current_Z = update_Z(current_thetas, current_pi, a, b, a_0, b_0)
            current_Z[:, nTypes-1] = 0
    
            res_thetas = update_thetas(current_thetas, current_Z, current_over_lambdas + current_lambda_0, current_p_g, current_n_cells, C, st_thetas, a, b, a_0, b_0)
            current_thetas = res_thetas[0]
            sum_dec_thetas = sum_dec_thetas + res_thetas[1]
    
            if ( (i>= burn_in) & (i % thinning_param == 0)   ):
                sum_thetas  = sum_thetas + current_thetas
                counter = counter + 1
    
            ###################### lambda_0     ########################################################################
            
            row_sums = current_thetas.sum(axis=1)
            current_h = current_thetas/row_sums[:, np.newaxis]
            res_lambda_0 = update_lambda_0(C, current_h, current_lambda_0, current_over_lambdas , current_p_g, current_n_cells, st_lambda_0)     
            current_lambda_0 = res_lambda_0[0]
            sum_dec_lambda_0 = sum_dec_lambda_0 + res_lambda_0[1]
    
            current_lambdas = current_over_lambdas + current_lambda_0
            
            ######################  p_g    #############################################################################
            
            res_p_g = update_p_g(C, current_h, current_lambdas, current_p_g, current_n_cells, st_p_g)
            current_p_g =  res_p_g[0]
            sum_dec_p_g =  sum_dec_p_g + res_p_g[1]
    
            ######################  n_cells   ##########################################################################
            
            if (mode_number_of_cells=="UNKNOWN" or mode_number_of_cells=="ASPRIORS"):
                res_n_cells = est_n_cells_MH(current_n_cells, C, current_h, current_p_g, current_lambdas, st_n_cells, n_cells, mode_number_of_cells, n_cells_prior_strength)
                current_n_cells = res_n_cells[0]
                sum_dec_n_cells = sum_dec_n_cells + res_n_cells[1]
    
            #############  over_lambdas ################################################################################
        
            res_over_lambdas = update_over_lambdas(C, B, current_h, current_over_lambdas, current_lambda_0, current_p_g, current_n_cells, st_over_lambdas, xx0, yy0, xx1, yy1, ind, prior_lambdas, how_many)
            current_over_lambdas = res_over_lambdas[0]
            sum_dec_over_lambdas = sum_dec_over_lambdas + res_over_lambdas[1]
    
    
            if( (i % how_often_step_size == 0) & (i>0) & (i<burn_in) ):
                
                thetas_acceptance_ratio = sum_dec_thetas/how_often_step_size
                lambda_0_acceptance_ratio = sum_dec_lambda_0/how_often_step_size
                over_lambdas_acceptance_ratio = sum_dec_over_lambdas/how_often_step_size
                p_g_acceptance_ratio = sum_dec_p_g/how_often_step_size
    
                sum_dec_thetas = np.zeros(nSpots)
                sum_dec_lambda_0 = 0
                sum_dec_over_lambdas = np.zeros(nGens)
                sum_dec_p_g = np.zeros(nGens)
                
                where_small_accpt = thetas_acceptance_ratio <= 0.23
                
                st_thetas[where_small_accpt] = st_thetas[where_small_accpt]*0.95
                st_thetas[~where_small_accpt] = st_thetas[~where_small_accpt]*1.05
    
                if (lambda_0_acceptance_ratio <= 0.23):
                    st_lambda_0 =  st_lambda_0*0.99
                else:
                    st_lambda_0 =  st_lambda_0*1.01
                    
                where_small_accpt_lambdas = over_lambdas_acceptance_ratio <= 0.23
                
                st_over_lambdas[where_small_accpt_lambdas] = st_over_lambdas[where_small_accpt_lambdas]*0.95
                st_over_lambdas[~where_small_accpt_lambdas] = st_over_lambdas[~where_small_accpt_lambdas]*1.05
                
                where_small_accpt_p_g  = p_g_acceptance_ratio <= 0.23
                
                st_p_g[where_small_accpt_p_g] = st_p_g[where_small_accpt_p_g]*0.95
                st_p_g[~where_small_accpt_p_g] = st_p_g[~where_small_accpt_p_g]*1.05
                
                
            if(i % how_often_drop == 0):
                pd.DataFrame((current_h).reshape(nSpots*nTypes)).T.to_csv( address_results + "/result_h.csv", header=False, index=False, mode="a")
                pd.DataFrame(current_n_cells).T.to_csv(address_results + "/result_n_cells.csv", header=False, index=False, mode="a")
                pd.DataFrame(current_p_g).T.to_csv(address_results + "/result_p_g.csv", header=False, index=False, mode="a")
    
    
        theta_est = sum_thetas/counter
        pd.DataFrame(theta_est).T.to_csv(address_results  +"/thetas_est.csv", header=False, index=False)
    
    how_many_chains = int(how_many_chains)
    
    start_time = time.perf_counter()
    result = Parallel(n_jobs=how_many_chains)(delayed(Gibbs)(i) for i in range(1,  how_many_chains+1   ))
    finish_time = time.perf_counter()
    print(f"Program finished in {finish_time-start_time} seconds")
    print(result)



ad_data = sys.argv[1]
ad_results = sys.argv[2]
how_many_chains = sys.argv[3]
Celloscope(ad_data, ad_results, how_many_chains)

