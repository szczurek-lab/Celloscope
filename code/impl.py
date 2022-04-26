import numpy as np
import pandas as pd
import scipy.special as sc
from scipy.stats import uniform
from scipy.stats import nbinom
from scipy.stats import gamma
from scipy.stats import norm
from scipy.stats import beta
from scipy.stats import binom
from scipy.stats import poisson
from multiprocess import Pool

import json
import os
import sys


def calosc(adres_data, adres_results, how_many_chains):
    
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
    
    def proposal_p_g(curr_p_g, step_size_p_g):
        nGens = len(curr_p_g)
        prop =  np.random.normal(loc=curr_p_g, scale=step_size_p_g, size=nGens)
        return(prop%1)
    
    def log_lik_p_g(C, h, lambdas, p_g, n_cells):
        nSpots = len(n_cells)
        pg_factor =(1-p_g)/p_g
        p_g_m = np.tile(1-p_g,(nSpots,1) )
        mu_gs = np.matmul(lambdas, np.transpose(h))
        mu_gs = mu_gs * n_cells
        r_gs = ((mu_gs.T)*pg_factor).T
        return np.sum( my_logpmf_nb(C,r_gs,p_g_m.T), axis=1) 
        
    def update_p_g(C, curr_h, curr_lambdas, curr_p_g, curr_n_cells, step_size_p_g):
        nGens = len(curr_p_g)
        prop_p_g = proposal_p_g(curr_p_g ,step_size_p_g) 
        log_curr_lik = log_lik_p_g(C, curr_h, curr_lambdas, curr_p_g, curr_n_cells)
        log_prop_lik = log_lik_p_g(C, curr_h, curr_lambdas, prop_p_g, curr_n_cells)  
        r = log_prop_lik - log_curr_lik
        los = uniform.rvs(size=nGens)
        decision = r> np.log(los)
        curr_p_g[decision] = prop_p_g[decision]
        
        return (curr_p_g, decision)  
    
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
        
        
        return (curr_thetas, decision )
    
    
    
    
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
        
        r = log_prop_lik - log_curr_lik + bias_curr - bias_prop  + norm.logpdf(prop_lambda_0, loc=1, scale=5) - norm.logpdf(curr_lambda_0, loc=1, scale=5)
    
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
    
    def tNCDF(x, mu, sigma):
        return 1 - ( (1 - norm.cdf(  (x - mu )/sigma  ))/( norm.cdf(mu/sigma)))
    
    def density_ceil_tNorm(x, mu, sigma):
        return tNCDF(x, mu, sigma) - tNCDF(x-1, mu, sigma)
    
    def est_n_cells_MH(curr_n_cells, C, h, p_g, lambdas, step_size_n_cells, n_cells, n_cells_prior_strength):
        
        nSpots = len(curr_n_cells)
        prop_n_cells = np.ceil(proposal_n_cells(curr_n_cells, step_size_n_cells))
        bias_prop = np.log(  density_ceil_tNorm( prop_n_cells, curr_n_cells, step_size_n_cells ))
        bias_curr = np.log(   density_ceil_tNorm( curr_n_cells, prop_n_cells,   step_size_n_cells))
        
        log_prior_curr = norm.logpdf(curr_n_cells, loc=n_cells, scale= n_cells_prior_strength)
        log_prior_prop = norm.logpdf(prop_n_cells, loc=n_cells, scale= n_cells_prior_strength)
            
        log_curr_lik = log_lik_n_cells(C, h, lambdas, p_g, curr_n_cells)
        log_prop_lik = log_lik_n_cells(C, h, lambdas, p_g, prop_n_cells)
        
        r_n_cells = log_prop_lik - log_curr_lik + bias_curr - bias_prop   + log_prior_prop - log_prior_curr
    
        los = uniform.rvs(size=nSpots)
        decision = r_n_cells> np.log(los)
    
        curr_n_cells[decision] = prop_n_cells[decision]   
        
        return (curr_n_cells, decision)
    
    
    def log_lik_over_lambdas(C, curr_h, curr_over_lam, curr_lambda_0, curr_p_g, curr_n_cells):
        nSpots = C.shape[1]
        curr_lambdas =  curr_over_lam + curr_lambda_0
        pg_factor =(1-curr_p_g)/curr_p_g
        p_g_m = np.tile(1-curr_p_g,(nSpots,1) )
        mu_gs = np.matmul(curr_lambdas, np.transpose(curr_h))
        mu_gs = mu_gs * curr_n_cells
        r_gs = ((mu_gs.T)*pg_factor).T
        
        return  np.sum(my_logpmf_nb(C,r_gs,p_g_m.T), axis=1)   
    
    def update_over_lambdas(C, B, curr_h, curr_over_lam, curr_lambda_0, curr_p_g, curr_n_cells, step_size_over_lambdas, xx0, yy0, xx1, yy1, ind, ile_po_kolei):
        nTypes = B.shape[1]
        nGens = B.shape[0]
         
        vector_curr_over_lambdas = curr_over_lam[ind]
    
        step_size = np.repeat(step_size_over_lambdas,  ile_po_kolei)
        vector_proposal_over_lambdas = proposal_over_lambdas(vector_curr_over_lambdas, step_size)
    
        prop_over_lam = np.zeros((nGens,nTypes))
        prop_over_lam[xx1,yy1] = vector_proposal_over_lambdas
        log_curr_lik = log_lik_over_lambdas(C, curr_h, curr_over_lam, curr_lambda_0, curr_p_g, curr_n_cells)
        log_prop_lik = log_lik_over_lambdas(C, curr_h, prop_over_lam, curr_lambda_0, curr_p_g, curr_n_cells)
        
        bias_prop = np.zeros((nGens,nTypes), dtype=float)
        bias_curr = np.zeros((nGens,nTypes), dtype=float)
         
        bias_prop[xx1, yy1] = norm.logcdf(vector_proposal_over_lambdas, loc=0, scale=step_size)
        bias_curr[xx1, yy1] = norm.logcdf(vector_curr_over_lambdas, loc=0, scale=step_size)  
        
        bias_prop = bias_prop.sum(axis=1)
        bias_curr = bias_curr.sum(axis=1)
        
        r = log_prop_lik - log_curr_lik   + bias_curr - bias_prop 
        los = uniform.rvs(size=nGens)
        decision = r> np.log(los)
    
        curr_over_lam[decision,] = prop_over_lam[decision,]
        vector_curr_over_lambdas = curr_over_lam[ind]
            
        return (curr_over_lam, decision, vector_curr_over_lambdas) 
    
    
    
    def Gibbs(set_up_number, chain_nr, params):
        
   
        address_data = adres_data
        address_results = adres_results
    
        np.random.seed(int(chain_nr)+1)
    
        burn_in = params['burn in']
        number_of_iterations = params['number of iterations']
        mode_number_of_cells = params['mode number of cells']
        a = params['a']
        b = params['b']
        a_0 = params['a_0']
        b_0 = params['b_0']
        alpha = params['alpha']
        nr_of_spots_to_follow = params['nr_of_spots_to_follow']
        thinning_param = params['thinning_parameter']
        n_cells_prior_strength = params['number of cells prior strength']
        how_often_step_size = params['how often update step size']
        how_often_drop = params['how often drop']
        
        
       
        set_up_number  = str("%02d" %  set_up_number)
        chain_nr = str("%02d" % chain_nr)
        
        #############################################################################
        ## results addresses  #######################################################
        #############################################################################
        
        adress_res_base = address_results + "/chain" + chain_nr + "/"
        #address_res_Z =  adress_res_base +"result_Z.csv"
        #address_res_lambda_0 = adress_res_base + "result_lambda_0.csv"
        #address_res_n_cells = adress_res_base   +"result_n_cells.csv"
        address_res_p_g = adress_res_base   +"result_p_g.csv"
    
        
        #############################################################################
        ### loading model input #####################################################
        #############################################################################
        
        C = pd.read_csv(address_data  +   '/C_gs.csv', index_col=0).to_numpy()
        B = pd.read_csv(address_data  +   '/matB.csv', index_col=0).to_numpy()
        xx1, yy1 = np.where(B == 1)   
        xx0, yy0 = np.where(B == 0)
        how_many_each_type = B.sum(axis=1)
        ind = np.array(B, dtype=bool)
        
        
        n_cells = pd.read_csv(address_data +   '/n_cells.csv')
        n_cells = n_cells['no.of.nuclei'].to_numpy().flatten()
    
        nSpots = C.shape[1]
        nTypes = B.shape[1]
        nGens =  C.shape[0]
        
        #nr_of_spots_to_follow
        #number_of_entries = nTypes * nSpots
        #number_of_entries = nTypes * 200
        number_of_entries = nTypes * nr_of_spots_to_follow
        
        sum_dec_thetas = np.zeros(nSpots)
        sum_dec_n_cells = np.zeros(nSpots)
        sum_dec_over_lambdas = np.zeros(nGens)
        sum_dec_p_g = np.zeros(nGens)
        sum_dec_lambda_0 = 0
        
        st_thetas = np.repeat(params['step size thetas'], nSpots)
        st_n_cells =  np.repeat(params['step size number of cells'], nSpots)
        st_over_lambdas = np.repeat( params['step size over lambdas'], nGens)
        
        st_p_g = np.repeat(params['step size p_g'], nGens)
        st_lambda_0 = params['step size lambda_0']
        
    
        if not os.path.exists(address_results +  "//chain" + chain_nr ):
            os.makedirs(address_results +  "//chain" + chain_nr)
        
        thetas_acceptance_ratio = pd.DataFrame()
        #thetas_acceptance_ratio.to_csv(adress_res_base +"thetas_acceptance_ratio.csv", header=False, index=False)
        
        lambda_0_acceptance_ratio = pd.DataFrame()
        #lambda_0_acceptance_ratio.T.to_csv(adress_res_base +"lambda_0_acceptance_ratio.csv", header=False, index=False)
        
        over_lambdas_acceptance_ratio = pd.DataFrame()
        #over_lambdas_acceptance_ratio.T.to_csv(adress_res_base +"over_lambdas_acceptance_ratio.csv", header=False, index=False)
        
        p_g_acceptance_ratio = pd.DataFrame()
        #p_g_acceptance_ratio.T.to_csv(adress_res_base +"p_g_acceptance_ratio.csv", header=False, index=False)
        
        #n_cells_acceptance_ratio = pd.DataFrame()
        #n_cells_acceptance_ratio.T.to_csv(adress_res_base +"n_cells_acceptance_ratio.csv", header=False, index=False)
    
        #pd.DataFrame(st_thetas).T.to_csv(adress_res_base +"st_thetas.csv", header=False, index=False)
        #pd.DataFrame(st_p_g).T.to_csv(adress_res_base +"st_p_g.csv", header=False, index=False)
        #pd.DataFrame(st_over_lambdas).T.to_csv(adress_res_base +"st_over_lambdas.csv", header=False, index=False)
        
        current_pi = uniform.rvs(loc=0, scale=1, size=(nSpots,nTypes ), random_state=None)
        current_Z = binom.rvs(1, current_pi)
        #result_Z = pd.DataFrame(current_Z)
        #result_Z.T.to_csv(address_res_Z, header=False, index=False)
        
        
        ############ thetas #####################################################################################################################
        
        current_thetas =  gamma.rvs(np.ones((nSpots,nTypes)))
        
       
        suma_thetas =  0
        #suma_Z = 0
    
        
        h_output = pd.DataFrame()
        h_output.to_csv(adress_res_base +"h_output.csv", header=False, index=False)
        
        over_lambdas_output = pd.DataFrame()
        over_lambdas_output.to_csv(adress_res_base +"over_lambdas_output.csv", header=False, index=False)
        
        n_cells_output = pd.DataFrame()
        n_cells_output.to_csv(adress_res_base +"n_cells_output.csv", header=False, index=False)
        
    
       
        ############# lambda_0 ######################################################################################################################
        
        current_lambda_0 = 2*uniform.rvs(size=1)  
        #pd.DataFrame().to_csv(address_res_lambda_0 , header=False, index=False)
        
        #############################################################################################################################################
        
        ############# p_g ###########################################################################################################################
        
        current_p_g = uniform.rvs(size=nGens)
        #pd.DataFrame().to_csv(address_res_p_g, header=False, index=False)
        
        #############################################################################################################################################
        
        ############# n_cells  ######################################################################################################################
        
        current_n_cells = poisson.rvs(size=nSpots, mu=50) 
        
        if (mode_number_of_cells=="KNOWN"):
            current_n_cells = n_cells     
        
        #############################################################################################################################################
        
        #############  over_lambdas #################################################################################################################
        
        current_over_lambdas = np.zeros((nGens, nTypes), dtype=float)
        current_over_lambdas[xx1,yy1] = 10*uniform.rvs(size=np.sum(B))
        
        ##############################################
    
    
        counter = 0
        for i in range(number_of_iterations):
            
            print(i)
    
            current_pi = beta.rvs( alpha/nTypes +current_Z, 2-current_Z)
            current_Z = update_Z(current_thetas, current_pi, a, b, a_0, b_0)
    
    
            res_thetas = update_thetas(current_thetas, current_Z, current_over_lambdas + current_lambda_0, current_p_g, current_n_cells, C, st_thetas, a, b, a_0, b_0)
            current_thetas = res_thetas[0]
            sum_dec_thetas = sum_dec_thetas + res_thetas[1]
    
    
            if ( (i>= burn_in) & (i % thinning_param == 0)   ):
                suma_thetas  = suma_thetas + current_thetas
                #suma_Z = suma_Z + current_Z
                counter = counter + 1
    
    
            ###################### lambda_0     ########################################################################
            
            row_sums = current_thetas.sum(axis=1)
            curr_h = current_thetas/row_sums[:, np.newaxis]
            res_lambda_0 = update_lambda_0(C, curr_h, current_lambda_0, current_over_lambdas , current_p_g, current_n_cells, st_lambda_0)     
            current_lambda_0 = res_lambda_0[0]
            
            sum_dec_lambda_0 = sum_dec_lambda_0 + res_lambda_0[1]
            
            ############################################################################################################
            
            current_lambdas = current_over_lambdas + current_lambda_0
            
            ######################  p_g    #############################################################################
            
            res_p_g = update_p_g(C, curr_h, current_lambdas, current_p_g, current_n_cells, st_p_g)
            current_p_g =  res_p_g[0]
            sum_dec_p_g =  sum_dec_p_g + res_p_g[1]
                 
            ############################################################################################################
            
            ######################  n_cells   ##########################################################################
            
    
            if (mode_number_of_cells=="ASPRIORS"):
                res_n_cells = est_n_cells_MH(current_n_cells, C, curr_h, current_p_g, current_lambdas, st_n_cells, n_cells, n_cells_prior_strength)
                current_n_cells = res_n_cells[0]
                sum_dec_n_cells = sum_dec_n_cells + res_n_cells[1]
    
    
                if(i % thinning_param == 0):
                    n_cells_output = n_cells_output.append(pd.DataFrame(current_n_cells).T)
                    
                if(i % how_often_drop  == 0):
                    n_cells_output.to_csv(adress_res_base +"n_cells_output.csv", header=False, index=False, mode="a")
                    n_cells_output = pd.DataFrame()
                    
                    
                    
            ############################################################################################################
            
            #############  over_lambdas ################################################################################
            
    
            
            res_over_lambdas = update_over_lambdas(C, B, curr_h, current_over_lambdas, current_lambda_0, current_p_g, current_n_cells, st_over_lambdas, xx0, yy0, xx1, yy1, ind, how_many_each_type)
            current_over_lambdas = res_over_lambdas[0]
            current_over_lambdas_vector = pd.DataFrame(res_over_lambdas[2]).T
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
                
                pd.DataFrame(st_thetas).T.to_csv(adress_res_base +"st_thetas.csv", header=False, index=False, mode="a")
                
                if (lambda_0_acceptance_ratio <= 0.23):
                    st_lambda_0 =  st_lambda_0*0.99
                else:
                    st_lambda_0 =  st_lambda_0*1.01
                    
                where_small_accpt_lambdas = over_lambdas_acceptance_ratio <= 0.23
                
                st_over_lambdas[where_small_accpt_lambdas] = st_over_lambdas[where_small_accpt_lambdas]*0.95
                st_over_lambdas[~where_small_accpt_lambdas] = st_over_lambdas[~where_small_accpt_lambdas]*1.05
                
                pd.DataFrame(st_over_lambdas).T.to_csv(adress_res_base +"st_over_lambdas.csv", header=False, index=False, mode="a")
                
                where_small_accpt_p_g  = p_g_acceptance_ratio <= 0.23
                
                st_p_g[where_small_accpt_p_g] = st_p_g[where_small_accpt_p_g]*0.95
                st_p_g[~where_small_accpt_p_g] = st_p_g[~where_small_accpt_p_g]*1.05
                
                #pd.DataFrame(st_p_g).T.to_csv(adress_res_base +"st_p_g.csv", header=False, index=False, mode="a")
                
                #pd.DataFrame(thetas_acceptance_ratio).T.to_csv(adress_res_base +"thetas_acceptance_ratio.csv", header=False, index=False, mode="a")
                #pd.DataFrame(lambda_0_acceptance_ratio).T.to_csv(adress_res_base +"lambda_0_acceptance_ratio.csv", header=False, index=False, mode="a")
                #pd.DataFrame(over_lambdas_acceptance_ratio).T.to_csv(adress_res_base +"over_lambdas_acceptance_ratio.csv", header=False, index=False, mode="a")
                #pd.DataFrame(p_g_acceptance_ratio).T.to_csv(adress_res_base +"p_g_acceptance_ratio.csv", header=False, index=False, mode="a")
                
    
            if(i % thinning_param == 0):
                vector_h_curr = pd.DataFrame(curr_h[0:nr_of_spots_to_follow,:].reshape(number_of_entries)).T
                h_output = h_output.append(vector_h_curr)
                over_lambdas_output = over_lambdas_output.append(current_over_lambdas_vector)
                
                result_p_g = pd.DataFrame(current_p_g)
                result_p_g.T.to_csv(address_res_p_g, header=False, index=False, mode="a")
                
                
            
            if(i % how_often_drop  == 0):
                h_output.to_csv(adress_res_base +"h_output.csv", header=False, index=False, mode="a")
                h_output = pd.DataFrame()
                
                over_lambdas_output.to_csv(adress_res_base +"over_lambdas_output.csv", header=False, index=False, mode="a")
                over_lambdas_output = pd.DataFrame()
    
    
        theta_est = suma_thetas/counter
        row_sums = theta_est.sum(axis=1)
        h_est = theta_est/row_sums[:, np.newaxis]
        pd.DataFrame(theta_est).T.to_csv(adress_res_base  +"thetas_est.csv", header=False, index=False)
        pd.DataFrame(h_est).T.to_csv(adress_res_base  +"h_est.csv", header=False, index=False)
        #pd.DataFrame(Z_est).T.to_csv(adress_res_base +"Z_est.csv", header=False, index=False)
    
    
                                  
    

    with open(adres_data + '/params.txt') as f:
        dd = f.read()
        
    my_params = json.loads(dd) 

    
    def multi_setup(chain_nr):
        Gibbs(1, chain_nr,  my_params)
        
    how_many_chains = int(how_many_chains)    
    
    if __name__ == '__main__':
        with Pool(how_many_chains) as p: 
            print(p.map(multi_setup, [i+1 for i in range(how_many_chains)]))



ad_data = sys.argv[1]
ad_results = sys.argv[2]
how_many_chains = sys.argv[3]
calosc(ad_data, ad_results, how_many_chains)
