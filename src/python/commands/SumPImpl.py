import os,sys,math,random
from phycas import *
from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions

class VarianceZeroError(Exception):
    def __init__(self):
        self.msg = 'Cannot calculate autocorrelation because variance is zero'
    def __str__(self):
        return self.msg
    
class VarianceUndefinedError(Exception):
    def __init__(self):
        self.msg = 'Cannot calculate standard deviation because sample size is less than 2'
    def __str__(self):
        return self.msg
    
class InvalidNumberOfColumnsError(Exception):
    def __init__(self, nparts, nexpected, line_num):
        self.msg = 'Number of values (%d) on line %d inconsistent with number of column headers (%d)' % (nparts, line_num, nexpected)
    def __str__(self):
        return self.msg

class ParamSummarizer(CommonFunctions):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Summarizes parameter sample contained in the param file, which is one
    of the two files output from an mcmc analysis (the other being the 
    tree file).
    
    """
    def __init__(self, opts):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes ParamSummarizer object.
        
        """
        CommonFunctions.__init__(self, opts)

    def interpolate(self, xx, x, y):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Given x and y, which define three reference points, (x[0],y[0]), 
        (x[1],y[1]) and (x[2],y[2]), find the point on the interpolating 
        polynomial corresponding to x-coordinate xx (scalar).
        
        """
        term0 = (xx - x[1])*(xx - x[2])/((x[0] - x[1])*(x[0] - x[2]))
        term1 = (xx - x[0])*(xx - x[2])/((x[1] - x[0])*(x[1] - x[2]))
        term2 = (xx - x[0])*(xx - x[1])/((x[2] - x[0])*(x[2] - x[1]))
        retval = term0*y[0] + term1*y[1] + term2*y[2]
        return retval
        
    def cumLagrange(self, which, x, y):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Given x and y, which define three reference points, (x[0],y[0]), 
        (x[1],y[1]) and (x[2],y[2]), find the integral under the interpolating
        polynomial for the first (which=1) or both (which=2) segments.
        
        """
        xx = x[which]
        
        x0 = x[0]
        x1 = x[1]
        x2 = x[2]
        
        y0 = y[0]
        y1 = y[1]
        y2 = y[2]
        
        psi0 = (x0 - x1)*(x0 - x2)
        psi1 = (x1 - x0)*(x1 - x2)
        psi2 = (x2 - x0)*(x2 - x1)
        
        xterm0 = (xx**3.0 - x0**3.0)/3.0
        xterm1 = (xx**2.0 - x0**2.0)/2.0
        xterm2 = xx - x0
        
        term0 = xterm0*(y0/psi0 + y1/psi1 + y2/psi2)
        term1 = xterm1*(y0*(x1 + x2)/psi0 + y1*(x0 + x2)/psi1 + y2*(x0 + x1)/psi2)
        term2 = xterm2*(y0*x1*x2/psi0 + y1*x0*x2/psi1 + y2*x0*x1/psi2)
        
        cum = term0 - term1 + term2
        return -cum
        
    def ss_simpsons(self, betas, means):
        """
        This approach uses Simpson's method to interpolate between beta values using the
        interpolation polynomial in Lagrange form. Simpson' method is described in most
        calculus textbooks, and (as used here) fits a parabola to every three consecutive
        points, using the area under the parabola as an approximation to the integral.
        This approach provides two estimates for the integral corresponding to each segment
        except for the first and last. We use the simple average of the two estimates when
        they are available.
        
        """
        nbetas = len(betas)
        if nbetas < 3:
            raise Exception("Must have at least 3 beta values to compute path sampling using Simpson's rule")
        marginal_like = 0.0
        for i in range(nbetas - 2):
            x = [betas[i], betas[i+1], betas[i+2]]
            y = [means[i], means[i+1], means[i+2]]
            # print '\ni = %d (dx = %.5f, dy = %.5f)' % (i,x[0] - x[2], y[0] - y[2])
            # print 'x%d <- c(%s)' % (i,','.join(['%.6f' % z for z in x]))
            # print 'y%d <- c(%s)' % (i,','.join(['%.6f' % z for z in y]))
            # print 'plot(x%d,y%d,type="b")' % (i,i)
            if i == 0:
                # no averaging on the first segment
                a = self.cumLagrange(1, x, y)
                marginal_like += a
                b = self.cumLagrange(2, x, y)
                marginal_like += (b - a)/2.0
            elif i == nbetas - 3:
                # no averaging on the last segment
                a = self.cumLagrange(1, x, y)
                marginal_like += a/2.0
                b = self.cumLagrange(2, x, y)
                marginal_like += (b - a)
            else:
                # average two estimates for each of the middle segments
                a = self.cumLagrange(1, x, y)
                b = self.cumLagrange(2, x, y)
                marginal_like += b/2.0
                # if i == 7 or i == 8:
                #     xxx = []
                #     yyy = []
                #     for z in range(21):
                #         xxx_value = x[0] + (x[2] - x[0])*float(z)/20.0
                #         yyy_value = self.interpolate(xxx_value, x, y)
                #         xxx.append(xxx_value)
                #         yyy.append(yyy_value)
                #     print 'xx%d <- c(%s)' % (i,','.join(['%.5f' % z for z in xxx]))
                #     print 'yy%d <- c(%s)' % (i,','.join(['%.5f' % z for z in yyy]))
                #     print 'plot(xx%d,yy%d,type="b")' % (i,i)
                #     print 'a     = %.5f' % a
                #     print 'b - a = %.5f' % (b - a)
                #     print '(x[0] - x[1])*(y[0] + y[1])/2 = ',(x[0] - x[1])*(y[0] + y[1])/2.0
                #     print '(x[1] - x[2])*(y[1] + y[2])/2 = ',(x[1] - x[2])*(y[1] + y[2])/2.0
        self.output(" %.8f Thermodynamic integration method (using Simpson's rule)" % (marginal_like))

    def ss_trapezoid(self, betas, means):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This method approximates the integral under the curve defined by betas
        (x-coordinates) and means (y-coordinates) using the trapezoid method
        (straight line interpolation). This is the method advocated by in the
        Lartillot and Phillippe (2006) paper that introduced the thermodynamic
        integration method to phylogenetics.
        
        """
        nbetas = len(betas)
        marginal_like = 0.0
        for i in range(nbetas):
            if i == 0:
                before = betas[0]
            else:
                before = betas[i-1]
            if i == nbetas - 1:
                after = betas[nbetas - 1]
            else:
                after = betas[i+1]
            diff = before - after
            if diff < 0.0:
                raise Exception('Phycas does not currently support path sampling from prior toward posterior')
            #print 'mean %d = %f (diff = %f)' % (i, means[i], diff)
            marginal_like += means[i]*diff/2.0
        self.output(' %.8f Path sampling method (using trapezoid rule)' % (marginal_like))
                
    def ss(self, betas, likes):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This method estimates the marginal likelihood using the product of 
        ratios (the stepping stones) bridging the gap between the posterior
        and the prior. Each ratio is estimated using importance sampling, with
        the importance distribution being the power posterior defined by 
        the smaller of the two beta values in the ratio. See Xie et al. (2009;
        Systematic Biology; submitted Jan. 2009) for details.
        
        """
        # Calculate marginal likelihood using Stepping Stone method
        # betas is a list of beta values ordered from the first sampled to the last sampled
        # likes is a map: i.e. likes[b] holds list of log-likelihoods sampled for beta value b
        # Assumes that betas[0] = 1.0 and betas[-1] = 0.0
        lnR = 0.0
        seR = 0.0
        nbetas = len(betas)
        if not self._betasSortedCorrectly(betas):
            #if betas[0] != 1.0 or betas[-1] != 0.0:
            #raise Exception('Stepping Stone method requires beta values to be ordered from 1.0 (first) to 0.0 (last)')
            raise Exception('Stepping Stone method requires beta values to be sorted from highest to lowest (0.0)')
        for i in range(1,nbetas):
            # find the difference between the two beta values for ratio i
            blarger = betas[i - 1]
            bsmaller = betas[i]
            beta_incr = blarger - bsmaller
            
            # find the maximum loglike
            loglikes = likes[bsmaller]
            n = len(loglikes)
            Lmax = max(loglikes)
            
            # find the log of the ratio of normalizing constants for ratio i
            tmp = 0.0
            for lnL in loglikes:
                tmp += math.exp(beta_incr*(lnL - Lmax))
            tmp /= float(n)
            lnRk = beta_incr*Lmax + math.log(tmp)
            lnR += lnRk
        
            # standard error calculation
            tmp1 = 0.0
            for lnL in loglikes:
                aa = math.exp(beta_incr*(lnL - Lmax))/tmp
                tmp1 += math.pow((aa - 1.0),2.0)
            seR += tmp1
        seR /= math.pow(float(n), 2.0)
        self.output(' %.8f Stepping-stone (SS) method using prior as reference distribution (se = %.8f)' % (lnR, seR))
        return lnR
        
    def _betasSortedCorrectly(self, betas):
        # beta values should not start at 1 (if so, probably using a file from former version of phycas)
        if betas[0] == 1.0:
            return False
        # last beta value should equal 0
        if betas[-1] != 0.0:
            return False
        # beta[i] should be smaller than beta[i-1] for i = 1, 2, ..., nbetas-1
        nbetas = len(betas)
        prev_beta = betas[0]
        for i in range(1,nbetas):
            curr_beta = betas[i]
            if curr_beta > prev_beta:
                return False
        return True
        
    def gss(self, betas, p):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This generalized Stepping Stone (SS) method estimates the 
        marginal likelihood using the product of ratios (the stepping stones) 
        bridging the gap between the posterior and the prior. Each ratio is 
        estimated using importance sampling, with the importance distribution
        being the power posterior defined by the smaller of the two beta
        values in the ratio. It differs from original SS in making use of 
        a reference distribution that may not be equivalent to the prior.
        
        """
        # Calculate marginal likelihood using Stepping Stone method.
        # betas is a list of beta values ordered from the first sampled to the last sampled.
        # likes is a map: i.e. likes[b] holds list of log-likelihoods sampled for beta value b.
        # priors is a map: i.e. priors[b] holds list of log-priors sampled for beta value b.
        # ref_dists is a map: i.e. ref_dists[b] holds list of log-working-priors 
        # sampled for beta value b. Assumes that betas[0] = 1.0 and betas[-1] = 0.0.
        likes = p['lnL']
        priors = p['lnPrior']
        ref_dists = p['lnRefDens']
        lnR = 0.0
        nbetas = len(betas)
        if not self._betasSortedCorrectly(betas):
            #if betas[0] != 1.0 or betas[-1] != 0.0:
            #raise Exception('Stepping Stone method requires beta values to be ordered from 1.0 (first) to 0.0 (last)')
            raise Exception('Stepping Stone method requires beta values to be sorted from highest to lowest (0.0)')
            
        self.output(' %10s %10s %10s %15s %15s' % ('b_(k-1)','beta_incr','n','lnRk','lnR(cum)'))
        for i in range(0,nbetas):
            #self.output('\nk = %d:' % i)
            
            # find the difference between the two beta values for ratio i
            if i > 0:
                blarger = betas[i - 1]
            else:
                blarger = 1.0
            bsmaller = betas[i]
            beta_incr = blarger - bsmaller
            
            # find the maximum term (lnL + lnp - lnwp)
            loglikes = likes[bsmaller]
            logpriors = priors[bsmaller]
            logwpriors = ref_dists[bsmaller]
            n = len(loglikes)
            etak = max([(lnL + lnp - lnwp) for lnL,lnp,lnwp in zip(loglikes,logpriors,logwpriors)])
            
            # find the log of the ratio of normalizing constants for ratio i
            tmp = 0.0
            for lnL,lnp,lnwp in zip(loglikes,logpriors,logwpriors):
                log_term = beta_incr*(lnL + lnp - lnwp - etak)
                #self.output('--> %15.5f %15.5f %15.5f %15.5f %15.5f' % (lnL,lnp,lnwp,etak,log_term))
                tmp += math.exp(log_term)
            tmp /= float(n)
            lnRk = beta_incr*etak + math.log(tmp)
            lnR += lnRk
            self.output(' %10.3f %10.3f %10d %15.6f %15.6f' % (bsmaller, beta_incr, n, lnRk, lnR))
        
        self.output(' %.8f Generalized Stepping Stone method' % lnR)
        return lnR
        
    def autocorr_ess(self, values):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Computes the lag-1 autocorrelation (r) for the supplied values list.
        Also computes the effective sample size (ess). Returns the tuple 
        (r,ess).
        
        """
        nvalues = len(values)
        n = float(nvalues)
        
        # calculate the mean
        m = sum(values)/n
        
        # calculate the variance
        ss = 0.0
        for x in values:
            ss += x**2.0
        var = (ss - n*m*m)/(n - 1.0)
        
        # calculate the covariance
        cov = 0.0
        for i in range(nvalues - 1):
            x = values[i] - m
            y = values[i+1] - m
            cov += x*y
        cov /= n - 1.0
        if var <= 0.0:
            raise VarianceZeroError()
        r = cov/var
        ess = n*(1.0 - r)/(1.0 + r)
        return (r,ess)
    
    def marginal_likelihood(self, headers, lines, burnin):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Estimates the marginal likelihood for the Stepping Stone/Thermodynamic
        Integration case and outputs the autocorrelation and effective sample 
        size for each parameter/beta combination. The supplied list headers 
        holds the column headers from the param file. The parameter names are
        taken from this. The supplied lines is a list of lines from the param
        file. The supplied burnin indicates the number of initial lines to 
        ignore (usually will be 1 because the starting values of all 
        parameters is always output and should always be skipped).
        
        """
        marglike = None
        
        # Create params dictionary of dictionaries such that, for example,
        # params['lnL'][1.0] is a list of all log-likelihood values sampled for 
        # beta = 1.0 (in the order in which they were sampled)
        betas = []
        params = {}
        for h in headers:
            params[h] = {}
        row_start = 2 + burnin  # first line holds the ID, second line are the headers
        curr_beta = None
        for i,line in enumerate(lines[row_start:]):
            parts = line.split()
            
            # each line (except first, ID, line) should have the same number of items
            if len(parts) != len(headers):
                raise InvalidNumberOfColumnsError(len(parts), len(headers), i + row_start + 1)
                
            # expecting beta to be the second item on each line (cycle is first item on each line)
            beta = float(parts[1])
            if curr_beta is None or curr_beta != beta:
                # start new beta value
                curr_beta = beta
                betas.append(beta)
                for h,x in zip(headers,parts):
                    params[h][beta] = [float(x)]
            else:
                # continuing previous beta value
                for h,x in zip(headers,parts):
                    params[h][beta].append(float(x))
                                        
        # Output first-order autocorrelation for each parameter (as well as the log-likelihood,
        # prior and, if present, the working prior) for each beta value separately
        self.output('\nAutocorrelations (lag 1):\n')
        self.output('%15s%s' % ('beta ->',' '.join(['%12.5f' % b for b in betas])))
        #print "@@@@@@@@@@@@@@@@@@@@ headers: |",["%s|" % h for h in headers[2:]]
        for h in headers[2:]:   # skip 'Gen' and 'beta' headers
            s = []
            for b in betas:
                p = params[h][b]
                #print 'trying h = %s, b = %g...' % (h,b)
                try:
                    r,ess = self.autocorr_ess(p)
                except VarianceZeroError:
                    s.append('%12s' % '---')
                except OverflowError:
                    s.append('%12s' % '---')
                else:
                    s.append('%12.5f' % r)
            self.output('%15s%s' % (h,' '.join(s)))

        # Output effective sample size for each parameter (as well as the log-likelihood,
        # prior and, if present, the working prior) for each beta value separately
        actual_sample_sizes = [len(params['lnL'][b]) for b in betas]
        self.output('\nEffective and actual sample sizes:\n')
        self.output('%15s%s' % ('beta ->', ' '.join(['%12.5f' % b for b in betas])))
        self.output('%15s%s' % ('actual', ' '.join(['%12d' % n for n in actual_sample_sizes])))
        sample_size_discrepancy = False
        for h in headers[2:]:   # skip 'Gen' and 'beta' headers
            s = []
            for b_index,b in enumerate(betas):
                p = params[h][b]
                ss_ok = (len(p) == actual_sample_sizes[b_index])
                if ss_ok:
                    asterisk = ''
                else:
                    sample_size_discrepancy = True
                    asterisk = '*'
                try:
                    r,ess = self.autocorr_ess(p)
                except VarianceZeroError:
                    s.append('%12s%s' % ('---',asterisk))
                else:
                    s.append('%12.1f%s' % (ess,asterisk))
            self.output('%15s%s' % (h,' '.join(s)))
            
        if sample_size_discrepancy:
            self.warning('* indicates discrepancy between actual sample size for lnL and at least one other parameter')         

        if headers[4] == 'lnRefDens':
            # Estimate marginal likelihood using generalized Stepping Stone (SS) method
            self.output('\nMarginal likelihood estimate:')
            try:
                marglike = self.gss(betas, params)
            except Exception,e:
                self.output(' %s' % e.message)
        else:
            # Estimate marginal likelihood using Thermodynamic Integration (TI) and classical
            # Stepping Stone (SS) method.
            
            # Compute means of log-likelihoods for each beta value (used for ps calculation)
            means = []
            for b in betas:
                p = params['lnL'][b]
                m = sum(p)/float(len(p))
                means.append(m)
            self.output('\nMean log-likelihood for each value of beta used\nfor marginal likelihood estimation:\n')
            self.output('%12s %12s' % ('beta','mean lnL'))
            for b,m in zip(betas, means):
                s = ''
                if b == 1.0:
                    s = '(posterior)'
                elif b == 0.0:
                    s = '(prior)'
                self.output('%12.5f %12.5f %s' % (b,m,s))

            self.output('\nMarginal likelihood estimates:')
            try:
                marglike = self.ss(betas, params['lnL'])
            except Exception,e:
                self.output(' %s' % e.message)
            try:
                self.ss_trapezoid(betas, means)
            except Exception,e:
                self.output(' %s' % e.message)
            #try:
            #    self.ss_simpsons(betas, means)
            #except Exception,e:
            #    self.output(' %s' % e.message)
        return marglike
                    
    def harmonic_mean(self, v):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Calculate marginal likelihood using the harmonic mean method using 
        log-likelihoods supplied in v.
        
        """
        nignored = 0
        n = len(v)
        min_lnL = min(v)
        sum_diffs = 0.0
        for lnl in v:
            diff = lnl - min_lnL
            if diff < 500.0:
                sum_diffs += math.exp(-diff)
            else:
                nignored += 1
        log_harmonic_mean = math.log(n) + min_lnL - math.log(sum_diffs)
        if nignored > 0:
            self.warning('ignoring %d sampled log-likelihoods in harmonic mean calculation' % nignored)
        self.output('Log of marginal likelihood (harmonic mean method) = %f' % log_harmonic_mean)
        return log_harmonic_mean
            
    def summary_stats(self, v, cutoff=95):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Computes the following summary statistics for the supplied vector v:
        first-order autocorrelation (lag=1), effective sample size, lower 
        credible interval bound, upper credible interval bound, minimum, 
        maximum, sample mean, and sample standard deviation (i.e. divide by 
        n-1). The value of cutoff is the percentage to use for the credible
        interval. If v is None, returns tuple of header strings. If v is not
        None, returns tuple of summary statistics.
        
        """
        if v is None:
            h = ('autocorr', 'ess', 'lower %d%%' % int(cutoff), 'upper %d%%' % int(cutoff), 'min', 'max', 'mean', 'stddev')
            return h
            
        if len(v) < 2:
            raise VarianceUndefinedError()
        s = []
        
        # compute autocorr and ess
        try:
            r,ess = self.autocorr_ess(v)
        except VarianceZeroError:
            raise VarianceZeroError()
        else:
            s.extend((r,ess))
            
        # compute lower and upper
        v.sort()
        n = float(len(v))
        p = float(cutoff)/100.0
        lower_at = int(math.ceil(n*(1 - p)))
        upper_at = int(math.ceil(n*p))
        lower = v[lower_at]
        upper = v[upper_at]
        s.extend((lower, upper))
        
        # compute min and max
        s.extend((v[0], v[-1]))
        
        # compute mean and stddev
        mean = sum(v)/n
        ss = sum([x**2 for x in v])
        var = (ss - n*mean**2)/(n - 1.0)
        sd = math.sqrt(var)
        s.extend((mean, sd))
        
        return tuple(s)
    
    def std_summary(self, headers, lines, burnin):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Produces a table of summary statistics for each parameter using the
        data in a param file. This is the companion to marginal_likelihood for
        the standard mcmc case. The supplied list headers holds the column 
        headers from the param file. The parameter names are taken from this.
        The supplied lines is a list of lines from the param file. The 
        supplied burnin indicates the number of initial lines to ignore 
        (usually will be 1 because the starting values of all parameters is 
        always output and should always be skipped). See documentation for
        the summary_stats function for a list of summary statistics computed.
        
        """
        # Create params dictionary such that, for example, params['lnL'] is a list of 
        # all log-likelihood values (in the order in which they were sampled)
        params = {}
        for h in headers:
            params[h] = []
        row_start = 2 + burnin  # first line ID, second line headers
        for i,line in enumerate(lines[row_start:]):
            parts = line.split()
            if len(parts) != len(headers):
                raise InvalidNumberOfColumnsError(len(parts), len(headers), i + row_start + 1)
            for h,x in zip(headers,parts):
                params[h].append(float(x))

        # Output summary statistics for each parameter (and the log-likelihood)
        self.output('\nSummary statistics:\n')
        stats_headers = ('param','n') + self.summary_stats(None)
        sz = len(stats_headers)
        gss = '%20s' + '%15s'*(sz-1) + '\n'
        self.output(gss % stats_headers)
        for h in headers[1:]:   # skip 'Gen' header
            v = params[h]
            try:
                stats = (h,len(v)) + self.summary_stats(v)
                gss = '%20s' + '%15d' + '%15.5f'*(sz - 2)
                self.output(gss % stats)
            except (VarianceZeroError,VarianceUndefinedError):
                gss = '%20s' + '%15s'*(sz-1) + '\n'
                sub = tuple(['---']*sz)
                self.output(gss % sub)
        marglike = None
        self.output()
        marglike = self.harmonic_mean(params['lnL'])
        return marglike
        
    def calcLogHM(self, vect_of_log_values):

        logn = math.log(float(len(vect_of_log_values)))
        logLmin = min(vect_of_log_values)
        sum_log_diffs = 0.0
        for logx in vect_of_log_values:
            sum_log_diffs += math.exp(logLmin - logx)
        loghm = logn + logLmin - math.log(sum_log_diffs)
    
        return loghm

    def _cpoOpenRFile(self):
        if self._cpoRFile is None:
            sp = self.optsout.cpoplot 
            self._cpoRFilename = sp._getFilename()
            self._cpoRFile = sp.open(self.stdout)
        return self._cpoRFile

    def _cpoOpenInfoFile(self):
        if self._cpoInfoFile is None:
            sp = self.optsout.cpoinfo
            self._cpoInfoFilename = sp._getFilename()
            self._cpoInfoFile = sp.open(self.stdout)
        return self._cpoInfoFile

    def cpo_summary(self, lines, burnin):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
            Produces an R file containing commands for producing a plot of CPO
            (Conditional Predictive Ordinates). This plot has site position as the
            x-coordinate and site CPO as the y-coordinate. The title of the plot
            contains the summary measure (sum of CPO over all sites).
            
            """
        self.output('\nCPO analysis')
        
        # Each line in lines comprises nsites log-site-likelihood values (one sample from the chain)
        nsites = len(lines[0].split())

        # Create nsites lists, each of which holds all log-site-likelihood values sampled from one site
        loglikes = [[] for i in range(nsites)]
        for line in lines[burnin:]:
            parts = line.split()
            assert len(parts) == nsites
            for i,logx in enumerate(parts):
                loglikes[i].append(float(logx))

        # Create the default partition if none has been defined
        partition.validate(nsites)
        
        # Compute the log-CPO measure for each site, and the sum over all sites
        # Sites that have been excluded will have lognm = 0.0 and thus will not
        # contribute to total_cpo
        cpovect = [0.0]*nsites  # This vector has sites grouped by partition subset!
        cpomap = {}             # key is site index, value is log(cpo)
        total_cpo = 0.0
        k = 0
        subset_names     = []  # list of subset names
        subset_sitelists = []  # list of subset sitelists
        for subset_index,(subset_name,subset_sitelist,subset_model) in enumerate(partition.subset):
            subset_names.append(subset_name)
            subset_sitelists.append(subset_sitelist)
            print 'CPO: processing subset = %s...' % subset_name
            for i in subset_sitelist:
                loghm = self.calcLogHM(loglikes[i-1])
                total_cpo += loghm
                cpovect[k] = loghm
                cpomap[i] = loghm
                k += 1 

        self.output('LPML (Log PseudoMarginal Likelihood) = %.5f' % total_cpo)
        
        # Identify the worst sites in terms of CPO, placing these in cpo_worst
        # The list cpo_worst lists sites by partition subset, not in original order
        cpo_worst = []
        for i in range(nsites):
            if cpovect[i] < 0.0:
                cpo_worst.append((i,cpovect[i]))
        cpo_worst.sort(cmp=lambda x,y: cmp(x[1], y[1]))
        nincluded = len(cpo_worst)
        last_of_worst = int(math.ceil(self.opts.cpo_cutoff*nincluded))
        cpo_worst[last_of_worst:] = []
        
        # Identify the worst sites again, placing these in the list cpo_worst_orig
        # The list cpo_worst_orig lists sites in their original order, not by partition subset
        # Sites with cpo higher than best_of_the_worst will have value None
        best_of_the_worst = cpo_worst[-1][1]
        cpo_worst_orig = [None]*nsites
        for i in cpomap.keys():
            if cpomap[i] <= best_of_the_worst:
                cpo_worst_orig[i-1] = cpomap[i]
        
        # Create a mask showing which sites are in the worst category
        # This mask string can be aligned to the original data matrix making it easy to identify the worst sites
        mask = ['-']*nsites
        for i,v in enumerate(cpo_worst_orig):
            if v is not None:
                mask[i] = '*'
        maskstr = ''.join(mask)

        # Create R and info files. The R file will produce a plot (using R) that plots log CPO values as impulses
        # and sorted by partition subset (with sites in the worst CPO category colored red). The info file is
        # tab-delimited and designed to contain the same information in way that could be easily read and processed
        # for other programs/uses
        self._cpoRFilename = None
        self._cpoRFile = None
        self._cpoInfoFilename = None
        self._cpoInfoFile = None
        nsubsets = len(subset_names)
        if bool(self.optsout.cpoplot):
            try:
                # Create CPO info file
                self._cpoOpenInfoFile()
                
                # Write mask
                self._cpoInfoFile.write('Mask showing worst (lowest %.1f%% CPO values) sites as *\n' % (100.0*self.opts.cpo_cutoff,))
                self._cpoInfoFile.write('  (order of sites is that of the original data matrix)\n')
                self._cpoInfoFile.write('BEGIN_MASK\n')
                self._cpoInfoFile.write('  %s\n' % maskstr)
                self._cpoInfoFile.write('END_MASK\n')

                subset_logCPOs = []
                for i in range(nsubsets):
                    subset_logCPOs.append([])
                tally = [0]*nsubsets
                
                # Write table of log(CPO) values
                self._cpoInfoFile.write('\nBEGIN_LOG_CPO_TABLE\n')
                self._cpoInfoFile.write('%12s\t%12s\t%12s\t%12s\n' % ('site', 'log(CPO)', 'subset', 'worst'))
                for site_index in range(nsites):
                    site_number = site_index+1
                    
                    inworst = 0
                    if cpo_worst_orig[site_index] is not None:
                        inworst = 1
                    
                    which_subset = None
                    for subset_index in range(nsubsets):
                        if site_number in subset_sitelists[subset_index]:
                            which_subset = subset_names[subset_index]
                            subset_logCPOs[subset_index].append(cpomap[site_number])
                            tally[subset_index] += inworst
                            break
                    assert which_subset is not None
                    self._cpoInfoFile.write('%12d\t%12.5f\t%12s\t%12d\n' % (site_number,cpomap[site_number], which_subset, inworst))
                self._cpoInfoFile.write('END_LOG_CPO_TABLE\n')

                # Write subset summary
                self._cpoInfoFile.write('\nBEGIN_SUBSET_SUMMARY\n')
                self._cpoInfoFile.write('%12s\t%12s\t%12s\t%12s\t%12s\n' % ('subset', 'nsites', 'nworst', 'mean log(CPO)', 'log(mean CPO)'))
                for i in range(nsubsets):
                    n = len(subset_logCPOs[i])
                    mean_log = sum(subset_logCPOs[i])/n
                    max_logCPO = max(subset_logCPOs[i])
                    sum_diffs = sum([math.exp(logCPO - max_logCPO) for logCPO in subset_logCPOs[i]])
                    log_mean = max_logCPO + math.log(sum_diffs) - math.log(n)
                    self._cpoInfoFile.write('%12s\t%12d\t%12d\t%12.5f\t%12.5f\n' % (subset_names[i], n, tally[i], mean_log, log_mean))
                self._cpoInfoFile.write('END_SUBSET_SUMMARY\n')
                
                self._cpoInfoFile.close()
                
                #numWorst = len(cpo_worst)
                #meanLogCPOWorst = 0.0
                #for i in range(numWorst):
                #    meanLogCPOWorst += cpo_worst[i][1]
                #meanLogCPOWorst /= numWorst
                #print 'Number of worst sites:',numWorst
                #print 'Mean of log(CPO) of worst sites:',meanLogCPOWorst

                # Create CPO R plot file
                self._cpoOpenRFile()
                self._cpoRFile.write('# Plot log(CPO) for each site (sites grouped by partition subset)\n')
                self._cpoRFile.write('x = c(%s)\n' % ','.join(['%d' % x for x in range(nsites)]))
                self._cpoRFile.write('y = c(%s)\n' % ','.join(['%g' % y for y in cpovect]))
                self._cpoRFile.write('colvec = rep("black",%d)\n' % nsites)
                self._cpoRFile.write('z = c(%s)\n' % ','.join(['%d' % (zz[0]+1) for zz in cpo_worst]))
                self._cpoRFile.write('colvec[z] = "red"\n')
                self._cpoRFile.write("plot(x, y, type='h', col=colvec, main='Overall CPO = %.5f', xlab='Site', ylab = 'CPO')\n" % total_cpo)
                
                cum = 0
                self._cpoRFile.write('abline(v=1)\n')
                for subset_index,(subset_name,subset_sitelist,subset_model) in enumerate(partition.subset):
                    cum += len(subset_sitelist)
                    self._cpoRFile.write('abline(v=%d)\n' % cum)
            
            finally:
                if self._cpoRFile:
                    self.optsout.cpoplot.close()
                if self._cpoInfoFile:
                    self.optsout.cpoinfo.close()
    
    def cpo_summary_obsolete(self, lines, burnin):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Produces an R file containing commands for producing a plot of CPO
        (Conditional Predictive Ordinates). This plot has site position as the
        x-coordinate and site CPO as the y-coordinate. The title of the plot
        contains the summary measure (sum of CPO over all sites).
        
        """
        self.output('\nCPO analysis')
        
        # Each line in lines comprises nsites log-site-likelihood values (one sample from the chain)
        nsites = len(lines[0].split())
        loglikes = [[] for i in range(nsites)]

        # Create nsites lists, each of which holds all log-site-likelihood values sampled from one site
        for line in lines:
            parts = line.split()
            for i,logx in enumerate(parts):
                loglikes[i].append(float(logx))
            
        # Compute the log-CPO measure for each site, and the sum over all sites
        # Sites that have been excluded will have lognm = 0.0 and thus will not
        # contribute to total_cpo
        yvect = []
        total_cpo = 0.0
        for i in range(nsites):
            loghm = self.calcLogHM(loglikes[i])
            total_cpo += loghm
            yvect.append(loghm)
        self.output('Model CPO = %.5f' % total_cpo)
        
        # Identify the worst sites in terms of CPO, placing these in cpo_worst
        cpo_worst = []
        for i in range(nsites):
            if yvect[i] < 0.0:
                cpo_worst.append((i,yvect[i]))
        cpo_worst.sort(cmp=lambda x,y: cmp(x[1], y[1]))
        nincluded = len(cpo_worst)
        last_of_worst = int(math.ceil(self.opts.cpo_cutoff*nincluded))
        cpo_worst[last_of_worst:] = []
        
        # Create a mask showing which sites are in the worst category
        mask = ['-']*nsites
        for i,j in cpo_worst:
            mask[i] = '*'
        for i in range(nsites):
            if yvect[i] == 0.0:
                mask[i] = 'x'
        maskstr = ''.join(mask)
        # begin again here: need to print out maskstr in output...

        # Use nearest neighbor smoother to produce a more easily-interpreted plot
        
        # Smoothing algorithm:
        #
        #    ysmoothed[i] = (sum_j w_j y_j) / (sum_j w_j)
        #
        # where: j is the set [x_{i-m}, ..., x_{i-1}, x_i, x_{i+1}, ..., x_{i+m}]
        #        w_j = exp{-(x_i - x_j)^2/(2*var)}
        #        x_i is the position of site i (i.e. for the 5th site, x_i = 5.0)
        # note: setting m to nsites seems to work pretty well, so might as well
        #       just use this vale of m and save the user having yet another
        #       setting to worry about.
        
        m = nsites  # m is the number of neighbors to each side of site i 
        sd = self.opts.cposmooth    # use standard deviation sd for Gaussian weights
        var = math.pow(sd,2.0)
        denom = 2.0*var
        ysmoothed = []
        for i in range(nsites):
            low = i - m
            high = i + m
            if low < 0:
                low = 0
            if high > nsites - 1:
                high = nsites - 1
            yvalues = yvect[low:high+1]
            xvalues = range(low, high+1)
            assert len(xvalues) == len(yvalues), '%d != %d for i = %d' % (len(xvalues),len(yvalues),i)
            weights = [math.exp(-math.pow(-(float(x) - float(i)), 2.0)/denom) for x in xvalues]
            wy = [w*y for w,y in zip(weights,yvalues)]
            yavg = sum(wy)/sum(weights)
            ysmoothed.append(yavg)
            
        self._cpoRFilename = None
        self._cpoRFile = None

        if bool(self.optsout.cpoplot):
            try:
                self._cpoOpenRFile()
                self._cpoRFile.write('# plot of log(CPO) across sites\n')
                #self._cpoRFile.write('quartz(bg="white")\n')
                self._cpoRFile.write('x = c(%s)\n' % ','.join(['%d' % x for x in range(nsites)]))
                self._cpoRFile.write('y = c(%s)\n' % ','.join(['%g' % y for y in yvect]))
                #self._cpoRFile.write("plot(x, y, type='l', main='Overall CPO = %.5f', xlab='Site', ylab = 'CPO')\n" % total_cpo)
                #self._cpoRFile.write('quartz(bg="white")\n')
                self._cpoRFile.write('colvec = rep("black",%d)\n' % nsites)
                self._cpoRFile.write('z = c(%s)\n' % ','.join(['%d' % (zz[0]+1) for zz in cpo_worst]))
                self._cpoRFile.write('colvec[z] = "red"\n')
                self._cpoRFile.write("plot(x, y, type='h', col=colvec, main='Overall CPO = %.5f', xlab='Site', ylab = 'CPO')\n" % total_cpo)
                self._cpoRFile.write('\n# plot of estimated probability density of log(CPO) values\n')
                self._cpoRFile.write('quartz(bg="white")\n')
                self._cpoRFile.write('plot(density(y))\n')
                self._cpoRFile.write('\n# Smoothed plot of log(CPO) across sites\n')
                self._cpoRFile.write('quartz(bg="white")\n')
                self._cpoRFile.write('ysmoothed = c(%s)\n' % ','.join(['%g' % y for y in ysmoothed]))
                self._cpoRFile.write("plot(x, ysmoothed, type='l', main='Smoothed using standard deviation %g', xlab='Site', ylab = 'CPO')\n" % sd)
            finally:
                if self._cpoRFile:
                    self.optsout.cpoplot.close()
        
    def handleFile(self, fn):
        burnin = self.opts.burnin
        lines = open(fn, 'r').readlines()
        marglike = None
        if len(lines) < 3 + burnin:
            self.output("File '%s' does not look like a parameter file (too few lines)")
        else:
            headers = lines[1].split()
            if headers[1] == 'beta':
                try:
                    marglike = self.marginal_likelihood(headers, lines, burnin)
                except InvalidNumberOfColumnsError, e:
                    print e
            else:
                marglike = self.std_summary(headers, lines, burnin)
        return marglike
                
    def handleCPOFile(self, cpofn):
        burnin = self.opts.burnin
        lines = open(cpofn, 'r').readlines()
        if len(lines) < 3 + burnin:
            self.output("File '%s' does not look like a site likelihood file (too few lines)" % cpofn)
        else:
            self.cpo_summary(lines, burnin)
        
    def run(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Reads the contents of the param file and decides whether to use the 
        marginal_likelihood or std_summary functions to summarize the data
        therein. The marginal_likelihood function is called if the second
        header is "beta"; otherwise, std_summary is called. If a cpofile
        name has been specified, computes conditional predictive ordinates.
        
        """
        fn = self.opts.file
        cpofn = self.opts.cpofile
        self.phycassert(len(fn) + len(cpofn) > 0, "Must specify either a parameter file ('file') or file of site likelihoods ('cpofile') when invoking the sump command")
        if len(fn) > 0:
            marglike = self.handleFile(fn)
        if len(cpofn) > 0:
            self.handleCPOFile(cpofn)
        return marglike
