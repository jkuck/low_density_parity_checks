#Code for checking SAT solver speed when using various
#low density parity check matrix constructions

#import subprocess
import commands
import numpy as np
import copy
from decimal import Decimal
from itertools import combinations
import operator as op
import time

import matplotlib
matplotlib.use('Agg') #prevent error running remotely
import matplotlib.pyplot as plt
#Directory the open-wbo_static executable is in
#Install wbo from http://sat.inesc-id.pt/open-wbo/installation.html
WBO_DIRECTORY = '/atlas/u/jkuck/software/open-wbo'

#Directory with the MaxHS executable
#Download from http://www.maxhs.org/downloads.html
MAX_HS_DIRECTORY = '/atlas/u/jkuck/software/MaxHS-3.0/build/release/bin'

#Directory containing the cryptominisat5 executable
#installation instructions: https://github.com/msoos/cryptominisat
CRYPTOMINISAT5_DIRECTORY = '/Users/jkuck/software/cryptominisat-5.0.1/build'
#CRYPTOMINISAT5_DIRECTORY = '/atlas/u/jkuck/software/cryptominisat-5.0.1/build'

SAT_SOLVER = "MAX_HS"
#SAT_SOLVER = "CRYPTOMINISAT5"
#SAT_SOLVER = "WBO" #something strange seems to happen using WBO, also seems slower than MAX_HS


def nCr(n, r):
    '''
    Compute n choose r
    https://stackoverflow.com/questions/4941753/is-there-a-math-ncr-function-in-python
    '''
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, xrange(n, n-r, -1))
    denom = reduce(op.mul, xrange(1, r+1))
    return numer//denom



class SAT_problem:
    def __init__(self, nbvar, nbclauses, clauses):
        '''
        Represent an unweighted SAT problem  

        Inputs:
        - nbvar: int, number of variables 
        - nbclauses: int, the number of clauses
        - clauses: list of strings, each string represents a clause as a sequence of non-zero
            integer numbers between -nbvar and nbvar and ends with 0. Positive numbers denote 
            the corresponding variables. Negative numbers denote the negations of the 
            corresponding variables. (as in http://www.maxhs.org/docs/wdimacs.html)
        '''
        self.nbvar = nbvar
        self.nbclauses = nbclauses
        self.clauses = clauses


    def convert_XOR_to_chunked_OR(XOR_constraint, max_clause_size):
        '''
    
        Inputs:
        - XOR_constraint: list of ints, representing an XOR constraint of variables specified by ints
        - max_clause_size: int, maximum number of original variables allowed per OR clause
    
        Outputs:
        - OR_clauses: list of list of ints, each list of ints represents an OR clause over the 
            variables referenced by the ints it contains
    
        '''
        return None

        

    def add_parity_constraints(self, m, f, use_XOR_clauses):
        '''
        Add low density parity constraints to the SAT problem, randomly hashing satisfying 
        solutions into m dimensional bins (m<n) by setting every element in the matrix A to 1
        with probabiliy f as described here: https://cs.stanford.edu/~ermon/papers/SparseHashing-revised.pdf

        Inputs:
        - m: int, the dimension of the output hash space (rows in matrix A)
        - f: float, probability of setting every entry in A to 1
        - use_XOR_clauses: bool, if True encode XOR clauses directl
        Outputs:

        '''
        parity_clauses = []
        parity_clause_count = 0
        assert(m < self.nbvar)
        for row in range(m): #add new clauses (parity constraints) for each row in A
            #sample which elements in the row are 1
            col_indices_of_ones = []
            ones_count = 0
            for col in range(1, self.nbvar + 1):
                if np.random.rand() < f:
                    col_indices_of_ones.append(col)
                    ones_count += 1
            #DEBUG
            for i in range(0, ones_count+1, 2):
                parity_clause_count += nCr(ones_count, i)

            #END DEBUG

            #require that an odd number of variables corresponding to col_indices_of_ones are true
            #e.g. not (x1 and x2 and x3 and x4 and not x5) = (not x1 or not x2 or not x3 or not x4 or x5)
            #enumerate all the even numbers less than or equal to the number of 1s in this row
            for number_ones in range(0, len(col_indices_of_ones) + 1, 2):
                must_be_0_constraints = combinations(col_indices_of_ones, number_ones)
                for constraint in must_be_0_constraints:
                    #constraint is a list of integers, each in 1 to n (n=self.nbvar)
                    #that specifies negated variables in the clause
                    clause = ''
                    neg_var_length = 0
                    pos_var_length = 0
                    for var in range(1, self.nbvar + 1):
                        if var in constraint:
                            clause += str(-var) + ' '
                            neg_var_length += len(str(var)) + 2
                        else:
                            clause += str(var) + ' '
                            pos_var_length += len(str(var)) + 1
                    clause += '0'
#                    assert(neg_var_length + pos_var_length == self.nbvar)
                    assert(len(clause) == pos_var_length + neg_var_length + 1), (len(clause), pos_var_length + neg_var_length + 1)
                    parity_clauses.append(clause)

        assert(len(parity_clauses) == parity_clause_count)
        self.nbclauses += parity_clause_count
        self.clauses.extend(parity_clauses)




def read_SAT_problem(problem_filename):
    '''
    Read in a SAT problem specified in WDIMACS format
    (http://www.maxhs.org/docs/wdimacs.html)

    Inputs:
    - problem_filename: string, the filename of the problem in WDIMACS format

    Outputs:
    - sat_problem: type SAT_problem, the SAT problem we read from problem_filename
    '''
    check_clause_count = 0
    clauses = []
    
    f = open(problem_filename, 'r')
    for line in f:
        if line[0] == 'p': #parameters line
            params = line.split()
            assert(params[1] == 'cnf') #we should be reading an unweighted SAT problem
            nbvar = int(params[2]) #number of variables
            nbclauses = int(params[3]) #number of clauses
        elif line[0] != 'c': #line isn't commented out, should be a clause
            clauses.append(line)
            check_clause_count += 1
        else:
            assert(line[0] == 'c')
    f.close()
    
    assert(check_clause_count == nbclauses), (check_clause_count, nbclauses) #did we read the correct number of clauses?
    sat_problem = SAT_problem(nbvar, nbclauses, clauses)
    return sat_problem

def write_SAT_problem(problem_filename, sat_problem, problem_type):
    '''
    Write a SAT problem to a file in WDIMACS format
    (http://www.maxhs.org/docs/wdimacs.html)

    Inputs:
        - problem_filename: string, the filename we will write the SAT problem to
        - sat_problem: type SAT_problem, the SAT problem we will write to problem_filename
        - problem_type: 'SAT' or 'weighted_MaxSat'

    Outputs:
        None, but we will write a file. Note that if a file with problem_filename exists,
        it will be erased and overwritten.
    '''
    assert(problem_type in ['SAT', 'weighted_MaxSat'])
    f = open(problem_filename, 'w')
    #write parameters line
    if problem_type == 'SAT':
        f.write("p cnf %d %d\n" % (sat_problem.nbvar, sat_problem.nbclauses))       
    else:
        f.write("p wcnf %d %d %d\n" % (sat_problem.nbvar, sat_problem.nbclauses, sat_problem.top))

    for clause in sat_problem.clauses:
        if problem_type == 'SAT':
            f.write(clause)
        else:
            f.write("%d %s" % (sat_problem.top, clause))

    if problem_type == 'weighted_MaxSat':
        if len(sat_problem.soft_clauses) > 0:
            f.write("\n")
        for soft_clause in sat_problem.soft_clauses:
            f.write(soft_clause)

    f.close()
    return None


def solve_SAT(sat_problem):
    '''
    Call a SAT solver to solve the specified SAT problem

    Inputs:
        - sat_problem: type SAT_problem, the SAT problem to solve

    Outputs:
        - satisying_solution: list of ints, each entry is either 1 or -1. 
        satisying_solution[i] is the value that variable i takes in the 
        satisfying solution to the SAT problem (1 indicates True, -1 False)
        - float: time in seconds spent solving the SAT problem
    '''
    satisying_solution = []

    write_SAT_problem('./temp_SAT_file.txt', sat_problem, problem_type='SAT')

    t0 = time.time()
    if SAT_SOLVER == "WBO":
        (status, output) = commands.getstatusoutput("%s/open-wbo_static ./temp_SAT_file.txt" % WBO_DIRECTORY)
    elif SAT_SOLVER == "CRYPTOMINISAT5":
        (status, output) = commands.getstatusoutput("%s/cryptominisat5 --verb 0 ./temp_SAT_file.txt" % CRYPTOMINISAT5_DIRECTORY)        
    else:
        assert(SAT_SOLVER == "MAX_HS")
        print 'hi'
        print "%s/maxhs ./temp_SAT_file.txt" % MAX_HS_DIRECTORY
        (status, output) = commands.getstatusoutput("%s/maxhs ./temp_SAT_file.txt" % MAX_HS_DIRECTORY)
    t1 = time.time()
    solver_time = t1 - t0

    print output

    zero_count = 0 #cryptominisat5 appends a 0 at the end of the solution, make sure we only get one
    for line in output.splitlines():
        if line[0] == 'v': #find the line in the output containing variable values in the solution
            params = line.split()
#            assert(len(params) == sat_problem.nbvar + 1), (len(params), sat_problem.nbvar+1, params)
            for i in range(1, len(params)):
                if int(params[i]) > 0:
                    satisying_solution.append(1)
                elif(int(params[i]) < 0):
                    satisying_solution.append(-1)
                else:
                    assert(int(params[i]) == 0)
                    zero_count+=1

    assert(zero_count == 0 or zero_count == 1)
    assert(len(satisying_solution) == sat_problem.nbvar + 1), (len(satisying_solution), sat_problem.nbvar+1)

    return satisying_solution, solver_time



if __name__=="__main__":
    sat_problem = read_SAT_problem("SAT_problems_cnf/%s" % "c499.isc")
#    sat_problem.add_parity_constraints(m = 40, f=0.0, use_XOR_clauses=True)
    solution, time = solve_SAT(sat_problem)
    print time
    sleep(1)


    #write_SAT_problem('temp_SAT_file.txt', sat_problem)
    #solve_SAT(sat_problem, problem_type='SAT')

    #dictionary with key: values of model_name: number of satisfying solutions
    model_counts = {
        "log-1.cnf" : 564153552511417968750,
        "log-2.cnf" : 32334741710,
        "log-3.cnf" : 279857462060,
        "log-4.cnf" : 23421510324076617565622131248,
        "log-5.cnf" : 724152621485436659540387630662916505600,
        "tire-1.cnf": 726440820,
        "tire-2.cnf": 738969640920,
        "tire-3.cnf": 222560409176,
        "tire-4.cnf": 103191650628000,
        "ra.cnf" : 18739277038847939886754019920358123424308469030992781557966909983211910963157763678726120154469030856807730587971859910379069462105489708001873004723798633342340521799560185957916958401869207109443355859123561156747098129524433371596461424856004227854241384374972430825095073282950873641,
        "rb.cnf" : 538812462750928282721716308734898413194103864553832956073815148020987917365241105816807625188941823391012398976196112157887068449390989186368113511090801920577156367304804512491926215360520651047719401241944845883406098779319363314639309655779343140788696509960447536163784266937815828202118895534508004061478961506257883130142920912100543747226035966976598909666696626176,
        "rc.cnf" : 7711354164502494341321469992586296215389019368675210425477306574721979716330369834365339212526497517240811559116147742518536494403648202692367307372347375124735195570200184068550084511307308544710567552927267754358098889434805018720727997741618028500536976980307056282336258719799038253686515232663203157545908220322265455973957544442057437359833452997837270984970131866759705201073382727090176,
        "sat-grid-pbl-0010.cnf" : 593962746002226256572855,
        "sat-grid-pbl-0015.cnf" : 3012964503482414330783936367006634039953704207876657607,
        "sat-grid-pbl-0020.cnf" : 505529009203800755681036389586231130590495376930061598744188591803177007627164988944437560726719,
        "sat-grid-pbl-0025.cnf" : 18051755963842913831469740173876709710679368179331841451265906666609267268260641283484985491991777996404498069889539875259775880457350865091313349862141,
        "sat-grid-pbl-0030.cnf" : 154089430409193243326541334620745040441316292191902104056995281410078040886572580007404908279631836580121336881940079952338956656599622521005239507469477568002534440349139077306892061020210022834318422387583588123648727,
        "c432.isc" : 68719476736,
        "c499.isc" : 2199023255552,
        "c880.isc" : 1152921504606846976,
        "c1355.isc" : 2199023255552,
        "c1908.isc" : 8589934592,
        "c2670.isc" : 13803492693581127574869511724554050904902217944340773110325048447598592,
#        "c7552.isc" : 205688069665150755269371147819668813122841983204197482918576128
    }

    for model_txt_file, Z in model_counts.iteritems():
        print "starting on model:", track_progress, " name:", model_txt_file
        sat_problem = read_SAT_problem("SAT_problems_cnf/%s" % model_txt_file)
        estimators = estimate_sharp_sat(sat_problem, gumbel_trials=NUM_GUMBEL_UPPER_BOUND_PERTURBATIONS, k=OUR_K)

        exact_log_Zs.append(Decimal(Z).ln())
        our_estimators.append(estimators['barv_estimator'])
        gumbel_upper_bounds.append(estimators['gumbel_upper'])


        #print "log(z) =", Decimal(Z).ln()#np.log(Z_sgp0010)
        #print "our estimator =", estimators['barv_estimator']
        #print "gumbel_upper bound =", estimators['gumbel_upper']

        track_progress+=1


    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(range(len(model_counts)), exact_log_Zs, 'go', label='exact log(Z)')
    ax.plot(range(len(model_counts)), our_estimators, 'r*', label='our estimator, k=%d' %OUR_K)
    ax.plot(range(len(model_counts)), gumbel_upper_bounds, 'bs', label='gumbel upper bound num_perturb=%d' % NUM_GUMBEL_UPPER_BOUND_PERTURBATIONS)


    plt.title('#SAT Model Count upper Bounds and Estimates')
    plt.xlabel('model index, not meaningful')
    plt.ylabel('log(Z)')
    lgd = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    fig.savefig('cur_test', bbox_extra_artists=(lgd,), bbox_inches='tight')    
    plt.close()   
