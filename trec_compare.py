#! /usr/bin/python

import os
import sys
import re
sys.path.append ('/home/bemike/lustre/_PR_tools/processors/')
import random,math
from optparse import OptionParser

_N_PERMUTE = 1000
_EPSILON = 5

def int2bin(i, n):
    """Convert decimal integer i to n-bit binary number (string).

    >>> int2bin(0, 8)
    '00000000'

    >>> int2bin(123, 8)
    '01111011'

    >>> int2bin(123L, 8)
    '01111011'

    >>> int2bin(15, 2)
    Traceback (most recent call last):
    ValueError: Value too large for given number of bits.

    """
    hex2bin = {'0': '0000', '1': '0001', '2': '0010', '3': '0011',
               '4': '0100', '5': '0101', '6': '0110', '7': '0111',
               '8': '1000', '9': '1001', 'a': '1010', 'b': '1011',
               'c': '1100', 'd': '1101', 'e': '1110', 'f': '1111'}
    # Convert to hex then map each hex digit to binary equivalent.
    result = ''.join([hex2bin[x] for x in hex(i).lower().replace('l','')[2:]])
                      
    # Shrink result to appropriate length.
    # Raise an error if the value is changed by the truncation.
    if '1' in result[:-n]:
        raise ValueError("Value too large for given number of bits.")
    result = result[-n:]
    # Zero-pad if length longer than mapped result.
    result = '0'*(n-len(result)) + result
    return result
    
def random_vec(bits, max_value=None):
    """Generate a random binary vector of length bits and given max value."""

    vector = ""
    for _ in range(int(bits / 10) + 1):
        i = int((2**10) * random.random())
        vector += int2bin(i, 10)

    if max_value and (max_value < 2 ** bits - 1):
        vector = int2bin((int(vector, 2) / (2 ** bits - 1)) * max_value, bits)
    
    return vector[0:bits]

def mean(seq):
    return float(sum(seq))/len(seq)
    
def randomize(bd,cd):
    methA = []
    methB = []    
    for q in bd:
        methA.append(bd[q])
        methB.append(cd[q])
    ntrials = len(methA)
    true_diff = abs(mean(methA) - mean(methB))
    p_val = 0.0
    for k in range(options.permutations):
        assigments = []
        permutedA = []
        permutedB = []        
        assigments = random_vec(ntrials)
        for i in range(ntrials):
            if assigments[i] == '0':
                permutedA.append( methA[i] )
                permutedB.append( methB[i] )
            else:
                permutedA.append( methB[i] )
                permutedB.append( methA[i] )
        p_val += float( abs( mean( permutedA ) - mean( permutedB ) ) >= true_diff ) 
    return p_val/options.permutations
    
def getRange(r):
    if r == "":
        return (-1,-1)
    (s,e) = r.split('-')    
    return (int(s),int(e))
    
def readTrecEval(f):
    qrange = getRange(options.range)
    qryDict = dict()
    with open(f) as teval:
        for l in teval:
            (metric,qid,value) = l.strip().split()
            if qid != 'all':
                if qrange[0] < 0:
                    if metric == options.metric:
                        qryDict[ qid ] = float(value)
                else:
                    if metric == options.metric and int(qid) in range(*qrange):
                        qryDict[ qid ] = float(value)
    return qryDict

def computeMean(qd):
    return sum(qd.values())/len(qd.values())

def computePrctImprv(bd,cd):
    bMean = computeMean( bd )
    cMean = computeMean( cd )
    return (cMean/bMean - 1 )*100

def computeIH(bd,cd):
    nImpr = 0.0
    nHurt = 0.0
    total = len( bd )
    for q in bd:
        if bd[q] == 0:
            if cd[q] == 0:
                continue
            else:
                nImpr += 1
        elif abs( cd[q]/bd[q] - 1 )*100 > options.threshold:
            if bd[q] > cd[q]:
                nHurt += 1
            elif bd[q] < cd[q]:
                nImpr += 1
    return (100*nImpr/total,100*nHurt/total)
    
def detailedImprv(bd,cd):
    details = []
    for q in bd:
        details.append( (q,cd[q] - bd[q]) )
    details.sort( key=lambda x: x[1], reverse=True )
    return details

    
def robustness(bd,cd):
    bins = ( range(0,25), range(25,50), range(50,75), range(75,100) ) 
    pos = [0]*5
    neg = [0]*5
    for q in bd:
        if bd[q] == 0:
            if cd[q] == 0:
                pos[0] += 1
            else:
                pos[4] += 1
        else:
            delta = int(( cd[q]/bd[q] - 1 )*100)
            for i,b in enumerate(bins):
                if abs( delta ) in b:
                    if delta >= 0:
                        pos[i] += 1
                        break
                    else:
                        neg[i] += 1
                        break
            else:
                if delta > 0:
                    pos[4] += 1
                else:
                    neg[4] += 1
    return (pos,neg)
                        

    
if __name__ == '__main__':
    usage = "usage: %prog [options] <baseline_run> <candidate_run>"
    parser = OptionParser(usage=usage)    
    parser.add_option("-m", "--metric", type="string", dest="metric", default="map", help="Metric to use for comparison")
    parser.add_option("-r", "--range", type="string", dest="range", default="", help="Range of queries to compare. Input as: start-end")
    parser.add_option("-p", "--permutations", type="int", dest="permutations", default= _N_PERMUTE, help="Number of permutations")
    parser.add_option("-t", "--threshold", type="int", dest="threshold", default= _EPSILON, help="Improvement threshold")
    (options, args) = parser.parse_args()
    if len(args) < 2:
        parser.print_help()
        sys.exit(0)
    else:
        baseDict = readTrecEval( args[0] )
        candDict = readTrecEval( args[1] )
        # Output results
        if len( baseDict.keys() ) == 0 or len( candDict.keys() ) == 0:
            print '** Error: Invalid trec_eval file as input!'
            sys.exit(0)
        if len(baseDict.keys()) != len(candDict.keys()):
            print '** Error: Base and candidate runs have different number of queries!'
            sys.exit(0)
        else:
            print '*'*50
            print 'Using metric *%s*' %(options.metric)
            print '*'*50
            print 'Number of queries\t%s' %(len(baseDict.keys()))
            print 'Mean baseline metric\t%.4f' %(computeMean(baseDict))
            print 'Mean candidate metric\t%.4f' %(computeMean(candDict))
            print '%% improvement:\t\t%.1f' %(computePrctImprv(baseDict, candDict))
            ih = computeIH(baseDict, candDict)
            print '%% queries improved\t%.1f (by more than %s%%)' %(ih[0], options.threshold)
            print '%% queries hurt\t\t%.1f (by more than %s%%)' %(ih[1], options.threshold)
            print 'P-val (%s iter)\t%.3f' %(options.permutations,randomize(baseDict,candDict))            
            print '*'*50
            rb = robustness(baseDict,candDict)
            print 'Robustness\t\t+[0,25):%s, +[25,50):%s, +[50,75):%s, +[75,100):%s, +[100,):%s' %(rb[0][0],rb[0][1],rb[0][2],rb[0][3],rb[0][4])
            print '\t\t\t-(0,25):%s, -[25,50):%s, -[50,75):%s, -[75,100):%s, -[100,):%s' %(rb[1][0],rb[1][1],rb[1][2],rb[1][3],rb[1][4])
            print '*'*50
            print 'Query identifiers\t(%s)' %(' '.join(sorted(baseDict.keys())))
            print 'Per query %s changes\t(%s)' %(options.metric, ' '.join(['%s:%.2f' %(x) for x in detailedImprv(baseDict,candDict)]) )
            print '*'*50
            
            
            