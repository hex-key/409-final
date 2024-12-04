#!/usr/bin/env python3
"""
Non-neural baseline system for the SIGMORPHON 2020 Shared Task 0.
Author: Mans Hulden
Modified by: Tiago Pimentel
Modified by: Jordan Kodner
Modified by: Omer Goldman
Last Update: 22/03/2021
"""

import sys, os, getopt, re
from functools import wraps
from glob import glob

"""
Finds the Hamming distance of two strings s and t (number of places where the two characters differ)
unts default starting from word initial and disregards any trailing characters that don't have corresponding place in the other word (so max possible value is the len of shorter word)
ex. "alone" and "lone" has a distance of 4 /// "app" and "appendix" 0 ///
"""
def hamming(s,t):
    return sum(1 for x,y in zip(s,t) if x != y)


"""
Code tries different alignment positions and returns the best one

Example: s = BCD, t = ABCDE, it will try
(lpad part)
BCD_____    _BCD____    __BCD___    etc.
___ABCDE    ___ABCDE    ___ABCDE

then (rpad part)
_____BCD    _____BCD    _____BCD    etc.
___ABCDE    __ABCDE_    _ABCDE__  

then at the end it will (i think) end up with 
bs = ____BCD_
bt = ___ABCDE

then s_pad and t_pad are trimmed versions (itll remove the trailing whitespace at the beginning)
so output = _BCD_, ABCDE (as expected)

used in the code as halign(lemma, form)
"""
def halign(s,t):
    """Align two strings (s and t) by Hamming distance."""
    slen, tlen = len(s), len(t)
    bs = bt = "" # best s and t
    minscore = slen + tlen + 1 # start with the worst possible score

    # try increasing amounts of left-padding for s
    for lpad in range(0, tlen+1):
        s_padded = '_' * lpad + s + (tlen - lpad) * '_'
        t_padded = slen * '_' + t
        score = hamming(s_padded, t_padded)
        if score < minscore:
            bs, bt = s_padded, t_padded
            minscore = score

    # try increasing amounts of right-padding for t
    for rpad in range(0, slen+1):
        s_padded = tlen * '_' + s
        t_padded = (slen - rpad) * '_' + t + '_' * rpad
        score = hamming(s_padded, t_padded)
        if score < minscore:
            bs, bt = s_padded, t_padded
            minscore = score

    # trim whitespaces if they're whitespace in both bs and bt
    zipped = list(zip(bs,bt))
    s_pad  = ''.join(s_ch for s_ch, t_ch in zipped if s_ch != '_' or t_ch != '_')
    t_pad = ''.join(t_ch for s_ch, t_ch in zipped if s_ch != '_' or t_ch != '_')
    return s_pad, t_pad

# understanding newin & newout syntax:
    # str.join() method concatenates elements of an iterable---in this case, a list---into a single string
    # ''.join(...) means that the characters generated will be concatenated without any separator.
    # LOGIC:
        #for each (i, o) in zipped:  
            #if i is not equal to '_' OR o is not equal to '_':  
                #newin = newin + i  
                #newout = newout + o 


def levenshtein(s, t, inscost = 1.0, delcost = 1.0, substcost = 1.0): #takes strings s and t, and initializes manipulation costs
    """ Calculates the levenshtein distance between two strings, returning the aligned representations with underscores. 
      The levenshtein distance considers insertions, deletions, and substitutions. 
      
      An inner function, lrec, tracks the alignment using variables for past and remaining characters, as well as the total cost of all manipulations. 
      
      The function returns the aligned version of both strings as well as the total edit distance. """
    
    @memolrec
    
    #internal function, takes params:
        #spast & tpast = aligned portions of the strings at the current point
        #srem & trem = the remaining unaligned portions of the strings at the current point

    def lrec(spast, tpast, srem, trem, cost): 
       
        # if the remaining part of s is empty,
            # align the remaining part of t by padding spast with underscores and adding the remaining characters of t to tpast.
            # and add the length of trem to cost.
        if len(srem) == 0:
            return spast + len(trem) * '_', tpast + trem, '', '', cost + len(trem)
        
        # if the remaining part of t is empty,
            # align the remaining part of s by padding tpast with underscores and adding the remaining characters of s to spast.
            # and add the length of srem to cost.
        if len(trem) == 0:
            return spast + srem, tpast + len(srem) * '_', '', '', cost + len(srem)

        addcost = 0
        # if the current remaining first characters of s and t are different,
            # there is a cost to substitute them
        if srem[0] != trem[0]:
            addcost = substcost

        # the function considers all possible edits of substitutions, insertions, and deletions. It compares their costs and returns the tuple with the smallest cost value.
        return min((lrec(spast + srem[0], tpast + trem[0], srem[1:], trem[1:], cost + addcost),
                   lrec(spast + '_', tpast + trem[0], srem, trem[1:], cost + inscost),
                   lrec(spast + srem[0], tpast + '_', srem[1:], trem, cost + delcost)),
                   key = lambda x: x[4])
    
    # the main function calls lrec with empty strings as spast and t as tpast (they haven't been traversed yet)
    # s as srem, t as trem
    # and initial cost of 0
    answer = lrec('', '', s, t, 0)
    
    # from the minimum cost tuple returned by lrec, it returns the 0th, 1st, and 4th index of the tuple. 
    # or, the aligned version of s, aligned version of t, and the total cost.n
    return answer[0],answer[1],answer[4]

def memolrec(func):
    """Memoizer for Levenshtein. This is used to cache results of recursive subproblems from the levenshtein method. 
    If a pair of remaining s and t strings is not in the cache already, call the levenshtein function and store the result
    after trimming unnecessary parts of sp and tp."""
    cache = {}
    @wraps(func)
    # sp & tp = already aligned string portions
    # sr & tr = remaining unaligned string portions
    # cost = total cost so far
    def wrap(sp, tp, sr, tr, cost):
        if (sr,tr) not in cache:
            res = func(sp, tp, sr, tr, cost)
            cache[(sr,tr)] = (res[0][len(sp):], res[1][len(tp):], res[4] - cost)
        return sp + cache[(sr,tr)][0], tp + cache[(sr,tr)][1], '', '', cost + cache[(sr,tr)][2]
    return wrap


"""
Used in prefix_suffix_rules_get
"""
def alignprs(lemma, form):
    """
    Break lemma/form into three parts:
    based on leading and trailing '_', doesn't take into account different characters
    ex. ('', 'demonstrate', '__', '', 'demonstrati', 'on')
    IN:  1 | 2 | 3
    OUT: 4 | 5 | 6
    1/4 are assumed to be prefixes, 2/5 the stem, and 3/6 a suffix.
    1/4 and 3/6 may be empty.
    """

    al = levenshtein(lemma, form, substcost = 1.1) # Force preference of 0:x or x:0 by 1.1 cost
    alemma, aform = al[0], al[1]
    # leading spaces
    lspace = max(len(alemma) - len(alemma.lstrip('_')), len(aform) - len(aform.lstrip('_')))
    # trailing spaces
    tspace = max(len(alemma[::-1]) - len(alemma[::-1].lstrip('_')), len(aform[::-1]) - len(aform[::-1].lstrip('_')))
    return alemma[0:lspace], alemma[lspace:len(alemma)-tspace], alemma[len(alemma)-tspace:], aform[0:lspace], aform[lspace:len(alemma)-tspace], aform[len(alemma)-tspace:]

"""
gets input from alignprs which means padded based on levenshtein distance and ignoring substitutions
    ex. ('', 'demonstrate', '__', '', 'demonstrati', 'on')

example for ^: 

lem_rs =  demonstrate__>
form_rs = demonstration>

s_rules = {(demonstrate__>, demonstration>)
           (emonstrate__>, emonstration>)
           (monstrate__>, monstration>)
           etc. }
and then it trims all the underscores 

"""
def prefix_suffix_rules_get(lemma, form):
    """Extract a number of suffix-change and prefix-change rules based on a given example lemma + inflected form."""
    lp,lr,ls,fp,fr,fs = alignprs(lemma, form)  # lemma prefix/root/suffix, form prefix/root/suffix

    def make_rules(l, f):
        rules = set()
        for i in range(min(len(l), len(f))):
            rules.add((l[i:], f[i:]))
        rules = {(x[0].replace('_',''), x[1].replace('_','')) for x in rules}
        return rules

    # Suffix rules (using root+suffix) 
    lem_rs  = lr + ls + ">"
    form_rs = fr + fs + ">"
    srules = make_rules(lem_rs, form_rs)

    # Prefix rules
    prules = set()
    if len(lp) + len(fp) >= 0: # if there is any prefix rule
        lem_p = "<" + lp
        form_p = "<" + fp
        prules = make_rules(lem_p, form_p)

    return prules, srules

"""
allprules and allsrules are in this form
{
    "msd" : {
        ("rule-in", "rule-out"): frequency
    },

    [etc etc for different msd-s]
}
"""
def apply_best_rule(lemma, msd, allprules, allsrules):
    """
    Applies the longest-matching suffix-changing rule given an input
    form and the MSD. Length ties in suffix rules are broken by frequency.
    For prefix-changing rules, only the most frequent rule is chosen.
    """

    bestrulelen = 0
    base = "<" + lemma + ">"
    if msd not in allprules and msd not in allsrules:
        return lemma # Haven't seen this inflection, so bail out

    # the "best rule" is found by max(applicablerules) where the thing being compared is the max length from rule[2], and the length of x[0] and x[1]. thats kind a wild ngl 

    if msd in allsrules:
        applicablerules = [(x[0],x[1],y) for x,y in allsrules[msd].items() if x[0] in base]
        
        if applicablerules:
            bestrule = max(applicablerules, key = lambda x: (len(x[0]), x[2], len(x[1])))
            base = base.replace(bestrule[0], bestrule[1])
            # CHANGE
            if "IMP" in msd:
                print("bestrule " + bestrule)

    if msd in allprules:
        applicablerules = [(x[0],x[1],y) for x,y in allprules[msd].items() if x[0] in base]
        

        if applicablerules:
            bestrule = max(applicablerules, key = lambda x: (x[2]))
            base = base.replace(bestrule[0], bestrule[1])

    base = base.replace('<', '')
    base = base.replace('>', '')
    return base


def numleadingsyms(s, symbol):
    return len(s) - len(s.lstrip(symbol))


def numtrailingsyms(s, symbol):
    return len(s) - len(s.rstrip(symbol))

###############################################################################


def main(argv):
    options, _ = getopt.gnu_getopt(argv[1:], 'ohp:', ['output','help','path='])
    TEST, OUTPUT, HELP, path = False, False, False, './data/'
    for opt, arg in options:
        if opt in ('-o', '--output'):
            OUTPUT = True
        if opt in ('-t', '--test'):
            TEST = True
        if opt in ('-h', '--help'):
            HELP = True
        if opt in ('-p', '--path'):
            path = arg

    if HELP:
            print("\n*** Baseline for the SIGMORPHON 2020 shared task ***\n")
            print("By default, the program runs all languages only evaluating accuracy.")
            print("To create output files, use -o")
            print("The training and dev-data are assumed to live in ./part1/development_languages/")
            print("Options:")
            print(" -o         create output files with guesses (and don't just evaluate)")
            print(" -t         evaluate on test instead of dev")
            print(" -p [path]  data files path. Default is ../data/")
            quit()

    totalavg, numlang = 0.0, 0
<<<<<<< HEAD
    for lang in [os.path.splitext(d)[0] for d in os.listdir(path) if '.trn' in d]:
        allprules, allsrules = {}, {}

        # RULES dictionaries are of the format
        # {'set of features': 'base suffix' 'modified suffix' '}
        if not os.path.isfile(path + lang +  ".trn"):
            continue
        lines = [line.strip() for line in open(path + lang + ".trn", "r", encoding='utf8') if line != '\n']
=======
    lang = "spa"
    
    allprules, allsrules = {}, {}
>>>>>>> master

    lines = [line.strip() for line in open(path + lang + ".trn", "r", encoding='utf8') if line != '\n']

    # First, test if language is predominantly suffixing or prefixing
    # If prefixing, work with reversed strings
    prefbias, suffbias = 0,0
    for l in lines:
        lemma, _, form = l.split(u'\t')
        aligned = halign(lemma, form)
        if ' ' not in aligned[0] and ' ' not in aligned[1] and '-' not in aligned[0] and '-' not in aligned[1]:
            prefbias += numleadingsyms(aligned[0],'_') + numleadingsyms(aligned[1],'_')
            suffbias += numtrailingsyms(aligned[0],'_') + numtrailingsyms(aligned[1],'_')

    # Read in training data lines and extract transformation rules from pairs
    for l in lines: 
        lemma, msd, form = l.split(u'\t')
        if prefbias > suffbias:
            lemma = lemma[::-1]
            form = form[::-1]
        prules, srules = prefix_suffix_rules_get(lemma, form)

        if msd not in allprules and len(prules) > 0:
            allprules[msd] = {}
        if msd not in allsrules and len(srules) > 0:
            allsrules[msd] = {}

        for r in prules:
            if (r[0],r[1]) in allprules[msd]:
                allprules[msd][(r[0],r[1])] += 1
            else:
                allprules[msd][(r[0],r[1])] = 1

        for r in srules:
            if (r[0],r[1]) in allsrules[msd]:
                allsrules[msd][(r[0],r[1])] += 1
            else:
                allsrules[msd][(r[0],r[1])] = 1

    # Read in trial data lines from either spa.tst or spa.dev
    if TEST:
        devlines = [line.strip() for line in open(path + lang + ".tst", "r", encoding='utf8') if line != '\n']
    else:
        devlines = [line.strip() for line in open(path + lang + ".dev", "r", encoding='utf8') if line != '\n']
    
    # Run eval on relevant set 
    numcorrect = 0
    numguesses = 0
    if OUTPUT:
        outfile = open(path + lang + ".out", "w", encoding='utf8')
    for l in devlines:
        lemma, msd, correct = l.split(u'\t')
        if prefbias > suffbias:
            lemma = lemma[::-1]
        outform = apply_best_rule(lemma, msd, allprules, allsrules)
        if prefbias > suffbias:
            outform = outform[::-1]
            lemma = lemma[::-1]
        if outform == correct:
            numcorrect += 1
        numguesses += 1
        if OUTPUT:
            outfile.write(lemma + "\t" + msd + "\t" + outform + "\n")

    if OUTPUT:
        outfile.close()

    totalavg += numcorrect/float(numguesses)

    print(lang + ": " + str(str(numcorrect/float(numguesses)))[0:7])



if __name__ == "__main__":
    main(sys.argv)
