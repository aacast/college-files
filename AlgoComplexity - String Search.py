# Python program for KMP and BMAlgorithm from GeeksforGeeks, Brute Force Algorithm from https://johnlekberg.com/
from time import process_time
from pathlib import Path
from collections import defaultdict
import sys
import string

def bruteForce(pat, txt):
    n, m = len(txt), len(pat)
    for i in range(1 + (n - m)):
        match = True
        for j in range(m):
            if txt[i + j] != pat[j]:
                match = False
                break
        if match:
            print("Found pattern at index " + str(i))
            break
        
def KMPSearch(pat, txt):
    M = len(pat)
    N = len(txt)
  
    # create lps[] that will hold the longest prefix suffix 
    # values for pattern
    lps = [0]*M
    j = 0 # index for pat[]
  
    # Preprocess the pattern (calculate lps[] array)
    computeLPSArray(pat, M, lps)
  
    i = 0 # index for txt[]
    while i < N:
        if pat[j] == txt[i]:
            i += 1
            j += 1
  
        if j == M:
            print ("Found pattern at index " + str(i-j))
            j = lps[j-1]
  
        # mismatch after j matches
        elif i < N and pat[j] != txt[i]:
            # Do not match lps[0..lps[j-1]] characters,
            # they will match anyway
            if j != 0:
                j = lps[j-1]
            else:
                i += 1
  
def computeLPSArray(pat, M, lps):
    len = 0 # length of the previous longest prefix suffix
  
    lps[0] # lps[0] is always 0
    i = 1
  
    # the loop calculates lps[i] for i = 1 to M-1
    while i < M:
        if pat[i]== pat[len]:
            len += 1
            lps[i] = len
            i += 1
        else:
            # This is tricky. Consider the example.
            # AAACAAAA and i = 7. The idea is similar 
            # to search step.
            if len != 0:
                len = lps[len-1]
  
                # Also, note that we do not increment i here
            else:
                lps[i] = 0
                i += 1
                
                
                
def preprocess_strong_suffix(shift, bpos, pat, m):
  
    # m is the length of pattern
    i = m
    j = m + 1
    bpos[i] = j
  
    while i > 0:
          
        '''if character at position i-1 is 
        not equivalent to character at j-1, 
        then continue searching to right 
        of the pattern for border '''
        while j <= m and pat[i - 1] != pat[j - 1]:
              
            ''' the character preceding the occurrence 
            of t in pattern P is different than the 
            mismatching character in P, we stop skipping
            the occurrences and shift the pattern 
            from i to j '''
            if shift[j] == 0:
                shift[j] = j - i
  
            # Update the position of next border
            j = bpos[j]
              
        ''' p[i-1] matched with p[j-1], border is found. 
        store the beginning position of border '''
        i -= 1
        j -= 1
        bpos[i] = j
  
# Preprocessing for case 2
def preprocess_case2(shift, bpos, pat, m):
    j = bpos[0]
    for i in range(m + 1):
          
        ''' set the border position of the first character 
        of the pattern to all indices in array shift
        having shift[i] = 0 '''
        if shift[i] == 0:
            shift[i] = j
              
        ''' suffix becomes shorter than bpos[0], 
        use the position of next widest border
        as value of j '''
        if i == j:
            j = bpos[j]
  

def search(text, pat):
  
    s = 0
    m = len(pat)
    n = len(text)
  
    bpos = [0] * (m + 1)
  
    shift = [0] * (m + 1)
 
    preprocess_strong_suffix(shift, bpos, pat, m)
    preprocess_case2(shift, bpos, pat, m)
  
    while s <= n - m:
        j = m - 1
        while j >= 0 and pat[j] == text[s + j]:
            j -= 1
        if j < 0:
            print("Found pattern at index = %d" % s)
            s += shift[0]
        else:
            s += shift[j + 1]
 
def BoyerMooreHorspool(pat, txt):
    m = len(pat)
    n = len(txt)
    if m > n: return -1
    skip = []
    for k in range(256): skip.append(m)
    for k in range(m - 1): skip[ord(pat[k])] = m - k - 1
    skip = tuple(skip)
    k = m - 1
    while k < n:
        j = m - 1; i = k
        while j >= 0 and txt[i] == pat[j]:
            j -= 1; i -= 1
        if j == -1:  print("Found pattern at index " + str(i + 1)) 
        k += skip[ord(txt[k])]
    return -1
    
def boyer_moore_horspool(pat, txt):
    m = len(pat)
    n = len(txt)

    if m > n:
        return -1

    skip = defaultdict(lambda: m)
    found_indexes = []

    for k in range(m - 1):
        skip[ord(pat[k])] = m - k - 1

    k = m - 1

    while k < n:
        j = m - 1
        i = k
        while j >= 0 and txt[i] == pat[j]:
            j -= 1
            i -= 1
        if j == -1:
            print("Found pattern at index " + str(i + 1)) 

        k += skip[ord(txt[k])]

    return found_indexes
    
    
    
    
d = 256
  
# pat  -> pattern
# txt  -> text
# q    -> A prime number
  
def KarpRabin(pat, txt, q):
    M = len(pat)
    N = len(txt)
    i = 0
    j = 0
    p = 0    # hash value for pattern
    t = 0    # hash value for txt
    h = 1
  
    # The value of h would be "pow(d, M-1)% q"
    for i in range(M-1):
        h = (h * d)% q
  
    # Calculate the hash value of pattern and first window
    # of text
    for i in range(M):
        p = (d * p + ord(pat[i]))% q
        t = (d * t + ord(txt[i]))% q
  
    # Slide the pattern over text one by one
    for i in range(N-M + 1):
        # Check the hash values of current window of text and
        # pattern if the hash values match then only check
        # for characters on by one
        if p == t:
            # Check for characters one by one
            for j in range(M):
                if txt[i + j] != pat[j]:
                    break
  
            j+= 1
            # if p == t and pat[0...M-1] = txt[i, i + 1, ...i + M-1]
            if j == M:
                print ("Pattern found at index " + str(i))
  
        # Calculate hash value for next window of text: Remove
        # leading digit, add trailing digit
        if i < N-M:
            t = (d*(t-ord(txt[i])*h) + ord(txt[i + M]))% q
  
            # We might get negative values of t, converting it to
            # positive
            if t < 0:
                t = t + q
  
    
text_file = open("string.txt", "r")
 
data = text_file.read()
 
text_file.close()
txt = data
pat = "ababbababab"


print("Brute Algorithm") 
t2_start = process_time()
bruteForce(pat, txt)
t2_stop = process_time()
print("Elapsed time:", t2_stop, t2_start) 
print("Elapsed time during the whole program in seconds:",
                            t2_stop-t2_start) 
print("")

print("KMP Algorithm")
t1_start = process_time()
KMPSearch(pat, txt)
t1_stop = process_time()
print("Elapsed time:", t1_stop, t1_start)
print("Elapsed time during the whole program in seconds:",
                            t1_stop-t1_start)

print("")
 
print("Boyer Moore Algorithm")
t3_start = process_time()
if __name__ == '__main__':
    search(txt, pat)
t3_stop = process_time()
print("Elapsed time:", t3_stop, t3_start)
print("Elapsed time during the whole program in seconds:",
                            t3_stop-t3_start)
                            
                            
print("")
 
print("Boyer Moore Horspool Algorithm")
t4_start = process_time()
if __name__ == '__main__':
    s = BoyerMooreHorspool(pat, txt)
t4_stop = process_time()
print("Elapsed time:", t4_stop, t4_start)
print("Elapsed time during the whole program in seconds:",
                            t4_stop-t4_start)
                            

print("")
 
print("Karp Rabin Algorithm")
q = 101
t5_start = process_time()
KarpRabin(pat, txt, q)
t5_stop = process_time()
print("Elapsed time:", t5_stop, t5_start)
print("Elapsed time during the whole program in seconds:",
                            t5_stop-t5_start)
    

