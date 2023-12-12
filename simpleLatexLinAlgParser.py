import re
import numpy as np
import math

UNICODE_CORRECTION_TRANSLATER = str.maketrans(
    {
        '−' : '-'
    }
)


def mathjaxToPyMatrix1(mathjax):
    matrix = []
    
    for raw_row in re.split(r"\\",mathjax):
        matrix.append(list(
            map(lambda s: float(s.replace("−","-").strip()),
                re.findall(r"[-\d.−]+",raw_row))
            )
        )
    return list(filter(None,matrix))

def pyToMatlabMatrix(matrix):
    return "[" + ' ; '.join(map(str,matrix)).replace("[","").replace("]","") + "]"

def texToPyCoeffAndBMatrices1(tex_sysofeq):
    A = []
    b = []
    row_ix = 0
    tex_sysofeq = "\n".join(tuple(filter(lambda srow: bool(re.search(r'\d', srow)),tex_sysofeq.split("\n"))))

    for raw_row in re.split(r"\\",tex_sysofeq.replace("- ","-").replace("+","")):
        parsed_row = re.split(r"\s&(?!\s*=)",raw_row)
        _bi = parsed_row.pop().strip()
        bi = re.match(r"[-\d.]+",_bi)
        if not bi:
            continue
        b.append([
            float(bi.group(0))
        ])
        A.append([])
        for raw_el in parsed_row:
            raw_el = raw_el.strip()
            el = re.match(r"(?<!\w)[-\d.]+",raw_el)
            if el is None:
                el = 0. if not re.match(r"[\w]+",raw_el) else 1.
            else:
                el = -1. if (el:=el.group(0).strip()) == '-' else float(el)
            A[row_ix].append(el)
        row_ix+= 1
    
    return  list(filter(None,A)), b 


if __name__ == '__main__':

     
    nums = np.fromiter(map(float,re.split(r"\n|\s","""2310 2320 2010 10800 2190 3360 5640 2540 3360 11800 2010
3430 10600 7370 2160 3200 2020 2850 3500 10200 8550 9500
2260 7730 2250""")),np.float_)
    sample_median = lambda nums : sorted(nums)[len(nums)//2]
    c_sample_sd = lambda nums: math.sqrt((sum(nums**2) - 
                                          (sum(nums)**2)/len(nums))
                                          /(len(nums) -1))
    
    
    print(len(nums[nums < 5000.])/len(nums))

    exit()

    mathjax = \
    """
    A = \left( \begin{matrix} 
    2 & 1 & 3 \\
    8 & 8 & 13 \\
    12 & 14 & 22
    \end{matrix} 
    \right)
    """
    tex = """
 3 & -2 & 4 \\
 21 & -15 & 33 \\
 -9 & 2 & 15
    """
    matrix = mathjaxToPyMatrix1(tex)
    octaveStr = pyToMatlabMatrix(matrix)
    print(np.matrix(matrix))
    print(octaveStr)
    
    texSystemOfEqs = \
        """
\begin{matrix} 
x_1  & -6x_2   & +7x_3  & -9x_4     &= & -44.1 \\
x_1  & -5x_2   &             &                & = & -6.8 \\
        &    x_2 & -5x_3     &                &= &  12.9 \\
        &          &        x_3  & -5x_4      & = & -13.3 \\
\end{matrix}
        """.translate(UNICODE_CORRECTION_TRANSLATER)
        
    A,b = texToPyCoeffAndBMatrices1(texSystemOfEqs)
    print(np.matrix(A))
    print(np.matrix(b))
    # _n = len(A) // 2
    # for i in range(len(A)):
    #     print("|",end="")
    #     print(*A[i],sep=" ",end= "| ")
    #     if i == _n:
    #         print("=  ",end="|")
    #     else:
    #         print("   ",end="|")
    #     print(*b[i],end="|\n\r")
            
    print(pyToMatlabMatrix(A))
    print(pyToMatlabMatrix(b))


    